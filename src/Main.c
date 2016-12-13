#include "XSbench_header.h"
#include "Types.h"
#include "CL/opencl.h"
#include "AOCL_Utils.h"

using namespace aocl_utils; 

#ifdef MPI
#include<mpi.h>
#endif

//inputs
Inputs in; 
int num_dpoints = 0;

//OPENCL config
cl_platform_id platform = NULL;
cl_device_id device = NULL;
cl_context context = NULL;
cl_command_queue queue = NULL;
cl_kernel kernel = NULL;
cl_program program = NULL;
cl_mem macro_xs_vector_buf;
cl_mem num_nucs_buf, concs_buf, energy_grid_buf, nuclide_grids_buf, mats_buf; 

//static functions
static bool init_ocl();
static int get_num_datapoints(int *);
static void cleanup();

int main( int argc, char* argv[] )
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================
	int version = 13;
	int mype = 0;
	//int max_procs = omp_get_num_procs();
	int i;	// thread, mat;
	//unsigned long seed;
	double omp_start, omp_end, p_energy;
	unsigned long long vhash = 0;
	int nprocs;

	#ifdef MPI
	MPI_Status stat;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	#endif
	
	// rand() is only used in the serial initialization stages.
	// A custom RNG is used in parallel portions.
	#ifdef VERIFICATION
	srand(26);
	#else
	srand(time(NULL));
	#endif

	// Process CLI Fields -- store in "Inputs" structure
	//Inputs in = read_CLI( argc, argv );
	in = read_CLI( argc, argv );

	// Set number of OpenMP Threads
	omp_set_num_threads(in.nthreads); 

	// Print-out of Input Summary
	if( mype == 0 )
		print_inputs( in, nprocs, version );

	// =====================================================================
	// Prepare Nuclide Energy Grids, Unionized Energy Grid, & Material Data
	// =====================================================================

	// Allocate & fill energy grids
	#ifndef BINARY_READ
	if( mype == 0) printf("Generating Nuclide Energy Grids...\n");
	#endif
	
	NuclideGridPoint ** nuclide_grids = gpmatrix(in.n_isotopes,in.n_gridpoints);
	
	#ifdef VERIFICATION
	generate_grids_v( nuclide_grids, in.n_isotopes, in.n_gridpoints );	
	#else
	generate_grids( nuclide_grids, in.n_isotopes, in.n_gridpoints );	
	#endif
	
	// Sort grids by energy
	#ifndef BINARY_READ
	if( mype == 0) printf("Sorting Nuclide Energy Grids...\n");
	sort_nuclide_grids( nuclide_grids, in.n_isotopes, in.n_gridpoints );
	#endif

	// Prepare Unionized Energy Grid Framework
	#ifndef BINARY_READ
	GridPoint * energy_grid = generate_energy_grid( in.n_isotopes,
	                          in.n_gridpoints, nuclide_grids ); 	
	#else
	GridPoint * energy_grid = (GridPoint *)malloc( in.n_isotopes *
	                           in.n_gridpoints * sizeof( GridPoint ) );
	int * index_data = (int *) malloc( in.n_isotopes * in.n_gridpoints
	                   * in.n_isotopes * sizeof(int));
	for( i = 0; i < in.n_isotopes*in.n_gridpoints; i++ )
		energy_grid[i].xs_ptrs = &index_data[i*in.n_isotopes];
	#endif

	// Double Indexing. Filling in energy_grid with pointers to the
	// nuclide_energy_grids.
	#ifndef BINARY_READ
	set_grid_ptrs( energy_grid, nuclide_grids, in.n_isotopes, in.n_gridpoints );
	#endif

	#ifdef BINARY_READ
	if( mype == 0 ) printf("Reading data from \"XS_data.dat\" file...\n");
	binary_read(in.n_isotopes, in.n_gridpoints, nuclide_grids, energy_grid);
	#endif
	
	// Get material data
	if( mype == 0 )
		printf("Loading Mats...\n");
	int *num_nucs  = load_num_nucs(in.n_isotopes);
	num_dpoints = get_num_datapoints(num_nucs);

	int **mats     = load_mats(num_nucs, in.n_isotopes);

	#ifdef VERIFICATION
	double **concs = load_concs_v(num_nucs);
	#else
	double **concs = load_concs(num_nucs);
	#endif

	#ifdef BINARY_DUMP
	if( mype == 0 ) printf("Dumping data to binary file...\n");
	binary_dump(in.n_isotopes, in.n_gridpoints, nuclide_grids, energy_grid);
	if( mype == 0 ) printf("Binary file \"XS_data.dat\" written! Exiting...\n");
	return 0;
	#endif

	// =====================================================================
	// Cross Section (XS) Parallel Lookup Simulation Begins
	// =====================================================================

	// Outer benchmark loop can loop through all possible # of threads
/*
	#ifdef BENCHMARK
	for( int bench_n = 1; bench_n <=omp_get_num_procs(); bench_n++ )
	{
		in.nthreads = bench_n;
		omp_set_num_threads(in.nthreads);
 	#endif
*/
	if( mype == 0 )
	{
		printf("\n");
		border_print();
		center_print("SIMULATION", 79);
		border_print();
	}

//	omp_start = omp_get_wtime();
/**  
	//initialize papi with one thread (master) here
	#ifdef PAPI
	if ( PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT){
		fprintf(stderr, "PAPI library init error!\n");
		exit(1);
	}
	#endif	
**/
	//OpenCL initilization
	if( !init() )
		return -1;
 	
	double macro_xs_vector[5];
	double * xs = (double *) calloc(5, sizeof(double));
	cl_int status;
	cl_event write_events[5];
	cl_event kernel_event;
	cl_event finish_event;
	
	status = clEnqueuWriteBuffer(queue, num_nucs_buf, CL_FALSE, 0, 12*sizeof(int), num_nucs, 0, NULL, &write_events[0]);
	checkError(&status, "Failed to enqueue write buffer.\n");

	status = clEnqueuWriteBuffer(queue, concs_buf, CL_FALSE, 0, num_dpoints*sizeof(double), concs, 0, NULL, &write_events[1]);
        checkError(&status, "Failed to enqueue write buffer.\n");
	
	status = clEnqueuWriteBuffer(queue, energy_grid_buf, CL_FALSE, 0, in.n_isotopes*in.n_gridpoints*sizeof(GridPoint), energy_grids, 0, NULL, &write_events[2]);
        checkError(&status, "Failed to enqueue write buffer.\n");

	status = clEnqueuWriteBuffer(queue, nuclide_grids_buf, CL_FALSE, 0, in.n_isotopes*in.n_gridpoints*sizeof(NuclideGridPoint), nuclide_grids, 0, NULL, &write_events[3]);
        checkError(&status, "Failed to enqueue write buffer.\n");

	status = clEnqueuWriteBuffer(queue, mats_buf, CL_FALSE, 0, num_dpoints*sizeof(int), mats, 0, NULL, &write_events[4]);
        checkError(&status, "Failed to enqueue write buffer.\n");
	
	int arg = 0;	
	#ifdef VERIFICATION
		unsigned int verification = 1;
	#else
		unsigned int verification = 0;
	#endif
	status = clSetKernelArg(kernel, arg++, sizeof(unsigned int), &verification);
	checkError(status, "Failed to set arg verification");
	
	status = clSetKernelArg(kernel, arg++, sizeof(long), &in.n_isotopes);
	checkError(status, "Failed to set arg 0");

	status = clSetKernelArg(kernel, arg++, sizeof(long), &in.n_gridpoints);
	checkError(status, "Failed to set arg 1");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &num_nucs_buf);
	checkError(status, "Failed to set kernel arg 2");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &concs_buf);
	checkError(status, "Failed to set arg 3");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &energy_grid_buf);
	checkError(status, "Faile to set arg 4");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &nuclide_grids_buf);
	checkError(status, "Failed to set arg 5");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &mats_buf);
	checkError(status, "Failed to set arg 6");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &macro_xs_vector_buf);
	checkError(status, "failed to ser arg 7");

	int global_work_size = in.lookups; 
	status = clEnqueueNDRangeKernel(queue, kernel, 1, NULL,
       		 &global_work_size, NULL, 5, write_events, &kernel_event);
    	checkError(status, "Failed to launch kernel");
	
	status = clEnqueueReadBuffer(queue, macro_xs_vector_buf, CL_TRUE, 0, 5*sizeof(double), macro_xs_vector, &kernel_event, &finish_event);
	
	clReleaseEvent(write_event[0]);
    	clReleaseEvent(write_event[1]);
 	clReleaseEvent(write_event[2]);
    	clReleaseEvent(write_event[3]);
	clReleaseEvent(write_event[4]);
	
	memcpy(xs, macro_xs_vector, 5*sizeof(double));

	// Verification hash calculation
	// This method provides a consistent hash accross
	// architectures and compilers.

/*	#ifdef VERIFICATION
		char line[256];
		sprintf(line, "%.5lf %d %.5lf %.5lf %.5lf %.5lf %.5lf",
			p_energy, mat,
			macro_xs_vector[0],
			macro_xs_vector[1],
			macro_xs_vector[2],
			macro_xs_vector[3],
			macro_xs_vector[4]);
		unsigned long long vhash_local = hash(line, 10000);
		vhash += vhash_local;
	#endif	
*/	
	if( mype == 0)	
	{	
		printf("\n" );
		printf("Simulation complete.\n" );
	}

	
//	print_results( in, mype, omp_end-omp_start, nprocs, vhash );
/*
	#ifdef BENCHMARK
	}
	#endif
*/
	cleanup();

	return 0;
}

static bool init_ocl()
{
	cl_int status;
        if(!setCwdToExeDir())
		return false;
	platform = findPlatform("Altera");
	if(platform == NULL)
	{
		printf("ERROR: unable to find Altera platform.\n");
		return false;
	}
	
	scoped_array<cl_device_id> devices;
	cl_uint num_devices;
	devices.reset(getDevices(platform, CL_DEVICE_TYPE_ALL, &num_devices));
	device = devices[0];

	context = clCreateContext(NULL, 1, &device, NULL, NULL, &status);
	checkError(status, "Failed to create context.\n");

	queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
	checkError(status, "Failed to create command queue.\n");

	std::string binary_file = getBoardBinaryFile("CalculateXS", device);
	printf("Using AOCX: %s\n", binary_file.c_str());
	program = createProgramFromBinary(context, binary_file.c_str(), &device, 1);

  	status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
	checkError(status, "Failed to build program.\n");
	
	const char *kernel_name = "CalculateXS";
	kernel = clCreateKernel(program, kernel_name, &status);
	checkError(status, "Failed to create kernel.\n");

	//Give buffer 
	macro_xs_vector_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, 5*sizeof(double), NULL, &status);
	checkError(status, "Failed to create output buffer.\n");

       	num_nucs_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, 12*sizeof(int), NULL, &status);
	checkError(status, "Failed to create num_nucs input buffer.\n");

	concs_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, num_dpoints*sizeof(double), NULL, &status);
	checkError(status, "Failed to create concs input buffer.\n");
	
	energy_grid_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, in.n_isotopes*in.n_gridpoints*sizeof(GridPoint), NULL, &status);
	checkError(status, "Failed to create energy_grid input buffer.\n");

	nuclide_grids_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, in.n_isotopes*in.n_gridpoints*sizeof(NuclideGridPoint), NULL, &status);
	checkError(status, "Failed to create nuclide_grids input buffer.\n");

	mats_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, num_dpoints*sizeof(int), NULL, &status);
	checkError(status, "Failed to create mats input buffer.\n");

	return true;
}

static void cleanup()
{
	if(kernel)
                clReleaseKernel(kernel);
        if(program)
                clReleaseProgram(program);
        if(queue)
                clReleaseCommandQueue(queue);
        if(context)
                clReleaseContext(context);
}

static int get_num_datapoints(int *num_nucs)
{
	int num_datapoints = 0;
	for(int i = 0; i < 12; i++)
		num_datapoints += num_nucs[i];

	return num_datapoints;
}

