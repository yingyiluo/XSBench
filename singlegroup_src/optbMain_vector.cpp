#include "Types.h"
#include "XSbench_header.h"

#include "CL/opencl.h"
#include "../common/inc/AOCL_Utils.h"

#include "etrace/rapl_reader.h"

#if defined(USE_INTELOPENCL)
#include <stdio.h>
#include <stdlib.h>
#endif

using namespace aocl_utils; 

#ifdef MPI
#include<mpi.h>
#endif

//inputs
Inputs in; 
int num_dpoints = 0;
int n_iso_grid;
CacheData* firstN = NULL;
int N = 1023;
int STAGE = 10; 

//OPENCL config
unsigned int perthread_lookups;
cl_platform_id platform = NULL;
cl_device_id device = NULL;
cl_context context = NULL;
cl_command_queue queue = NULL;
cl_kernel kernel = NULL;
cl_program program = NULL;
cl_mem vhash_buf;
cl_mem concs_buf, energy_grid_buf, energy_grid_xs_buf, nuclide_grids_buf, mats_buf, firstN_buf; 
//time measurements
double start_ts_sec;
double ocl_start_ts_sec;
double ocl_kernel_ts_sec;
double ocl_post_ts_sec;
double ocl_end_ts_sec;
double end_ts_sec;
cl_ulong profiled_kernel_time_ns;

//functions
static bool init_ocl();
static int get_num_datapoints(cl_int16);
void cleanup();
extern double gettimesec();
double bw();
long grid_search(long, double, GridPoint *);
void calculate_micro_xs( double, int, long, long, GridPoint*, NuclideGridPoint **, int, double * );	

void snap_energy();
void report_energy();

int main( int argc, char* argv[] )
{
	// =====================================================================
	// Initialization & Command Line Read-In
	// =====================================================================
	start_ts_sec = gettimesec();
	
	int version = 13;
	int mype = 0;
	//int max_procs = omp_get_num_procs();
	int i;	// thread, mat;
	//unsigned long seed;
	double omp_start, omp_end, p_energy;
	//unsigned long long vhash = 0;
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
 	n_iso_grid = in.n_isotopes*in.n_gridpoints;
	//firstN.reset(N);
	posix_memalign((void**)&firstN, AOCL_ALIGN, N*sizeof(CacheData));

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
	
	generate_grids( nuclide_grids, in.n_isotopes, in.n_gridpoints );	
	
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
	GridPoint * energy_grid = (GridPoint *)malloc( n_iso_grid * sizeof( GridPoint ) );
	int * index_data = (int *) malloc( n_iso_grid
	                   * in.n_isotopes * sizeof(int));
	for( i = 0; i < n_iso_grid; i++ )
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
	int* num_nucs_p  = load_num_nucs(in.n_isotopes);
	cl_int16 num_nucs;
	if( in.n_isotopes == 68 )
                num_nucs.s0  = 34; // HM Small is 34, H-M Large is 321
        else
                num_nucs.s0  = 321; // HM Small is 34, H-M Large is 321

        num_nucs.s1  = 5;
        num_nucs.s2  = 4;
        num_nucs.s3  = 4;
        num_nucs.s4  = 27;
        num_nucs.s5  = 21;
        num_nucs.s6  = 21;
        num_nucs.s7  = 21;
        num_nucs.s8  = 21;
        num_nucs.s9  = 21;
        num_nucs.sA = 9;
        num_nucs.sB = 9;

	num_dpoints = get_num_datapoints(num_nucs);

	int **mats     = load_mats(num_nucs_p, in.n_isotopes);
	//for(i = 0; i < num_dpoints; i++)
	//	printf("%d ", (*mats)[i]);
	
	double **concs = load_concs(num_nucs_p);

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

	//OpenCL initilization
	if( !init_ocl() )
		return -1;
 	
	unsigned long *vhash = NULL;
	posix_memalign((void**)&vhash, AOCL_ALIGN, in.nthreads*sizeof(unsigned long));
	double *energy = NULL;
	posix_memalign((void**)&energy, AOCL_ALIGN, n_iso_grid*sizeof(double));
	int *energy_xs = NULL;
	posix_memalign((void**)&energy_xs, AOCL_ALIGN, n_iso_grid*in.n_isotopes*sizeof(int));	
	int sample = 0;
	//printf("n_iso_grid is %ld\n", n_iso_grid);
	for(long il = 0; il < n_iso_grid; il++){
		energy[il] = energy_grid[il].energy;
		if( il == (long)((sample+1)*((double)n_iso_grid/(N+1))) )
		{
			firstN[sample].data = energy[il];
                        firstN[sample].index = il;
			//printf("firstN[sample].index is %ld.\n", firstN[sample].index);
                        sample++;
		}
		long j; 
		long index = il * in.n_isotopes;
		for(j = 0; j < in.n_isotopes; j++)
			energy_xs[index + j] = energy_grid[il].xs_ptrs[j];
	}
	
	cl_int status;
	cl_event write_events[5];
	cl_event kernel_event;
	cl_event finish_event;
	
	ocl_start_ts_sec = gettimesec();	

	status = clEnqueueWriteBuffer(queue, firstN_buf, CL_TRUE, 0, N*sizeof(CacheData), firstN, 0, NULL,&write_events[0]);
        checkError(status, "Failed to enqueue write buffer");	
/*	
	status = clEnqueueWriteBuffer(queue, num_nucs_buf, CL_TRUE, 0, 12*sizeof(int), num_nucs, 0, NULL, &write_events[0]);
	checkError(status, "Failed to enqueue write buffer.\n");
*/	
	status = clEnqueueWriteBuffer(queue, concs_buf, CL_TRUE, 0, num_dpoints*sizeof(double), *concs, 0, NULL, &write_events[1]);
        checkError(status, "Failed to enqueue write buffer.\n");

	status = clEnqueueWriteBuffer(queue, energy_grid_buf, CL_TRUE, 0, n_iso_grid*sizeof(double), energy, 0, NULL, &write_events[2]);
        checkError(status, "Failed to enqueue write buffer.\n");
	//printf("enqueued energy_grid_buf.\n");
	//printf("n_iso_grid = %d\n", n_iso_grid);
	//printf("n_isotopes = %d\n", in.n_isotopes);
	//for(i = 0; i < n_iso_grid; i++){

	//printf("before enqueue grid_xs.\n");	
	status = clEnqueueWriteBuffer(queue, energy_grid_xs_buf, CL_TRUE, 0, n_iso_grid*in.n_isotopes*sizeof(int), energy_xs, 0, NULL, NULL);
	//printf("i = %d\n", i);
	//}
	//printf("enqueued energy_grid_xs_buf.\n");
	//printf("before enqueue nuclide_grid_buf.\n");
	status = clEnqueueWriteBuffer(queue, nuclide_grids_buf, CL_TRUE, 0, n_iso_grid*sizeof(NuclideGridPoint), *nuclide_grids, 0, NULL, &write_events[3]);
        checkError(status, "Failed to enqueue write buffer.\n");
	//printf("before enqueue mats_buf.\n");
	status = clEnqueueWriteBuffer(queue, mats_buf, CL_TRUE, 0, num_dpoints*sizeof(int), *mats, 0, NULL, &write_events[4]);
        checkError(status, "Failed to enqueue write buffer.\n");

	//printf("after enqueue write buffer.\n");

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

	perthread_lookups = in.lookups/in.nthreads;
	status = clSetKernelArg(kernel, arg++, sizeof(unsigned int), &perthread_lookups);
	checkError(status, "Failed to set arg 1.1");

	status = clSetKernelArg(kernel, arg++, sizeof(int), &N);
        checkError(status, "Failed to set kernel arg 2: N");

        status = clSetKernelArg(kernel, arg++, sizeof(int), &STAGE);
        checkError(status, "Failed to set kernel arg 2: STAGE");

        status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &firstN_buf);
        checkError(status, "Failed to set kernel arg 2: firstN");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_int16), &num_nucs);
	checkError(status, "Failed to set kernel arg 2");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &concs_buf);
	checkError(status, "Failed to set arg 3");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &energy_grid_buf);
	checkError(status, "Faile to set arg 4");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &energy_grid_xs_buf);
	checkError(status, "Failed to set arg 4.1");
	
	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &nuclide_grids_buf);
	checkError(status, "Failed to set arg 5");

	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &mats_buf);
	checkError(status, "Failed to set arg 6");
//printf("before vhash buf arg setting.\n");
	status = clSetKernelArg(kernel, arg++, sizeof(cl_mem), &vhash_buf);
	checkError(status, "failed to ser arg 7");
	//flush queue	
	clFinish(queue);		
	
	ocl_kernel_ts_sec = gettimesec();
	snap_energy();

	//printf("after set kernel args\n");
	size_t global_work_size = in.nthreads; 
//	size_t local_work_size = in.lookups/in.nthreads;
//	size_t local_work_size = 1;

	status = clEnqueueNDRangeKernel(queue, kernel, 1, NULL,
       		 &global_work_size, NULL, 5, write_events, &kernel_event);
    	checkError(status, "Failed to launch kernel");
	//printf("enqueue kernel finish.\n");
	clFinish(queue);
	snap_energy();
	ocl_post_ts_sec = gettimesec();	
	//printf("kernel finished.\n");
	status = clEnqueueReadBuffer(queue, vhash_buf, CL_TRUE, 0, in.nthreads*sizeof(unsigned long), vhash, 1, &kernel_event, &finish_event);
			
	cl_uint t = 1;
        clWaitForEvents(t, &finish_event);
	ocl_end_ts_sec = gettimesec();	
	profiled_kernel_time_ns = getStartEndTime(kernel_event);

	clReleaseEvent(write_events[0]);
    	clReleaseEvent(write_events[1]);
 	clReleaseEvent(write_events[2]);
    	clReleaseEvent(write_events[3]);
	clReleaseEvent(write_events[4]);
		
	cleanup();
	end_ts_sec = gettimesec();

	#ifdef VERIFICATION 
	unsigned long * vhash_v = (unsigned long*)malloc(sizeof(unsigned long)*in.nthreads);
	for(i = 0; i < in.nthreads; i++)
	{
		ulong seed = (i+1)*19+17;
		double p_energy;
        	int mat;
       		double macro_xs_vector[5] = {0};
		unsigned int hash = 5381;       
        	vhash_v[i] = 0;

		for(int n = 0; n < perthread_lookups; n++)
 		{
			p_energy = rn(&seed);
        		mat = pick_mat(&seed);
        		double xs_vector[5];	
        		int p_nuc;	
			long idx = 0;   
        		double conc;
			idx = grid_search( in.n_isotopes * in.n_gridpoints, p_energy,
                        	   	   energy_grid);
		//	printf("idx_v is %ld\n", idx);

			for( int j = 0; j < num_nucs_p[mat]; j++ )
	       		{       
        	        	p_nuc = mats[mat][j];
                		conc = concs[mat][j];
        //			printf("Main: p_nuc is %d, conc is %f\n", p_nuc, conc );
	        		calculate_micro_xs( p_energy, p_nuc, in.n_isotopes,
                                	    	    in.n_gridpoints, energy_grid,
                                    		    nuclide_grids, idx, xs_vector );
                		for( int k = 0; k < 5; k++ )
                        		macro_xs_vector[k] += xs_vector[k] * conc;
        		}
		}
		//	printf("p_energy_v is %f, mat_v is %d\n", p_energy, mat);		
		hash = ((hash << 5) + hash) + (int)p_energy;
                hash = ((hash << 5) + hash) + (int)mat;
                for(int k = 0; k < 5; k++)
                       	hash = ((hash << 5) + hash) + macro_xs_vector[k];
                vhash_v[i] = hash % 10000;	
	}
	
//	for(i = 0; i < in.nthreads; i++){
//		printf("i: %d, vhash is %ld, vhash_v is %ld\n", i, vhash[i], vhash_v[i]);
//	}
	bool pass = true;
	for(i = 0; (i < in.nthreads) && pass; i++){
		if(vhash[i] != vhash_v[i]){
			printf("i: %d, vhash is %ld, vhash_v is %ld\n", i, vhash[i], vhash_v[i]);
			pass = false;
			printf("Veification FAIL.\n");
			break;
		}
	}
	if(pass == true)	
		printf("Verification PASS.\n");
	#endif 
//	memcpy(xs, macro_xs_vector, 5*sizeof(double));

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
	
//	print_results( in, mype, omp_end-omp_start, nprocs, vhash );
/*
	#ifdef BENCHMARK
	}
	#endif
*/

	//result print
	if (mype == 0) {
		printf("\n" );
                printf("Simulation complete.\n" );
		
		border_print();
                center_print("RESULTS", 79);
                border_print();	
		
		printf("START_TS_SEC=%f\n", start_ts_sec);
                printf("OCL_START_TS_SEC=%f\n", ocl_start_ts_sec);
		printf("OCL_KERNEL_TS_SEC=%f\n", ocl_kernel_ts_sec);
                printf("OCL_POST_TS_SEC=%f\n", ocl_post_ts_sec);
                printf("OCL_END_TS_SEC=%f\n", ocl_end_ts_sec);
                printf("END_TS_SEC=%f\n", end_ts_sec);
		printf("PROFILED_KERNEL_TIME_NSEC=%lu\n", profiled_kernel_time_ns);

               // printf("ARRAYINIT_SEC=%lf\n", nthrun, arrayinit_end_ts_sec - arrayinit_start_ts_sec);
                printf("OCL_BEFORE_SEC=%lf\n", ocl_kernel_ts_sec - ocl_start_ts_sec);
                printf("OCL_KERNEL_SEC=%lf\n", profiled_kernel_time_ns*1e-9);
                printf("OCL_AFTER_SEC=%lf\n",  ocl_end_ts_sec - ocl_kernel_ts_sec - profiled_kernel_time_ns*1e-9);
                printf("TOTAL_OCL_SEC=%lf\n", ocl_end_ts_sec - ocl_start_ts_sec);
                printf("TOTAL_RUNNING_SEC=%lf\n", end_ts_sec-start_ts_sec);
		printf("BANDWIDTH_LOOKUPS/Sec=%lf\n", bw());
                report_energy();
        }

	return 0;
}

static bool init_ocl()
{
	cl_int status;
        if(!setCwdToExeDir())
		return false;

#if defined(USE_INTELOPENCL)
        platform = findPlatform("Intel");
#else
        platform = findPlatform("Altera");
#endif	
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
/*
	cl_ulong buffer;
    	clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(buffer), &buffer, NULL);
    	printf("CL_DEVICE_MAX_WORK_GROUP_SIZE = %llu \n", (unsigned long long)buffer);
*/	
	queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &status);
	checkError(status, "Failed to create command queue.\n");
	
#if defined(USE_INTELOPENCL)
	// build ocl kernel directly
	char *kernel_text;
	size_t text_size;
	FILE *fp;
	char fn[1024];
	snprintf(fn, 1024, "optbsearch/optbCalculateXS.cl");

	fp = fopen(fn, "r");
	if (!fp) {
		printf("failed to open: %s\n", fn);
		exit(1);
	}

	fseek(fp, 0, SEEK_END);
        text_size = ftell(fp);
        fseek(fp, 0, SEEK_SET);
        kernel_text = (char *)malloc(text_size+1);
	if (!kernel_text) {
		printf("%s(%d)\n", __FILE__, __LINE__);
		exit(1);
	}
	if (fread(kernel_text, text_size, 1, fp) == 0) {
		printf("%s(%d)\n", __FILE__, __LINE__);
		exit(1);
	}
	fclose(fp);
	kernel_text[text_size] = 0;
	//printf("%s\n", kernel_text);

	cl_program program = clCreateProgramWithSource(context, 1, (const char**)&kernel_text, 0, &status);
	status = clBuildProgram(program, 1, &device, "", 0, 0);
	if (status != CL_SUCCESS) {
		size_t sz = 0, retsz;
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &sz);

		char *log = (char*)malloc(sz+1);
		if (log) {
			clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sz, log, &retsz);
			log[retsz] = 0;
			printf(log);
		}

		exit(1);
	}

#else
	std::string binary_file = getBoardBinaryFile("optbCalculateXS", device);
	printf("Using AOCX: %s\n", binary_file.c_str());
	program = createProgramFromBinary(context, binary_file.c_str(), &device, 1);

  	status = clBuildProgram(program, 0, NULL, "", NULL, NULL);
	checkError(status, "Failed to build program.\n");
#endif
	
	const char *kernel_name = "calculate_macro_xs";
	kernel = clCreateKernel(program, kernel_name, &status);
	checkError(status, "Failed to create kernel.\n");

	rapl_reader_init();
	rapl_reader_snap();

	//Give buffer 
	vhash_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, in.nthreads*sizeof(unsigned long), NULL, &status);
	checkError(status, "Failed to create output buffer.\n");
/*
       	num_nucs_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, 12*sizeof(int), NULL, &status);
	checkError(status, "Failed to create num_nucs input buffer.\n");
*/
	concs_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, num_dpoints*sizeof(double), NULL, &status);
	checkError(status, "Failed to create concs input buffer.\n");
	
	energy_grid_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, n_iso_grid*sizeof(double), NULL, &status);
	checkError(status, "Failed to create energy_grid input buffer.\n");
	
	energy_grid_xs_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, n_iso_grid*in.n_isotopes*sizeof(int), NULL, &status);
	checkError(status, "Failed to create input energy_grid_xs buffer. \n");

	nuclide_grids_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, n_iso_grid*sizeof(NuclideGridPoint), NULL, &status);
	checkError(status, "Failed to create nuclide_grids input buffer.\n");

	mats_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, num_dpoints*sizeof(int), NULL, &status);
	checkError(status, "Failed to create mats input buffer.\n");

	firstN_buf = clCreateBuffer(context, CL_MEM_READ_ONLY, N*sizeof(CacheData), NULL, &status);
        checkError(status, "Failed to create buffer: firstN_buf");
	
	return true;
}

static int get_num_datapoints(cl_int16 num_nucs)
{
//	int num_datapoints = 0;
//	for(int i = 0; i < 12; i++)
	int num_datapoints = num_nucs.s0+num_nucs.s1+num_nucs.s2+num_nucs.s3+num_nucs.s4+num_nucs.s5
			+num_nucs.s6+num_nucs.s7+num_nucs.s8+num_nucs.s9+num_nucs.sA+num_nucs.sB;

	return num_datapoints;
}

void cleanup()
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

long grid_search( long n, double quarry, GridPoint * A)
{
        long lowerLimit = 0;
        long upperLimit = n-1;
        long examinationPoint;
        long length = upperLimit - lowerLimit;

        while( length > 1 )
        {
                examinationPoint = lowerLimit + ( length / 2 );

                if( A[examinationPoint].energy > quarry )
                        upperLimit = examinationPoint;
                else
                        lowerLimit = examinationPoint;

                length = upperLimit - lowerLimit;
        }

        return lowerLimit;
}

// Calculates the microscopic cross section for a given nuclide & energy
void calculate_micro_xs(double p_energy, int nuc, long n_isotopes,
                        long n_gridpoints,
                        GridPoint * energy_grid,
                        NuclideGridPoint ** nuclide_grids,
                        int idx, double * xs_vector )
{	
	// Variables
	double f;
	NuclideGridPoint * low, * high;

	// pull ptr from energy grid and check to ensure that
	// we're not reading off the end of the nuclide's grid
	if( energy_grid[idx].xs_ptrs[nuc] == n_gridpoints - 1 )
		low = &nuclide_grids[nuc][energy_grid[idx].xs_ptrs[nuc] - 1];
	else
		low = &nuclide_grids[nuc][energy_grid[idx].xs_ptrs[nuc]];
	
	high = low + 1;
	
	// calculate the re-useable interpolation factor
	f = (high->energy - p_energy) / (high->energy - low->energy);

	// Total XS
	xs_vector[0] = high->total_xs - f * (high->total_xs - low->total_xs);
	
	// Elastic XS
	xs_vector[1] = high->elastic_xs - f * (high->elastic_xs - low->elastic_xs);
	
	// Absorbtion XS
	xs_vector[2] = high->absorbtion_xs - f * (high->absorbtion_xs - low->absorbtion_xs);
	
	// Fission XS
	xs_vector[3] = high->fission_xs - f * (high->fission_xs - low->fission_xs);
	
	// Nu Fission XS
	xs_vector[4] = high->nu_fission_xs - f * (high->nu_fission_xs - low->nu_fission_xs);	
}
/*	
double gettimesec(void)
{
	struct timeval tv;
        gettimeofday(&tv, 0);
        return (double)tv.tv_sec + (double)tv.tv_usec/1000.0/1000.0;
}
*/
double bw(void) 
{
	return 1e9*(double)in.lookups / (double)profiled_kernel_time_ns;
}

void snap_energy(void) {
	rapl_reader_snap();
}

void report_energy(void) {
        double delta;
        double pkg[8],pp0[8],pp1[8],mem[8];

        rapl_reader_get_energy(&delta, pkg, pp0, pp1, mem);
        printf("RAPL_CPU0_ENERGY_J=%lf # NOTE: RAPL conditionally works\n",
                       pkg[0] + mem[0]);
        printf("RAPL_CPU1_ENERGY_J=%lf\n",
                       pkg[1] + mem[1]);
        printf("RAPL_CPU_ENERGY_DELTA_SEC=%lf\n", delta);
}
