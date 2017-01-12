#include "../Types.h"
/*
typedef struct __attribute__((packed)) __attribute__((aligned(64))){
	double energy;
	double total_xs;
	double elastic_xs;
	double absorbtion_xs;
	double fission_xs;
	double nu_fission_xs;
} NuclideGridPoint;
*/
// Calculates the microscopic cross section for a given nuclide & energy
void calculate_micro_xs(   double p_energy, int nuc, long n_isotopes,
                           long n_gridpoints, 
			   __global int * restrict energy_grid_xs,
                           __global NuclideGridPoint * restrict nuclide_grids,
                           int idx, 
                           double * restrict xs_vector ){
	
	// Variables
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable	
	double f;
	NuclideGridPoint low, high;
	// pull ptr from energy grid and check to ensure that
	// we're not reading off the end of the nuclide's grid
        long index = nuc * n_gridpoints;
       	long index_xs = idx * n_isotopes + nuc;

	if(energy_grid_xs[index_xs] == n_gridpoints - 1 ){
		low = nuclide_grids[index + energy_grid_xs[index_xs] - 1];
	        high = nuclide_grids[index + energy_grid_xs[index_xs]];
	}else{
		low = nuclide_grids[index + energy_grid_xs[index_xs]];
                high = nuclide_grids[index + energy_grid_xs[index_xs] + 1];	
	}
	// calculate the re-useable interpolation factor
	f = (high.energy - p_energy) / (high.energy - low.energy);

	// Total XS
	xs_vector[0] = mad( -f, (high.total_xs - low.total_xs), high.total_xs );
	
	// Elastic XS
	xs_vector[1] = mad( -f, (high.elastic_xs - low.elastic_xs), high.elastic_xs );
	
	// Absorbtion XS
	xs_vector[2] = mad( -f, (high.absorbtion_xs - low.absorbtion_xs), high.absorbtion_xs );
	
	// Fission XS
	xs_vector[3] = mad( -f, (high.fission_xs - low.fission_xs), high.fission_xs );
	
	// Nu Fission XS
	xs_vector[4] = mad( -f, (high.nu_fission_xs - low.nu_fission_xs), high.nu_fission_xs );
	
}

// (fixed) binary search for energy on unionized energy grid
// returns lower index
long grid_search( long total, double quarry, __global double *restrict A, __constant CacheData *restrict firstNdata, const int N, const int STAGE)
{
	int quick_low = 0;
        int quick_high = N - 1;
        int quick_length = quick_high - quick_low;
        int quick_mid;
        int i = 0;

	long lowerLimit = 0;
	long upperLimit = total-1;
	long examinationPoint;
	long length;
	
	while( i < STAGE )
        {       
                quick_mid = quick_low + (quick_length>>1);
                CacheData exmain_data = firstNdata[quick_mid]; 
                if(exmain_data.data > quarry)
                {       
                        upperLimit = exmain_data.index;
                        quick_high = quick_mid;
                }
                else if(exmain_data.data == quarry)
                        return exmain_data.index;
                else
                {       
                        lowerLimit = exmain_data.index;
                        quick_low = quick_mid; 
                }
                quick_length = quick_high - quick_low;
                i++;
        }

        length = upperLimit - lowerLimit; 

	while( length > 1 )
	{
		examinationPoint = lowerLimit + (length >> 1);	

		if( A[examinationPoint] > quarry )
			upperLimit = examinationPoint;
		else
			lowerLimit = examinationPoint;
		
		length = upperLimit - lowerLimit;
	}
	
	return lowerLimit;
}

double rn(unsigned long * seed)
{
	double ret;
        unsigned long n1;
        unsigned long a = 16807;
        unsigned long m = 2147483647;
        n1 = ( a * (*seed) ) % m;
        *seed = n1;
        ret = (double) n1 / m;
        return ret;
}

// picks a material based on a probabilistic distribution
int pick_mat( unsigned long * seed )
{
	// I have a nice spreadsheet supporting these numbers. They are
	// the fractions (by volume) of material in the core. Not a 
	// *perfect* approximation of where XS lookups are going to occur,
	// but this will do a good job of biasing the system nonetheless.

	// Also could be argued that doing fractions by weight would be 
	// a better approximation, but volume does a good enough job for now.

	double dist[12];
	dist[0]  = 0.140;	// fuel
	dist[1]  = 0.052;	// cladding
	dist[2]  = 0.275;	// cold, borated water
	dist[3]  = 0.134;	// hot, borated water
	dist[4]  = 0.154;	// RPV
	dist[5]  = 0.064;	// Lower, radial reflector
	dist[6]  = 0.066;	// Upper reflector / top plate
	dist[7]  = 0.055;	// bottom plate
	dist[8]  = 0.008;	// bottom nozzle
	dist[9]  = 0.015;	// top nozzle
	dist[10] = 0.025;	// top of fuel assemblies
	dist[11] = 0.013;	// bottom of fuel assemblies
	
	double roll = rn(seed);
	// makes a pick based on the distro
	for( int i = 0; i < 12; i++ )
	{
		double running = 0;
		for( int j = i; j > 0; j-- )
			running += dist[j];
		if( roll < running )
			return i;
	}

	return 0;
}

// Calculates macroscopic cross section based on a given material & energy 
__kernel void calculate_macro_xs(
			 const uint verification,  
		         const long n_isotopes, 
		 	 const long n_gridpoints,
			 const uint per_thread_lookups,
                         const int N_cache_datapoints,
                         const int search_stages,
			 __constant CacheData *restrict first_N_datapoints,
			 __global int * restrict num_nucs,
                         __global double * restrict concs,
                         __global double * restrict energy_grid,
                         __global int * restrict energy_grid_xs,
			 __global NuclideGridPoint * restrict nuclide_grids,
                         __global int * restrict mats,
                         __global ulong * restrict vhash )
{
	int thread = get_global_id(0);
	ulong seed = (thread+1)*19+17;

	#pragma OPENCL EXTENSION cl_khr_fp64 : enable			
	double p_energy;
	int mat;
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
	double macro_xs_vector[5] = {0};
	
	for(int n = 0; n < per_thread_lookups; n++)
	{	
		p_energy = rn(&seed);
		mat = pick_mat(&seed);
	
		#pragma OPENCL EXTENSION cl_khr_fp64 : enable	
		double xs_vector[5];
 
		int p_nuc; // the nuclide we are looking up
		long idx = 0;	
		#pragma OPENCL EXTENSION cl_khr_fp64 : enable
		double conc; // the concentration of the nuclide in the material
	
		// binary search for energy on unionized energy grid (UEG)
		idx = grid_search( n_isotopes * n_gridpoints, p_energy,
	        	           energy_grid, first_N_datapoints, 
				   N_cache_datapoints, search_stages);		

        	//calculate the index for mats and concs from 2d array
		int index = 0;
        	for(int i = 0; i < mat; i++ )
        		index += num_nucs[i];
		// Once we find the pointer array on the UEG, we can pull the data
		// from the respective nuclide grids, as well as the nuclide
		// concentration data for the material
		// Each nuclide from the material needs to have its micro-XS array
		// looked up & interpolatied (via calculate_micro_xs). Then, the
		// micro XS is multiplied by the concentration of that nuclide
		// in the material, and added to the total macro XS array.	
		for( int j = 0; j < num_nucs[mat]; j++ )
		{ 
			int index_j = index + j;
			p_nuc = mats[index_j];
                	conc = concs[index_j];
			calculate_micro_xs( p_energy, p_nuc, n_isotopes,
		        	            n_gridpoints, energy_grid_xs,
		                	    nuclide_grids, idx, xs_vector );	
			#pragma unroll 
			for( int k = 0; k < 5; k++ )
 				macro_xs_vector[k] = mad( xs_vector[k], conc, macro_xs_vector[k] );
		}
	}

	if(verification == 1){
	//	printf("idx is %ld, p_energy is %f, mat is %d\n", idx, p_energy, mat);
		unsigned int hash = 5381;	
		vhash[thread] = 0;
		hash = ((hash << 5) + hash) + (int)p_energy;
		hash = ((hash << 5) + hash) + (int)mat;
		for(int k = 0; k < 5; k++)
			hash = ((hash << 5) + hash) + macro_xs_vector[k];
		vhash[thread] = hash % 10000;
	}	
}

