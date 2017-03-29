The codes that used for FCCM paper is in XSBench/singlegroup_src/ folder
1. Make files
1) Makefile //compile code to runnable for XSBench (no optimization)
2) Makefile.optb //compile code to runnable for optimized binary search and vector operation case (there is compile option VECTOR_ADD, you can set/unset it to have different binaries)
3) Makefile.intel //silimar to Makefile, but will generate runnable to run OpenCL code on CPU 
4) Makefile.intel.optb //similar to Makefile.optb, but will generate runnable to run OpenCL code on CPU 

2. AOCX binaries
There’re pre-compiled docx binaries with different optimizations in /var/tmp/XSBench_aocx/ folder on duteros.

3. XSBench
XSBench runnable has a series of input arguments that you can specify:
-t <threads>     Number of OpenCL work items to run
-s <size>        Size of H-M Benchmark to run (small, large, XL, XXL)
-g <gridpoints>  Number of gridpoints per nuclide (overrides -s defaults)
-l <lookups>     Number of Cross-section (XS) lookups
Default is equivalent to: -s large -l 15000000

4. Running the code
When you want to run XSBench, you 
1) compile the code using make files;
2) copy the corresponding aocx file into current location;
3) excute the runnable file. 
For example, if you want to run the base case for XSBench, do the following steps:
$ make
$ cp /var/tmp/XSBench_aocx/base/CalculateXS.aocx .
$./XSBench -t 15000000

5. runxs.sh
This script is used to automate the running process after you compiled the code with make files.
It also captures the power data of FPGA when running OpenCL kernel. 
