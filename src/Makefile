#===============================================================================
# User Options
#===============================================================================

COMPILER    = gnu
OPTIMIZE    = yes
DEBUG       = no
PROFILE     = no
MPI         = no
PAPI        = no
VEC_INFO    = no
VERIFY      = no
BENCHMARK   = no
BINARY_DUMP = no
BINARY_READ = no

#===============================================================================
# Program name & source code list
#===============================================================================

program = XSBench

ifeq ($(wildcard $(ALTERAOCLSDKROOT)),)
$(error Set ALTERAOCLSDKROOT to the root directory of the Altera SDK for OpenCL software installation)
endif
ifeq ($(wildcard $(ALTERAOCLSDKROOT)/host/include/CL/opencl.h),)
$(error Set ALTERAOCLSDKROOT to the root directory of the Altera SDK for OpenCL software installation.)
endif

COMMON_FILES = ../common/src/AOCL_Utils.cpp
AOCL_COMPILE_CONFIG=$(shell aocl compile-config) -I../common/inc
AOCL_LINK_CONFIG=$(shell aocl link-config) -lnalla_pcie_mmd

source = \
Main.cpp \
io.c \
GridInit.c \
XSutils.c \
Materials.c

obj = Main.o io.o GridInit.o XSutils.o Materials.o 
#obj = $(source:.c=.o .cpp=.o)
#===============================================================================
# Sets Flags
#===============================================================================

# Standard Flags
CFLAGS := -std=gnu99

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  CFLAGS += -fopenmp
endif

# Intel Compiler
ifeq ($(COMPILER),intel)
  CC = icc
  CFLAGS += -openmp 
endif

# BG/Q gcc Cross-Compiler
ifeq ($(MACHINE),bluegene)
  CC = mpicc
endif

# Debug Flags
ifeq ($(DEBUG),yes)
  CFLAGS += -g
  LDFLAGS  += -g
endif

# Profiling Flags
ifeq ($(PROFILE),yes)
  CFLAGS += -pg
  LDFLAGS  += -pg
endif

# Optimization Flags
ifeq ($(OPTIMIZE),yes)
  CFLAGS += -O3
endif

# Compiler Vectorization (needs -O3 flag) information
ifeq ($(VEC_INFO),yes)
  CFLAGS += -ftree-vectorizer-verbose=6
endif

# PAPI source (you may need to provide -I and -L pointing
# to PAPI depending on your installation
ifeq ($(PAPI),yes)
  source += papi.c
  CFLAGS += -DPAPI
  #CFLAGS += -I/soft/apps/packages/papi/papi-5.1.1/include
  #LDFLAGS += -L/soft/apps/packages/papi/papi-5.1.1/lib -lpapi
  LDFLAGS += -lpapi
endif

# MPI
ifeq ($(MPI),yes)
  CC = mpicc
  CFLAGS += -DMPI
endif

# Verification of results mode
ifeq ($(VERIFY),yes)
  CFLAGS += -DVERIFICATION
endif

# Adds outer 'benchmarking' loop to do multiple trials for
# 1 < threads <= max_threads
ifeq ($(BENCHMARK),yes)
  CFLAGS += -DBENCHMARK
endif

# Binary dump for file I/O based initialization
ifeq ($(BINARY_DUMP),yes)
  CFLAGS += -DBINARY_DUMP
endif

# Binary read for file I/O based initialization
ifeq ($(BINARY_READ),yes)
  CFLAGS += -DBINARY_READ
endif


#===============================================================================
# Targets to Build
#===============================================================================

#$(program): $(obj) Types.h XSbench_header.h
#	g++ -fopenmp $(obj) $(COMMON_FILES) -o $@ $(LDFLAGS) $(AOCL_COMPILE_CONFIG) $(AOCL_LINK_CONFIG)

#%.o: %.c %.cpp
#	$(CC) $(CFLAGS) -c $< -o $@ 

all:
	g++ -fopenmp $(source) $(COMMON_FILES) -o $(program) $(LDFLAGS) $(AOCL_COMPILE_CONFIG) $(AOCL_LINK_CONFIG)
clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) papi.c XSbench_header.h

run:
	./$(program)

bgqrun:
	qsub -t 10 -n 1 -O test XSBench
