SHELL=/bin/bash

MPI_CC=mpiicc

# The array size.
N= 10000 #need revise to meet the requirement

# The length where to use the openmpi
BEGIN_LENGTH_OMP= 1000 #need revise to meet the requirement
BEGIN_LENGTH_MPI= 0 #need revise to meet the requirement
# How many iterations to run
ITERS=20 #need revise to meet the requirement

# Boundary values.
TOP_BOUNDARY_VALUE=5.8
BOTTOM_BOUNDARY_VALUE=10.2
LEFT_BOUNDARY_VALUE=4.3
RIGHT_BOUNDARY_VALUE=9.2

# Number of threads.
NUM_THREADS_TAUB=12 #need revise to meet the requirement

# Compiler optimization level.
OPT_LEVEL=-O0 #-qopt-reportP#need revise to meet the requirement



##########################################
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
# DO NOT MODIFY ANYTHING BELOW THIS LINE #
##########################################

# Used for checking the results.
ERROR_THRESHOLD=1e-3

# The papi library location.

# Common program arguments.
COMMON_PROG_ARGS=-DN=$(N) \
			 		  		 -DITERS=$(ITERS) \
			        	 -DTOP_BOUNDARY_VALUE=$(TOP_BOUNDARY_VALUE) \
			        	 -DBOTTOM_BOUNDARY_VALUE=$(BOTTOM_BOUNDARY_VALUE) \
			        	 -DLEFT_BOUNDARY_VALUE=$(LEFT_BOUNDARY_VALUE) \
			        	 -DRIGHT_BOUNDARY_VALUE=$(RIGHT_BOUNDARY_VALUE) \
			        	 -DERROR_THRESHOLD=$(ERROR_THRESHOLD) \
								 -std=c99 \
								 -qopenmp#need revise to meet the requirement


# Program arguments for Taub.
TAUB_PROG_ARGS=$(COMMON_PROG_ARGS) \
  			     	 -DNUM_THREADS=$(NUM_THREADS_TAUB) \
  			     	 -DLENGTH_OMP=$(BEGIN_LENGTH_OMP) \
				 -DLENGTH_MPI=$(BEGIN_LENGTH_MPI) \

# Compilation command for Taub, no PAPI.
TAUB_MPICC=$(MPI_CC) $(TAUB_PROG_ARGS)

all: taub_hybrid

taub_hybrid: *.c
	$(TAUB_MPICC) *.c -o taub_hybrid

clean:
	rm -f taub*
