# Makefile for use with the USI ICS HPC cluster
#
# Tested as of 2016.11.30
#
# Required modules:
#   petsc

ALL: main

# We use the bash shell to allow for our checks below
SHELL=/bin/bash

OBJ = main.o system.o dump.o

main: $(OBJ) 
	${CLINKER} -o $@ $^ ${PETSC_LIB}

clean ::
	rm -f main $(OBJ) out.m output.*

# Some tests, to be run in an interactive session (see salloc)
test : run_test1 run_test2

# Some PETSc command line options
TEST_OPTIONS=-ts_monitor -snes_monitor -ksp_monitor -ts_view

run_test1: ${EXNAME}
	@rm -f test1.tmp
	@echo "Running Test 1"
	@mpiexec -n 1 ./main ${TEST_OPTIONS}  \
    2>&1 > test1.tmp
	@diff test1.tmp testref/test1.ref && \
    echo "Success" || \
    echo "Check: output does not match reference (see diff above)"
	@rm -f test1.tmp

run_test2: ${EXNAME}
	@rm -f test2.tmp
	@echo "Running Test 2"
	@mpiexec -n 2 ./main ${TEST_OPTIONS}  \
    2>&1 > test2.tmp
	@diff  test2.tmp testref/test2.ref && \
    echo "Success" || \
    echo "Check: output does not match reference (see diff above)"
	@rm -f test2.tmp

.PHONY: ALL allclean test run_test1 run_test2

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
