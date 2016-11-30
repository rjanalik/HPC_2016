EXNAME=runMe

all : ${EXNAME}

${EXNAME}: code.o  chkopts
	-${CLINKER} -o ${EXNAME} code.o  ${PETSC_LIB}

clean ::
	rm -f test *.o

test :
	mpirun -n 1 ./runMe
	mpirun -n 2 ./runMe

test_batch :
	sbatch prun1 runMe
	sbatch prun8 runMe

.PHONY : clean runMe test test_batch

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
