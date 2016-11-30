all : ex23 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

ex23: ex23.o  chkopts
	-${CLINKER} -o ex23 ex23.o  ${PETSC_KSP_LIB}
	${RM} ex23.o
