#include "petsc.h"

int main(int argc, char ** argv){

  PetscErrorCode ierr;
  PetscMPIInt rank,size;

  ierr = PetscInitialize(&argc,&argv,NULL,NULL);CHKERRQ(ierr);
  
  MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Hello, World!\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC_COMM_WORLD has %d ranks\n",size);CHKERRQ(ierr);

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD," Hello from rank %d\n",rank);CHKERRQ(ierr);
  ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);

  PetscFinalize();

  return 0;
}
