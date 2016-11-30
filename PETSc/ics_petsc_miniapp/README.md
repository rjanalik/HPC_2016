This is material adapted from the CSCS Summer School 2016, largely by Ben Cumming
and Gilles Fourestey.
https://github.com/eth-cscs/SummerSchool2016

This is a version of the miniapp which uses PETSc. It is different in spirit from
the other versions of the miniapp, in that it is not intended to be an 
essentially-identical piece of code with to-be-ported kernels, but is rather
an example of writing an equivalent code using a higher-level library.

For more information, see the notes in the .c files here

## Quick start on USI ICS Cluster
Load the PETSc module and the python module

    module load PETSc
    module load python

Navigate to this directory 


Build the executable
    make

Test on the login node
    make test

Test in an interactive session

    salloc
    module load petsc
    make test

You should see
    Running Test 1
    Success
    Running Test 2
    Success

Run your own experiments in the interactive session

    mpiexec -n 4 ./main -nx 99 -ny 88 -ts_monitor -snes_monitor -ksp_monitor 
    mpiexec -n 1 ./main -nx 16 -ny 16 -ts_monitor -snes_monitor -ksp_monitor -assemble 1 -pc_type gamg -dump 1


   To view the `.bov` file that is generated (only for single-processor runs with the `-dump` option), we borrow the procedure from the CSCS miniapp

    module load python
    python plotting.py

  You may get some harmless warnings about fonts.

  To view the generated `output.png`, copy it to your local machine and open it
