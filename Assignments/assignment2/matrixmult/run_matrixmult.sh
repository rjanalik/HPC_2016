#!/bin/bash -l

#SBATCH --job-name=matrixmult
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --output=matrixmult-%j.out
#SBATCH --error=matrixmult-%j.err

# load modules
module load gcc/6.1.0 openblas/0.2.18_gcc-6.1 gnuplot

./basic_dgemm && ./blas_dgemm && ./blocked_dgemm && gnuplot timing.gp
