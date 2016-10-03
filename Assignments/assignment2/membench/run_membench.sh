#!/bin/bash -l

#SBATCH --job-name=membench
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --output=membench-%j.out
#SBATCH --error=membench-%j.err

# load modules
module load gcc/6.1.0 gnuplot

./membench | sed -e '/:/	s//: /g' -e '/  */	s//	/g' | cut -f2,4,6 > generic.xxx && sed -e '/sarlacc/ s//generic/' gnuplot.template > generic.gp && gnuplot generic.gp
