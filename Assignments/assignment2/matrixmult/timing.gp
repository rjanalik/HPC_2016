set title "NxN matrix-matrix-multiplication on Quad-Core AMD Opteron(tm) Processor 2344@1.7GH"
set xlabel "Matrix size (N)"
set ylabel "Performance (MFlop/s)"
set grid

set terminal postscript color "Helvetica" 14
set output "timing.ps"

# set terminal png color "Helvetica" 14
# set output "timing.png"

# plot "timing.data" using 2:4 title "square_dgemm" with linespoints


# For performance comparisons

plot "timing_basic_dgemm.data"   using 1:2 title "Naive dgemm" with linespoints, \
     "timing_blocked_dgemm.data" using 1:2 title "Student dgemm" with linespoints, \
     "timing_blas_dgemm.data"   using 1:2 title "GOTO BLAS dgemm" with linespoints
