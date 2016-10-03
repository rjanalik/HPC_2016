#!/bin/sh
./basic_dgemm && ./blas_dgemm && ./blocked_dgemm && gnuplot timing.gp

