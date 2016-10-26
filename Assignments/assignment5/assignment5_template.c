/*******************************************************************************
*   High Performance Computing, Fall 2016
*   Assignment 6
*   Institute of Computational Science, Universita della Svizzera italiana
*   Olaf Schenk
********************************************************************************
*/
/*
   Example.
   ======================
 
   The program computes the solution to the system of linear
   equations with a square matrix A and multiple
   right-hand sides B, where A is the coefficient matrix:
 
   Exercise 1.
   ===========

     Use the following 5x5 matrix A 

           6.80 -0.11  0.66  1.97  1.23
          -0.11  3.30  0.36 -0.44  1.08
      A = -0.66  0.36  5.70  0.27  1.04
           1.97 -0.44  0.27  7.17  0.14
           1.23  1.08  1.04  0.14  6.87

     and B is the right-hand side matrix which consists of three
     right right-hand sides:
 
          4.02  -1.56   9.81 
          6.19   4.00  -4.09 
     B = -8.22  -8.67  -4.57 
         -7.57   1.75  -8.61 
         -3.03   2.86   8.99 
 
     to compute a solution of A*X=B. 

      (1.a) Use the Gaussian Elimination method (LAPACK Routine dgesv,
           column-major, high-level) to compute a solution of A*X=B 

      (1.b) Print the 2-norm of the residual for the right-hand side at 
            the end of the solution process. Use the BLAS routines
            DGEMV and res = sdot (n, x, incx, y, incy)

      Technical Details:
      ====================
 
      Description of DGESV (column-major)
      The LAPACK Routine dgesv solves for X the system of linear 
      equations A*X = B,  where A is an n-by-n matrix, the columns 
      of matrix B are individual right-hand sides, and the columns 
      of X are the corresponding solutions.

      The LU decomposition with partial pivoting and row interchanges
      is used to factor A as A = P*L*U, where P is a permutation matrix, 
      L is unit lower triangular, and U is upper triangular. The factored 
      form of A is then used to solve the system of equations A*X = B.

      Solution:

           6.80, -0.11,  0.66,  1.97,  1.23,
          -0.11,  3.30,  0.36, -0.44,  1.08,
      X = -0.66,  0.36,  5.70,  0.27,  1.04,
           1.97, -0.44,  0.27,  7.17,  0.14,
           1.23,  1.08,  1.04,  0.14,  6.87

      Details of LU factorization:
  
           6.80  -0.11  -0.66   1.97   1.23
          -0.02   3.30   0.35  -0.41   1.10
    L/U =  0.10   0.11   5.72   0.12   0.80
           0.29  -0.12   0.09   6.54  -0.15
           0.18   0.33   0.18  -0.02   6.13

      Pivot indices P = [ 1 2 3 4 5]
 
*/

 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

void * _smalloc(size_t  count,  size_t  size)
{
        void * new = malloc (count * size);
        if (new == NULL)
        {
                printf("ERROR during memory allocation!\n");
                exit (7);
        }
        return  new;
}

#define  mem_free(var)                  do { free(var); var = NULL; } while (0)
#define  mem_alloc(var, cnt, typ) var = (typ *) _smalloc  ((cnt), sizeof(typ))


/* DGESV prototype */
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                   double* b, int* ldb, int* info );

/* DGEEV prototype */
extern void dgeev_(char *jobv11l, char *jobvr, int* n, double* a, int* lda, 
                   double* wr, double* wi, double* vl, int* ldvl, double* vr, 
                   int* ldvr, double* work, int* lwork, int* info);

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, double* a, int lda );
extern void print_int_vector( char* desc, int n, int* a );
extern void print_double_vector( char* desc, int n, double* a );

/* Parameters */
#define N 5
#define R 3

/* Main program */
int main() 
{
	/* Locals */
	int   _n[11] = {5, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
        int _rhs[11] = {3,   1,    1,    1,    1,    1,    1,    1,    1,    1,    1}; 
        int nrhs;
        int n;                    //      n: size of A; 
        int lda, ldb;             //    ldx: leading dimensions of a and b
        int info;                 //   info: output parameter for SCALAPACK functions
        int ldvr, lwork;          //  

# ifdef _OPENMP
        double wt1,wt2;
# endif

	/* Local arrays */
        int    *ipiv 	= NULL;   //   ipiv: permutation applied on the columns of a
	double *a 	= NULL;   //      a: matrix
	double *a_copy	= NULL;   // a_copy: copy of a
	double *b 	= NULL;   //      b: Right-Hand-Sides
	double *b_copy 	= NULL;   // b_copy: copy of b
	double *res 	= NULL;   //    res: residual
	double *x 	= NULL;   //      x: solution of the system
	double *x1 	= NULL;   //     x1: solution of the system

        double *wr      = NULL;  
        double *wi      = NULL;  
        double *vr      = NULL;  
        double *work    = NULL;  

        /* SMALL-EXAMPLE: 5x5 MATRIX WITH 3 RHS */
        // NOTE: matrices are stored columnwise, i.e., the first n elementrs of
        // the array 'a' are the elements of the first column of the matrix a.
        // a = {a_11, a_21, ..., a_n1, a_12, a_22, ..., a_n2, a_31, a_32, ...

        double aa[N*N] = {
           6.80, -0.11,  0.66,  1.97,  1.23,
          -0.11,  3.30,  0.36, -0.44,  1.08,
          -0.66,  0.36,  5.70,  0.27,  1.04,
           1.97, -0.44,  0.27,  7.17,  0.14,
           1.23,  1.08,  1.04,  0.14,  6.87
        };
        double aa_copy[N*N] = {
           6.80, -0.11,  0.66,  1.97,  1.23,
          -0.11,  3.30,  0.36, -0.44,  1.08,
          -0.66,  0.36,  5.70,  0.27,  1.04,
           1.97, -0.44,  0.27,  7.17,  0.14,
           1.23,  1.08,  1.04,  0.14,  6.87
        };
        double bb[N*R] = {
            4.02,  6.19, -8.22, -7.57, -3.03,
           -1.56,  4.00, -8.67,  1.75,  2.86,
            9.81, -4.09, -4.57, -8.61,  8.99
        };
        double bb_copy[N*R] = {
            4.02,  6.19, -8.22, -7.57, -3.03,
            -1.56,  4.00, -8.67,  1.75,  2.86,
            9.81, -4.09, -4.57, -8.61,  8.99
        };
        double resres[N] = {
            0.00, 0.00, 0.00, 0.00, 0.00
        };
        double xx[N] = {
            0.00, 0.00, 0.00, 0.00, 0.00
        };
        double xx1[N] = {
            0.00, 0.00, 0.00, 0.00, 0.00
        };
        /* SMALL-EXAMPLE: END */

        int i, j, k, ndim;
        double dummy;

#ifdef _OPENMP
#   pragma omp parallel
       { 
#        pragma omp single 
         printf("OpenMP-parallel with %1d threads\n", omp_get_num_threads());
        } /* end omp parallel */
#endif
    
        for (ndim = 0; ndim < 10; ndim++)
        {
            n    = _n[ndim];
            nrhs = _rhs[ndim];
            lda  = n;
            ldb  = n;


            /* Generate matrices */

            mem_alloc(ipiv,   n,      int);
            mem_alloc(a,      n*n,    double); 
            mem_alloc(a_copy, n*n,    double);  
            mem_alloc(b,      n*nrhs, double);  
            mem_alloc(b_copy, n*nrhs, double);  
            mem_alloc(x,      n,      double);  
            mem_alloc(x1,     n,      double);  
            mem_alloc(res,    n,      double);


            if (ndim == 0)
            {
                // This is a small example.
                
                memcpy(a,           aa, n*n*sizeof(double));
                memcpy(a_copy, aa_copy, n*n*sizeof(double));
                memcpy(b,           bb, n*nrhs*sizeof(double));
                memcpy(b_copy, bb_copy, n*nrhs*sizeof(double));
                memcpy(res,     resres, n*sizeof(double));
                memcpy(x,           xx, n*sizeof(double));
                memcpy(x1,         xx1, n*sizeof(double));

            }
            else
            {
                                
                for (i = 0; i < n; i++) 
                    for (j = 0; j < n; j++) 
	                a[........] =  -1;
            
                for (j = 0; j < n; j++) 
	            a[........] =  1.1*n;
 
                for (i = 0; i < n; i++) 
                    for (j = 0; j < n; j++) 
	                a_copy [........] = a[........];

                for (j = 0; j < n; j++) 
	            b[j] =  1;

                for (j = 0; j < n; j++) 
	            b_copy[j] = b[j];
            }           

            /* Executable statements */
            printf( "\n\n\ndgesv (column-major, high-level), #of equations          %d \n", n );

# ifdef _OPENMP
            wt1=omp_get_wtime();
# endif
            /* Solve the equations A*X = B */
            /* dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info); */
            dgesv_( ........) ;
# ifdef _OPENMP
            wt2=omp_get_wtime();
            printf( "wall clock time factorization (omp_get_wtime)      %12.4g sec\n", wt2-wt1 );
# endif

            /* Check for the exact singularity */
            if (info > 0) 
            {
                printf("The diagonal element of the triangular factor of A,\n");
                printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
                printf("the solution could not be computed.\n");
                exit(1);
            }

            /* Print solution */
            if (ndim == 0)
            {

               /* Executable statements */
               printf( "\n\n\ndgesv (column-major, high-level) Example Program Results\n" );
               print_matrix( "Solution", n, nrhs, b, ldb );          
               /* Print details of LU factorization */
               print_matrix( "Details of LU factorization", n, n, a, lda );
               /* Print pivot indices */
               print_int_vector( "Pivot indices", n, ipiv );
               printf( "\n");
            }
            /* Print residual for each right-hand side */
            for (i = 0; i < nrhs; i++) 
            {
                /* Auxiliary variables */
                double alpha    =  1.0;
                double zero     =  0.0;
                double beta     = -1.0;
                double residual =  0.0;
                double normb    =  0.0;
                int    one_inc  =  1;

                normb = 0.0;
                for (j = 0; j < n; j++) 
                {
                    res[j] =  b_copy[i*n + j];                      // res contains b
                    normb +=  b_copy[i*n + j] * b_copy[i*n + j];
                }
                normb = sqrt(normb);

                dgemv_(....);

                for (j = 0; j < n; j++) 
                    residual += res[j]*res[j];      // residual = sum_{i=0}^{n} res[j]^2

                for (j = 0; j < n; j++) 
                    normb += b[i*n]*res[j];          

                printf( "Norm of residual: ||A x - b ||                     %17.5e \n", residual);
                printf( "Norm of residual residual: ||A x - b || // || b || %17.5e \n", sqrt(residual)/normb);
	    } 

       /*  
           Exercise 2 - Eigenvalue Solver 
       */
 
               for (i = 0; i < n; i++) 
                    for (j = 0; j < n; j++) 
	                a [i + j*n] = a_copy[i + j*n];


            /* Executable statements */
            printf( "dgeev (column-major, high-level), #of equations          %d \n", n );

# ifdef _OPENMP
            wt1=omp_get_wtime();
# endif
            /* Solve the equations A*V = W*V 
  
            Arguments of DGEEV
 	    ===================   

	    JOBVL   (input) CHARACTER*1   
	            = 'N': left eigenvectors of A are not computed;   
        	    = 'V': left eigenvectors of A are computed.   

  	    JOBVR   (input) CHARACTER*1   
   	            = 'N': right eigenvectors of A are not computed;   
   	            = 'V': right eigenvectors of A are computed.   

   	    N       (input) INTEGER   
   	            The order of the matrix A. N >= 0.   

      	    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
     	            On entry, the N-by-N matrix A.   
     	            On exit, A has been overwritten.   

          LDA       (input) INTEGER   
       	            The leading dimension of the array A.  LDA >= max(1,N).   

       	   WR       (output) DOUBLE PRECISION array, dimension (N)   
   	   WI       (output) DOUBLE PRECISION array, dimension (N)   
            	    WR and WI contain the real and imaginary parts,   
                    respectively, of the computed eigenvalues.  Complex   
                    conjugate pairs of eigenvalues appear consecutively   
                    with the eigenvalue having the positive imaginary part   
                    first.   

   	 VL         (output) DOUBLE PRECISION array, dimension (LDVL,N)   
                    If JOBVL = 'V', the left eigenvectors u(j) are stored one   
                    after another in the columns of VL, in the same order   
                    as their eigenvalues.   
                    If JOBVL = 'N', VL is not referenced.   
                    If the j-th eigenvalue is real, then u(j) = VL(:,j),   
                    the j-th column of VL.   
                    If the j-th and (j+1)-st eigenvalues form a complex   
                    conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and   
                    u(j+1) = VL(:,j) - i*VL(:,j+1).   

    	LDVL        (input) INTEGER   
                    The leading dimension of the array VL.  LDVL >= 1; if   
                    JOBVL = 'V', LDVL >= N.   

        VR          (output) DOUBLE PRECISION array, dimension (LDVR,N)   
                    If JOBVR = 'V', the right eigenvectors v(j) are stored one   
                    after another in the columns of VR, in the same order   
                    as their eigenvalues.   
                    If JOBVR = 'N', VR is not referenced.   
                    If the j-th eigenvalue is real, then v(j) = VR(:,j),   
                    the j-th column of VR.   
                    If the j-th and (j+1)-st eigenvalues form a complex   
                    conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and   
                    v(j+1) = VR(:,j) - i*VR(:,j+1).   

   	LDVR        (input) INTEGER   
                    The leading dimension of the array VR.  LDVR >= 1; if   
                    JOBVR = 'V', LDVR >= N.   

    	WORK        (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
                    On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    	LWORK       (input) INTEGER   
                    The dimension of the array WORK.  LWORK >= max(1,3*N), and   
                    if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good   
                    performance, LWORK must generally be larger.   

                    If LWORK = -1, then a workspace query is assumed; the routine   
                    only calculates the optimal size of the WORK array, returns   
                    this value as the first entry of the WORK array, and no error   
                    message related to LWORK is issued by XERBLA.   

    	INFO        (output) INTEGER   
                    = 0:  successful exit   
                    < 0:  if INFO = -i, the i-th argument had an illegal value.   
                    > 0:  if INFO = i, the QR algorithm failed to compute all the   
                    eigenvalues, and no eigenvectors have been computed;   
                    elements i+1:N of WR and WI contain eigenvalues which   
                    have converged.   

                    CALL DGEEV('No left vectors','Vectors (right)',N,A,LDA,WR,WI, DUMMY,1,VR,LDVR,WORK,LWORK,INFO) */

            mem_alloc(wr,      n,    double);
            mem_alloc(wi,      n,    double);
            mem_alloc(vr,   n*n,    double);
            lwork = 5*n;
            mem_alloc(work,  5*n,    double);
   
            j = 1;
            ldvr = n;
            dgeev_(  ........  );
# ifdef _OPENMP
            wt2=omp_get_wtime();
            printf( "wall clock time Eigenvalue (omp_get_wtime)        %12.4g sec\n", wt2-wt1 );
# endif  


            if (ndim == 0)
            {
              /* Executable statements */
               printf( "\n\n\n dgeev (column-major, high-level) Example Program Results");
               print_matrix( "Matrix", n, n, a_copy, lda );
               print_matrix( "Real    Eigenvalues", n, 1 , wr, n );          
               print_matrix( "Complex Eigenvalues", n, 1 , wi, n );          
               print_matrix( "Eigenvectors", n, n, vr, lda );
               printf( "\n");

            }

            mem_free(    wr);  
            mem_free(    wi);  
            mem_free(    vr);  
            mem_free(  work);  





            /* Check for the exact singularity */
            if (info > 0) 
            {
                printf("The diagonal element of the triangular factor of A,\n");
                printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
                printf("the solution could not be computed.\n");
                exit(1);
            }


 
            /*  
              Release memory. 
            */         
            mem_free(  ipiv);  
            mem_free(     a);  
            mem_free(a_copy);  
            mem_free(     b);  
            mem_free(b_copy);  
            mem_free(     x);  
            mem_free(    x1);  
            mem_free(   res);  

        }

        return 0;
} 
/* End of LAPACKE_dgesv Example */


/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, double* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, int n, int* a ) {
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
}

/* Auxiliary routine: printing a vector of double */
void print_double_vector( char* desc, int n, double* a ) {
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %e", a[j] );
	printf( "\n" );
}

