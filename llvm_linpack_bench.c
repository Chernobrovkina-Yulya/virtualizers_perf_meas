# include <float.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>

int main (int argc, char** argv);
double cpu_time ( );
void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy );
double ddot ( int n, double dx[], int incx, double dy[], int incy );
int dgefa ( double a[], int lda, int n, int ipvt[] );
void dgesl ( double a[], int lda, int n, int ipvt[], double b[], int job );
void dscal ( int n, double sa, double x[], int incx );
int idamax ( int n, double dx[], int incx );
double r8_random ( int iseed[4] );
double *r8mat_gen ( int lda, int n );
double r8mat_norm_li ( int lda, int m, int n, double A[] );
void timestamp ( );

/******************************************************************************/

int main (int argc, char** argv)

/******************************************************************************/
/*
  Purpose:

    linpack_bench() runs the LINPACK benchmark program.

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    06 February 2025

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart.
    This version by John Burkardt

  Local:

    int N: the order of the matrix.
*/
{
  double *A;
  double *ALU;
  double A_norm;
  double *b;
  double cray;
  double cray_ratio;
  double eps;
  int i;
  int info;
  int *ipvt;
  int j;
  int job;
  int lda;
  double mflops;
  int n;
  double ops;
  double *r;
  double r_norm;
  double ratio;
  double t;
  double t2;
  double unit;
  double *x;
  double *x_exact;
  double x_norm;

  timestamp ( );
  n = atoi(argv[1]);
  lda = n + 1;

  printf ( "\n" );
  printf ( "linpack_bench()\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  The LINPACK benchmark.\n" );
  printf ( "  Language: C\n" );
  printf ( "  Datatype: Double precision real\n" );
  printf ( "  Matrix order N               = %d\n", n );
  printf ( "  Leading matrix dimension LDA = %d\n", lda );
/*
  Allocate space for arrays.
*/
  A = r8mat_gen ( lda, n );
  ALU = ( double * ) malloc ( lda * n * sizeof ( double ) );
  b = ( double * ) malloc ( n * sizeof ( double ) );
  ipvt = ( int * ) malloc ( n * sizeof ( int ) );
  r = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );
  x_exact = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Copy A.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      ALU[i+j*lda] = A[i+j*lda];
    }
  }

  for ( i = 0; i < n; i++ )
  {
    x_exact[i] = 1.0;
  }

  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + A[i+j*lda] * x_exact[j];
    }
  }

  t = cpu_time ( );
  info = dgefa ( ALU, lda, n, ipvt );
  t = cpu_time ( ) - t;

  if ( info != 0 )
  {
    printf ( "\n" );
    printf ( "linpack_bench(): Fatal error!\n" );
    printf ( "  The matrix A is apparently singular.\n" );
    printf ( "  Abnormal end of execution.\n" );
    return 1;
  }

  job = 0;
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  t2 = cpu_time ( );
  dgesl ( ALU, lda, n, ipvt, x, job );
  t2 = cpu_time ( ) - t2;

  t = t + t2;
/*
  Compute residual r = b - A*x.
*/
  for ( i = 0; i < n; i++ )
  {
    r[i] = b[i];
    for ( j = 0; j < n; j++ )
    {
      r[i] = r[i] - A[i+j*lda] * x[j];
    }
  }
/*
  Compute infinity norms.
*/
  A_norm = r8mat_norm_li ( lda, n, n, A );

  x_norm = 0.0;
  for ( i = 0; i < n; i++ )
  {
    x_norm = fmax ( x_norm, fabs ( x[i] ) );
  }

  r_norm = 0.0;
  for ( i = 0; i < n; i++ )
  {
    r_norm = fmax ( r_norm, fabs ( r[i] ) );
  }
/*
  Report.
*/
  eps = DBL_EPSILON;
  ratio = r_norm / A_norm / x_norm / eps;

  ops = ( double ) ( 2 * n * n * n ) / 3.0 + 2.0 * ( double ) ( n * n );
  mflops = ops / ( 1000000 * t );
  unit = 2.0 / mflops;

  cray = 0.056;
  cray_ratio = t / cray;

  printf ( "\n" );
  printf ( "  Normalized residual = %14f \n", ratio );
  printf ( "  Residual norm       = %14f \n", r_norm );
  printf ( "  Machine epsilon     = %14f \n", eps );
  printf ( "  First X[]           = %14f \n", x[0] );
  printf ( "  Last X[]            = %14f \n", x[n-1] );
  printf ( "  Time in seconds     = %14f \n", t );
  printf ( "  MegaFLOPS           = %14f \n", mflops );
  printf ( "  Unit                = %14f \n", unit );
  printf ( "  Cray ratio          = %14f \n", cray_ratio );
/*
  Free memory.
*/
  free ( A );
  free ( ALU );
  free ( b );
  free ( ipvt );
  free ( r );
  free ( x );
  free ( x_exact );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "linpack_bench():\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

double cpu_time ( )

/******************************************************************************/
/*
  Purpose:

    cpu_time() returns the current reading on the CPU clock.

  Discussion:

    The CPU time measurements available through this routine are often
    not very accurate.  In some cases, the accuracy is no better than
    a hundredth of a second.  

  Licensing:

    This code is distributed under the MIT license. 

  Modified:

    06 June 2005

  Author:

    John Burkardt

  Output:

    double CPU_TIME, the current reading of the CPU clock, in seconds.
*/
{
  double value;

  value = ( double ) clock ( ) 
        / ( double ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/
/*?
{
    "type": "method",
    "virtualize":
    {
        "key1": {
            "type": "bool",
            "value": true
        }
    }
}
?*/
void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy )

/******************************************************************************/
/*
  Purpose:

    DAXPY computes constant times a vector plus a vector.

  Discussion:

    This routine uses unrolled loops for increments equal to one.

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    30 March 2007

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart.
    This version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539, 
    ACM Transactions on Mathematical Software, 
    Volume 5, Number 3, September 1979, pages 308-323.

  Input:

    int N, the number of elements in DX and DY.

    double DA, the multiplier of DX.

    double DX[*], the first vector.

    int INCX, the increment between successive entries of DX.

    double DY[*], the second vector.

    int INCY, the increment between successive entries of DY.

  Output:

    double DY[*]: replaced by DY[*] + DA * DX[*].

*/
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( da == 0.0 )
  {
    return;
  }
/*
  Code for unequal increments or equal increments
  not equal to 1.
*/
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dy[iy] = dy[iy] + da * dx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
/*
  Code for both increments equal to 1.
*/
  else
  {
    m = n % 4;

    for ( i = 0; i < m; i++ )
    {
      dy[i] = dy[i] + da * dx[i];
    }

    for ( i = m; i < n; i = i + 4 )
    {
      dy[i  ] = dy[i  ] + da * dx[i  ];
      dy[i+1] = dy[i+1] + da * dx[i+1];
      dy[i+2] = dy[i+2] + da * dx[i+2];
      dy[i+3] = dy[i+3] + da * dx[i+3];
    }
  }
  return;
}
/******************************************************************************/
/*?
{
    "type": "method",
    "virtualize":
    {
        "key1": {
            "type": "bool",
            "value": true
        }
    }
}
?*/
double ddot ( int n, double dx[], int incx, double dy[], int incy )

/******************************************************************************/
/*
  Purpose:

    ddot() forms the dot product of two vectors.

  Discussion:

    This routine uses unrolled loops for increments equal to one.

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    30 March 2007

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart.
    This version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539, 
    ACM Transactions on Mathematical Software, 
    Volume 5, Number 3, September 1979, pages 308-323.

  Input:

    int N, the number of entries in the vectors.

    double DX[*], the first vector.

    int INCX, the increment between successive entries in DX.

    double DY[*], the second vector.

    int INCY, the increment between successive entries in DY.

  Output:

    double DDOT, the sum of the product of the corresponding
    entries of DX and DY.
*/
{
  double dtemp;
  int i;
  int ix;
  int iy;
  int m;

  dtemp = 0.0;

  if ( n <= 0 )
  {
    return dtemp;
  }
/*
  Code for unequal increments or equal increments
  not equal to 1.
*/
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dtemp = dtemp + dx[ix] * dy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
/*
  Code for both increments equal to 1.
*/
  else
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      dtemp = dtemp + dx[i] * dy[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      dtemp = dtemp + dx[i  ] * dy[i  ] 
                    + dx[i+1] * dy[i+1] 
                    + dx[i+2] * dy[i+2] 
                    + dx[i+3] * dy[i+3] 
                    + dx[i+4] * dy[i+4];
    }
  }
  return dtemp;
}
/******************************************************************************/

int dgefa ( double a[], int lda, int n, int ipvt[] )

/******************************************************************************/
/*
  Purpose:

    dgefa() factors a real general matrix.

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    16 May 2005

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart.
    This version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
    LINPACK User's Guide,
    SIAM, (Society for Industrial and Applied Mathematics),
    3600 University City Science Center,
    Philadelphia, PA, 19104-2688.
    ISBN 0-89871-172-X

  Input:

    double A[LDA*N]: the matrix to be factored.
 
    int LDA, the leading dimension of A.

    int N, the order of the matrix A.

  Output:

    double A[LDA*N]: a compressed set of data defining an upper triangular 
    matrix and the multipliers used to obtain it.  The factorization can 
    be written A=L*U, where L is a product of
    permutation and unit lower triangular matrices, and U is upper triangular.

    int IPVT[N], the pivot indices.

    int DGEFA, singularity indicator.
    0, normal value.
    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
    but it does indicate that DGESL or DGEDI will divide by zero if called.
    Use RCOND in DGECO for a reliable indication of singularity.
*/
{
  int info;
  int j;
  int k;
  int l;
  double t;
/*
  Gaussian elimination with partial pivoting.
*/
  info = 0;

  for ( k = 1; k <= n-1; k++ )
  {
/*
  Find L = pivot index.
*/
    l = idamax ( n-k+1, a+(k-1)+(k-1)*lda, 1 ) + k - 1;
    ipvt[k-1] = l;
/*
  Zero pivot implies this column already triangularized.
*/
    if ( a[l-1+(k-1)*lda] == 0.0 )
    {
      info = k;
      continue;
    }
/*
  Interchange if necessary.
*/
    if ( l != k )
    {
      t = a[l-1+(k-1)*lda];
      a[l-1+(k-1)*lda] = a[k-1+(k-1)*lda];
      a[k-1+(k-1)*lda] = t;
    }
/*
  Compute multipliers.
*/
    t = -1.0 / a[k-1+(k-1)*lda];

    dscal ( n-k, t, a+k+(k-1)*lda, 1 );
/*
  Row elimination with column indexing.
*/
    for ( j = k+1; j <= n; j++ )
    {
      t = a[l-1+(j-1)*lda];
      if ( l != k )
      {
        a[l-1+(j-1)*lda] = a[k-1+(j-1)*lda];
        a[k-1+(j-1)*lda] = t;
      }
      daxpy ( n-k, t, a+k+(k-1)*lda, 1, a+k+(j-1)*lda, 1 );
    }

  }

  ipvt[n-1] = n;

  if ( a[n-1+(n-1)*lda] == 0.0 )
  {
    info = n;
  }

  return info;
}
/******************************************************************************/

void dgesl ( double a[], int lda, int n, int ipvt[], double b[], int job )

/******************************************************************************/
/*
  Purpose:

    dgesl() solves a real general linear system A * X = B.

  Discussion:

    DGESL can solve either of the systems A * X = B or A' * X = B.

    The system matrix must have been factored by DGECO or DGEFA.

    A division by zero will occur if the input factor contains a
    zero on the diagonal.  Technically this indicates singularity
    but it is often caused by improper arguments or improper
    setting of LDA.  It will not occur if the subroutines are
    called correctly and if DGECO has set 0.0 < RCOND
    or DGEFA has set INFO == 0.

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    16 May 2005

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart.
    This version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
    LINPACK User's Guide,
    SIAM, (Society for Industrial and Applied Mathematics),
    3600 University City Science Center,
    Philadelphia, PA, 19104-2688.
    ISBN 0-89871-172-X

  Input:

    double A[LDA*N], the output from DGECO or DGEFA.

    int LDA, the leading dimension of A.

    int N, the order of the matrix A.

    int IPVT[N], the pivot vector from DGECO or DGEFA.

    double B[N]: the right hand side vector.

    int JOB.
    0, solve A * X = B;
    nonzero, solve A' * X = B.

  Output:

    double B[N]: the solution vector.
*/
{
  int k;
  int l;
  double t;
/*
  Solve A * X = B.
*/
  if ( job == 0 )
  {
    for ( k = 1; k <= n-1; k++ )
    {
      l = ipvt[k-1];
      t = b[l-1];

      if ( l != k )
      {
        b[l-1] = b[k-1];
        b[k-1] = t;
      }

      daxpy ( n-k, t, a+k+(k-1)*lda, 1, b+k, 1 );

    }

    for ( k = n; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
      t = -b[k-1];
      daxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
    }
  }
/*
  Solve A' * X = B.
*/
  else
  {
    for ( k = 1; k <= n; k++ )
    {
      t = ddot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
      b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
    }

    for ( k = n-1; 1 <= k; k-- )
    {
      b[k-1] = b[k-1] + ddot ( n-k, a+k+(k-1)*lda, 1, b+k, 1 );
      l = ipvt[k-1];

      if ( l != k )
      {
        t = b[l-1];
        b[l-1] = b[k-1];
        b[k-1] = t;
      }
    }
  }
  return;
}
/******************************************************************************/
/*?
{
    "type": "method",
    "virtualize":
    {
        "key1": {
            "type": "bool",
            "value": true
        }
    }
}
?*/
void dscal ( int n, double sa, double x[], int incx )

/******************************************************************************/
/*
  Purpose:

    dscal() scales a vector by a constant.

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    30 March 2007

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart.
    This version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Input:

    int N, the number of entries in the vector.

    double SA, the multiplier.

    double X[*], the vector to be scaled.

    int INCX, the increment between successive entries of X.

  Output:

    double X[*]: the scaled vector.
*/
{
  int i;
  int ix;
  int m;

  if ( n <= 0 )
  {
  }
  else if ( incx == 1 )
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for ( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }
  }
  return;
}
/******************************************************************************/
/*?
{
    "type": "method",
    "virtualize":
    {
        "key1": {
            "type": "bool",
            "value": true
        }
    }
}
?*/
int idamax ( int n, double dx[], int incx )

/******************************************************************************/
/*
  Purpose:

    idamax() finds the index of the vector element of maximum absolute value.

  Discussion:

    This index is a 1-based index, not a 0-based index!

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    30 March 2007

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart.
    This version by John Burkardt

  Reference:

    Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979.

    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    Basic Linear Algebra Subprograms for Fortran Usage,
    Algorithm 539,
    ACM Transactions on Mathematical Software,
    Volume 5, Number 3, September 1979, pages 308-323.

  Input:

    int N, the number of entries in the vector.

    double X[*], the vector to be examined.

    int INCX, the increment between successive entries of SX.

  Output:

    int IDAMAX, the index of the element of maximum absolute value.
*/
{
  double dmax;
  int i;
  int ix;
  int value;

  value = 0;

  if ( n < 1 || incx <= 0 )
  {
    return value;
  }

  value = 1;

  if ( n == 1 )
  {
    return value;
  }

  if ( incx == 1 )
  {
    dmax = fabs ( dx[0] );

    for ( i = 1; i < n; i++ )
    {
      if ( dmax < fabs ( dx[i] ) )
      {
        value = i + 1;
        dmax = fabs ( dx[i] );
      }
    }
  }
  else
  {
    ix = 0;
    dmax = fabs ( dx[0] );
    ix = ix + incx;

    for ( i = 1; i < n; i++ )
    {
      if ( dmax < fabs ( dx[ix] ) )
      {
        value = i + 1;
        dmax = fabs ( dx[ix] );
      }
      ix = ix + incx;
    }
  }

  return value;
}
/******************************************************************************/

double r8_random ( int iseed[4] )

/******************************************************************************/
/*
  Purpose:

    r8_random() returns a uniformly distributed random number between 0 and 1.

  Discussion:

    This routine uses a multiplicative congruential method with modulus
    2**48 and multiplier 33952834046453 (see G.S.Fishman,
    'Multiplicative congruential random number generators with modulus
    2**b: an exhaustive analysis for b = 32 and a partial analysis for
    b = 48', Math. Comp. 189, pp 331-344, 1990).

    48-bit integers are stored in 4 integer array elements with 12 bits
    per element. Hence the routine is portable across machines with
    integers of 32 bits or more.

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    05 February 2025

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart.
    This version by John Burkardt

  Input:

    int ISEED(4): the seeds of the random number generator; the array
    elements must be between 0 and 4095, and ISEED(4) must be odd.

  Output:

    int ISEED(4): the updated seeds.

    double R8_RANDOM, the next pseudorandom number.
*/
{
  int ipw2 = 4096;
  int it1;
  int it2;
  int it3;
  int it4;
  int m1 = 494;
  int m2 = 322;
  int m3 = 2508;
  int m4 = 2549;
  double r = 1.0 / 4096.0;
  double value;
/*
  Multiply the seed by the multiplier modulo 2**48.
*/
  it4 = iseed[3] * m4;
  it3 = it4 / ipw2;
  it4 = it4 - ipw2 * it3;
  it3 = it3 + iseed[2] * m4 + iseed[3] * m3;
  it2 = it3 / ipw2;
  it3 = it3 - ipw2 * it2;
  it2 = it2 + iseed[1] * m4 + iseed[2] * m3 + iseed[3] * m2;
  it1 = it2 / ipw2;
  it2 = it2 - ipw2 * it1;
  it1 = it1 + iseed[0] * m4 + iseed[1] * m3 + iseed[2] * m2 + iseed[3] * m1;
  it1 = ( it1 % ipw2 );
/*
  Return updated seed
*/
  iseed[0] = it1;
  iseed[1] = it2;
  iseed[2] = it3;
  iseed[3] = it4;
/*
  Convert 48-bit integer to a real number in the interval (0,1)
*/
  value = 
      r * ( ( double ) ( it1 ) 
    + r * ( ( double ) ( it2 ) 
    + r * ( ( double ) ( it3 ) 
    + r * ( ( double ) ( it4 ) ) ) ) );

  return value;
}
/******************************************************************************/

double *r8mat_gen ( int lda, int n )

/******************************************************************************/
/*
  Purpose:

    r8mat_gen() generates a random R8MAT.

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    06 June 2005

  Author:

    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, Pete Stewart.
    This version by John Burkardt

  Input:

    integer LDA, the leading dimension of the matrix.

    integer N, the order of the matrix.

  Output:

    double R8MAT_GEN[LDA*N], the N by N matrix.
*/
{
  double *a;
  int i;
  int init[4] = { 1, 2, 3, 1325 };
  int j;

  a = ( double * ) malloc ( lda * n * sizeof ( double ) );

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= n; i++ )
    {
      a[i-1+(j-1)*lda] = r8_random ( init ) - 0.5;
    }
  }

  return a;
}
/******************************************************************************/

double r8mat_norm_li ( int lda, int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    r8mat_norm_li() returns the matrix L-oo norm of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The matrix L-oo norm is defined as:

      R8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).

    The matrix L-oo norm is derived from the vector L-oo norm,
    and satisifies:

      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).

  Licensing:

    This code is distributed under the MIT license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Input:

    int M, the number of rows in A.

    int N, the number of columns in A.

    double A[M*N], the matrix whose L-oo norm is desired.

  Output:

    double R8MAT_NORM_LI, the L-oo norm of A.
*/
{
  int i;
  int j;
  double row_sum;
  double value;

  value = 0.0;

  for ( i = 0; i < m; i++ )
  {
    row_sum = 0.0;
    for ( j = 0; j < n; j++ )
    {
      row_sum = row_sum + fabs ( a[i+j*lda] );
    }
    value = fmax ( value, row_sum );
  }
  return value;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    timestamp() prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the MIT license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
