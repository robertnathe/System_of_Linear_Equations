#include <stdio.h>
#include <math.h>		// required for fabs()
#define N 3
/******************************************************************************

                              Online C++ Compiler.
               Code, Compile, Run and Debug C++ program online.
Write your code in this editor and press "Run" button to compile and execute it.

*******************************************************************************/
////////////////////////////////////////////////////////////////////////////////
// File: LU_Solution.cpp                                                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// File: lower_triangular.c                                                   //
// Routines:                                                                  //
//    Lower_Triangular_Solve                                                  //
//    Lower_Triangular_Inverse                                                //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Lower_Triangular_Solve(double *L, double *B, double x[], int n)       //
//                                                                            //
//  Description:                                                              //
//     This routine solves the linear equation Lx = B, where L is an n x n    //
//     lower triangular matrix.  (The superdiagonal part of the matrix is     //
//     not addressed.)                                                        //
//     The algorithm follows:                                                 //
//                      x[0] = B[0]/L[0][0], and                              //
//     x[i] = [B[i] - (L[i][0] * x[0]  + ... + L[i][i-1] * x[i-1])] / L[i][i],//
//     for i = 1, ..., n-1.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *L   Pointer to the first element of the lower triangular       //
//                 matrix.                                                    //
//     double *B   Pointer to the column vector, (n x 1) matrix, B.           //
//     double *x   Pointer to the column vector, (n x 1) matrix, x.           //
//     int     n   The number of rows or columns of the matrix L.             //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix L is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], x[N];                                            //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     err = Lower_Triangular_Solve(&A[0][0], B, x, n);                       //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else printf(" The solution is \n");                                    //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Lower_Triangular_Solve (double *L, double B[], double x[], int n)
{
  int i, k;

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.

  for (k = 0; k < n; L += n, k++)
    {
      if (*(L + k) == 0.0)
	return -1;		// The matrix L is singular
      x[k] = B[k];
      for (i = 0; i < k; i++)
	x[k] -= x[i] * *(L + i);
      x[k] /= *(L + k);
    }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//  int Lower_Triangular_Inverse(double *L,  int n)                           //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the inverse of the lower triangular matrix L.  //
//     The superdiagonal part of the matrix is not addressed.                 //
//     The algorithm follows:                                                 //
//        Let M be the inverse of L, then L M = I,                            //
//     M[i][i] = 1.0 / L[i][i] for i = 0, ..., n-1, and                       //
//     M[i][j] = -[(L[i][j] M[j][j] + ... + L[i][i-1] M[i-1][j])] / L[i][i],  //
//     for i = 1, ..., n-1, j = 0, ..., i - 1.                                //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *L   On input, the pointer to the first element of the matrix   //
//                 whose lower triangular elements form the matrix which is   //
//                 to be inverted. On output, the lower triangular part is    //
//                 replaced by the inverse.  The superdiagonal elements are   //
//                 not modified.                                              //
//                 its inverse.                                               //
//     int     n   The number of rows and/or columns of the matrix L.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix L is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double L[N][N];                                                        //
//                                                                            //
//     (your code to create the matrix L)                                     //
//     err = Lower_Triangular_Inverse(&L[0][0], N);                           //
//     if (err < 0) printf(" Matrix L is singular\n");                        //
//     else {                                                                 //
//        printf(" The inverse is \n");                                       //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Lower_Triangular_Inverse (double *L, int n)
{
  int i, j, k;
  double *p_i, *p_j, *p_k;
  double sum;

//         Invert the diagonal elements of the lower triangular matrix L.

  for (k = 0, p_k = L; k < n; p_k += (n + 1), k++)
    {
      if (*p_k == 0.0)
	return -1;
      else
	*p_k = 1.0 / *p_k;
    }

//         Invert the remaining lower triangular matrix L row by row.

  for (i = 1, p_i = L + n; i < n; i++, p_i += n)
    {
      for (j = 0, p_j = L; j < i; p_j += n, j++)
	{
	  sum = 0.0;
	  for (k = j, p_k = p_j; k < i; k++, p_k += n)
	    sum += *(p_i + k) * *(p_k + j);
	  *(p_i + j) = -*(p_i + i) * sum;
	}
    }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// File: lower_triangular_lt.c                                                //
// Routines:                                                                  //
//    Lower_Triangular_Solve_lt                                               //
//    Lower_Triangular_Inverse_lt                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Lower_Triangular_Solve_lt(double *L, double *B, double x[], int n)    //
//                                                                            //
//  Description:                                                              //
//     This routine solves the linear equation Lx = B, where L is an n x n    //
//     lower triangular matrix.                                               //
//                                                                            //
//     The input matrix L is stored in lower triangular form, i.e. the        //
//     elements of L are stored sequentially starting with the first row which//
//     contains one element, followed by the second row which consists of two //
//     elements, and so on until the last row which consists of n elements.   //
//                                                                            //
//     The algorithm follows:                                                 //
//                      x[0] = B[0]/L[0][0], and                              //
//     x[i] = [B[i] - (L[i][0] * x[0]  + ... + L[i][i-1] * x[i-1])] / L[i][i],//
//     for i = 1, ..., n-1.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *L   Pointer to the first element of the lower triangular       //
//                 matrix.                                                    //
//     double *B   Pointer to the column vector, (n x 1) matrix, B.           //
//     double *x   Pointer to the column vector, (n x 1) matrix, x.           //
//     int     n   The number of rows or columns of the matrix L.             //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix L is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N * (N + 1) / 2], B[N], x[N];                                 //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     err = Lower_Triangular_Solve_lt(A, B, x, n);                           //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else printf(" The solution is \n");                                    //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Lower_Triangular_Solve_lt (double *L, double B[], double x[], int n)
{
  int i, k;

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix stored in lower triangular form.

  for (k = 0; k < n; L += ++k)
    {
      if (*(L + k) == 0.0)
	return -1;		// The matrix L is singular
      x[k] = B[k];
      for (i = 0; i < k; i++)
	x[k] -= x[i] * *(L + i);
      x[k] /= *(L + k);
    }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  int Lower_Triangular_Inverse_lt(double *L,  int n)                        //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the inverse of the lower triangular matrix L.  //
//                                                                            //
//     The input matrix L is stored in lower triangular form, i.e. the        //
//     elements of L are stored sequentially starting with the first row which//
//     contains one element, followed by the second row which consists of two //
//     elements, and so on until the last row which consists of n elements.   //
//                                                                            //
//     The algorithm follows:                                                 //
//        Let M be the inverse of L, then L M = I,                            //
//     M[i][i] = 1.0 / L[i][i] for i = 0, ..., n-1, and                       //
//     M[i][j] = -[(L[i][j] M[j][j] + ... + L[i][i-1] M[i-1][j])] / L[i][i],  //
//     for i = 1, ..., n-1, j = 0, ..., i - 1.                                //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *L   On input, the pointer to the first element of the matrix   //
//                 whose lower triangular elements form the matrix which is   //
//                 to be inverted. On output, the lower triangular part is    //
//                 replaced by the inverse.  The superdiagonal elements are   //
//                 not modified.                                              //
//     int     n   The number of rows and/or columns of the matrix L.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix L is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double L[N * (N - 1) / 2];                                             //
//                                                                            //
//     (your code to create the matrix L)                                     //
//     err = Lower_Triangular_Inverse_lt(L, N);                               //
//     if (err < 0) printf(" Matrix L is singular\n");                        //
//     else {                                                                 //
//        printf(" The inverse is \n");                                       //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Lower_Triangular_Inverse_lt (double *L, int n)
{
  int i, j, k;
  double *p_i, *p_j, *p_k;
  double sum;

//         Invert the diagonal elements of the lower triangular matrix L.

  for (k = 0, p_k = L; k < n; p_k += ++k + 1)
    {
      if (*p_k == 0.0)
	return -1;
      else
	*p_k = 1.0 / *p_k;
    }

//         Invert the remaining lower triangular matrix L row by row.

  for (i = 1, p_i = L + 1; i < n; p_i += ++i)
    {
      for (j = 0, p_j = L; j < i; p_j += ++j)
	{
	  sum = 0.0;
	  for (k = j, p_k = p_j; k < i; p_k += ++k)
	    sum += *(p_i + k) * *(p_k + j);
	  *(p_i + j) = -*(p_i + i) * sum;
	}
    }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// File: unit_lower_triangular.c                                              //
// Routines:                                                                  //
//    Unit_Lower_Triangular_Solve                                             //
//    Unit_Lower_Triangular_Inverse                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Unit_Lower_Triangular_Solve(double *L, double *B, double x[], int n) //
//                                                                            //
//  Description:                                                              //
//     This routine solves the linear equation Lx = B, where L is an n x n    //
//     unit lower triangular matrix.  (Only the subdiagonal part of the matrix//
//     is addressed.)  The diagonal is assumed to consist of 1's and is not   //
//     addressed.                                                             //
//     The algorithm follows:                                                 //
//                          x[0] = B[0], and                                  //
//            x[i] = B[i] - (L[i][0] * x[0]  + ... + L[i][i-1] * x[i-1]),     //
//     for i = 1, ..., n-1.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *L   Pointer to the first element of the unit lower triangular  //
//                 matrix.                                                    //
//     double *B   Pointer to the column vector, (n x 1) matrix, B.           //
//     double *x   Pointer to the column vector, (n x 1) matrix, x.           //
//     int     n   The number of rows or columns of the matrix L.             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], x[N];                                            //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     Unit_Lower_Triangular_Solve(&A[0][0], B, x, n);                        //
//     printf(" The solution is \n");                                         //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void
Unit_Lower_Triangular_Solve (double *L, double B[], double x[], int n)
{
  int i, k;

//         Solve the linear equation Lx = B for x, where L is a unit lower
//         triangular matrix.

  x[0] = B[0];
  for (k = 1, L += n; k < n; L += n, k++)
    for (i = 0, x[k] = B[k]; i < k; i++)
      x[k] -= x[i] * *(L + i);
}


////////////////////////////////////////////////////////////////////////////////
//  void Unit_Lower_Triangular_Inverse(double *L,  int n)                     //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the inverse of the unit lower triangular       //
//     matrix L.  Only the subdiagonal part of the matrix is addressed.       //
//     The diagonal is assumed to consist of 1's and is not addressed.        //
//     The algorithm follows:                                                 //
//        Let M be the inverse of L, then L M = I,                            //
//          M[i][j] = -( L[i][j] M[j][j] + ... + L[i][i-1] M[i-1][j] ),       //
//     for i = 1, ..., n-1, j = 0, ..., i - 1.                                //
//                                                                            //
//  Arguments:                                                                //
//     double *L   On input, the pointer to the first element of the matrix   //
//                 whose unit lower triangular elements form the matrix which //
//                 is to be inverted. On output, the lower triangular part is //
//                 replaced by the inverse.  The diagonal and superdiagonal   //
//                 elements are not modified.                                 //
//     int     n   The number of rows and/or columns of the matrix L.         //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double L[N][N];                                                        //
//                                                                            //
//     (your code to create the matrix L)                                     //
//     Unit_Lower_Triangular_Inverse(&L[0][0], N);                            //
//     printf(" The inverse is \n");                                          //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Unit_Lower_Triangular_Inverse (double *L, int n)
{
  int i, j, k;
  double *p_i, *p_j, *p_k;

//         Invert the subdiagonal part of the matrix L row by row where
//         the diagonal elements are assumed to be 1.0.

  for (i = 1, p_i = L + n; i < n; i++, p_i += n)
    {
      for (j = 0, p_j = L; j < i; p_j += n, j++)
	{
	  *(p_i + j) = -*(p_i + j);
	  for (k = j + 1, p_k = p_j + n; k < i; k++, p_k += n)
	    *(p_i + j) -= *(p_i + k) * *(p_k + j);
	}
    }
}

////////////////////////////////////////////////////////////////////////////////
// File: unit_lower_triangular_lt.c                                           //
// Routines:                                                                  //
//    Unit_Lower_Triangular_Solve_lt                                          //
//    Unit_Lower_Triangular_Inverse_lt                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Unit_Lower_Triangular_Solve_lt(double *L, double *B, double x[],     //
//                                                                  int n)    //
//                                                                            //
//  Description:                                                              //
//     This routine solves the linear equation Lx = B, where L is an n x n    //
//     unit lower triangular matrix.  The diagonal is assumed to consist of   //
//     1's and is not addressed.                                              //
//                                                                            //
//     The input matrix L is stored in lower triangular form, i.e. the        //
//     elements of L are stored sequentially starting with the first row which//
//     contains one element, followed by the second row which consists of two //
//     elements, and so on until the last row which consists of n elements.   //
//                                                                            //
//     The algorithm follows:                                                 //
//                          x[0] = B[0], and                                  //
//             x[i] = B[i] - (L[i][0] * x[0]  + ... + L[i][i-1] * x[i-1]),    //
//     for i = 1, ..., n-1.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *L   Pointer to the first element of the unit lower triangular  //
//                 matrix.                                                    //
//     double *B   Pointer to the column vector, (n x 1) matrix, B.           //
//     double *x   Pointer to the column vector, (n x 1) matrix, x.           //
//     int     n   The number of rows or columns of the matrix L.             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N * (N + 1) / 2], B[N], x[N];                                 //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     Unit_Lower_Triangular_Solve_lt(A, B, x, n);                            //
//     else printf(" The solution is \n");                                    //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void
Unit_Lower_Triangular_Solve_lt (double *L, double B[], double x[], int n)
{
  int i, k;

//         Solve the linear equation Lx = B for x, where L is a unit lower
//         triangular matrix stored in lower triangular form.

  for (k = 0; k < n; L += ++k)
    {
      x[k] = B[k];
      for (i = 0; i < k; i++)
	x[k] -= x[i] * *(L + i);
    }
}


////////////////////////////////////////////////////////////////////////////////
//  void Unit_Lower_Triangular_Inverse_lt(double *L,  int n)                  //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the inverse of the unit lower triangular       //
//     matrix L.  The diagonal is not addressed.                              //
//                                                                            //
//     The input matrix L is stored in lower triangular form, i.e. the        //
//     elements of L are stored sequentially starting with the first row which//
//     contains one element, followed by the second row which consists of two //
//     elements, and so on until the last row which consists of n elements.   //
//                                                                            //
//     The algorithm follows:                                                 //
//        Let M be the inverse of L, then L M = I,                            //
//          M[i][j] = -( L[i][j] M[j][j] + ... + L[i][i-1] M[i-1][j] ),       //
//     for i = 1, ..., n-1, j = 0, ..., i - 1.                                //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *L   On input, the pointer to the first element of the matrix   //
//                 whose unit lower triangular elements form the matrix which //
//                 is to be inverted. On output, the lower triangular part is //
//                 replaced by the inverse.  The diagonal elements are not    //
//                 modified.                                                  //
//     int     n   The number of rows and/or columns of the matrix L.         //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double L[N*(N+1)/2];                                                   //
//                                                                            //
//     (your code to create the matrix L)                                     //
//     err = Unit_Lower_Triangular_Inverse(L, N);                             //
//     printf(" The inverse is \n");                                          //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void
Unit_Lower_Triangular_Inverse_lt (double *L, int n)
{
  int i, j, k;
  double *p_i, *p_j, *p_k;

//         Invert the subdiagonal part of the matrix L row by row where
//         the diagonal elements are assumed to be 1.0.

  for (i = 1, p_i = L + 1; i < n; p_i += ++i)
    {
      for (j = 0, p_j = L; j < i; p_j += ++j)
	{
	  *(p_i + j) = -*(p_i + j);
	  for (k = j + 1, p_k = p_j + j + 1; k < i; p_k += ++k)
	    *(p_i + j) -= *(p_i + k) * *(p_k + j);
	}
    }
}

////////////////////////////////////////////////////////////////////////////////
// File: unit_upper_triangular.c                                              //
// Routines:                                                                  //
//    Unit_Upper_Triangular_Solve                                             //
//    Unit_Upper_Triangular_Inverse                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Unit_Upper_Triangular_Solve(double *U, double *B, double x[], int n)  //
//                                                                            //
//  Description:                                                              //
//     This routine solves the linear equation Ux = B, where U is an n x n    //
//     unit upper triangular matrix.  (Only the superdiagonal part of the     //
//     matrix is addressed.)  The diagonal is assumed to consist of 1's and   //
//     is not addressed.                                                      //
//     The algorithm follows:                                                 //
//                  x[n-1] = B[n-1], and                                      //
//       x[i] = B[i] - (U[i][i+1] * x[i+1]  + ... + U[i][n-1] * x[n-1]),      //
//     for i = n-2, ..., 0.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *U   Pointer to the first element of the upper triangular       //
//                 matrix.                                                    //
//     double *B   Pointer to the column vector, (n x 1) matrix, B.           //
//     double *x   Pointer to the column vector, (n x 1) matrix, x.           //
//     int     n   The number of rows or columns of the matrix U.             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], x[N];                                            //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     Unit_Upper_Triangular_Solve(&A[0][0], B, x, n);                        //
//     printf(" The solution is \n");                                         //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void
Unit_Upper_Triangular_Solve (double *U, double B[], double x[], int n)
{
  int i, k;

//         Solve the linear equation Ux = B for x, where U is an upper
//         triangular matrix.
  x[n - 1] = B[n - 1];
  for (k = n - 2, U += n * (n - 2); k >= 0; U -= n, k--)
    {
      x[k] = B[k];
      for (i = k + 1; i < n; i++)
	x[k] -= x[i] * *(U + i);
    }
}


////////////////////////////////////////////////////////////////////////////////
//  int Unit_Upper_Triangular_Inverse(double *U,  int n)                      //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the inverse of the unit upper triangular matrix//
//     U.  The subdiagonal part of the matrix is not addressed.               //
//     The diagonal is assumed to consist of 1's and is not addressed.        //
//     The algorithm follows:                                                 //
//        Let M be the inverse of U, then U M = I,                            //
//          M[i][j] = -( U[i][i+1] M[i+1][j] + ... + U[i][j] M[j][j] ),       //
//     for i = n-2, ... , 0,  j = n-1, ..., i+1.                              //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *U   On input, the pointer to the first element of the matrix   //
//                 whose unit upper triangular elements form the matrix which //
//                 is to be inverted. On output, the upper triangular part is //
//                 replaced by the inverse.  The subdiagonal elements are     //
//                 not modified.                                              //
//     int     n   The number of rows and/or columns of the matrix U.         //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double U[N][N];                                                        //
//                                                                            //
//     (your code to create the matrix U)                                     //
//     Unit_Upper_Triangular_Inverse(&U[0][0], N);                            //
//     printf(" The inverse is \n");                                          //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// File: unit_upper_triangular_ut.c                                           //
// Routines:                                                                  //
//    Unit_Upper_Triangular_Solve_ut                                          //
//    Unit_Upper_Triangular_Inverse_ut                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Unit_Upper_Triangular_Solve_ut(double *U, double *B, double x[],     //
//                                                                  int n)    //
//                                                                            //
//  Description:                                                              //
//     This routine solves the linear equation Ux = B, where U is an n x n    //
//     unit upper triangular matrix.  The diagonal is assumed to consist of   //
//     1's and is not addressed.                                              //
//                                                                            //
//     The input matrix U is stored in upper triangular form, i.e. the        //
//     elements of U are stored sequentially starting with the first row which//
//     contains n elements, followed by the second row which consists of n-1  //
//     elements, and so on until the last row which consists of 1 element.    //
//                                                                            //
//     The algorithm follows:                                                 //
//                       x[n-1] = B[n-1], and                                 //
//       x[i] = [B[i] - (U[i][i+1] * x[i+1]  + ... + U[i][n-1] * x[n-1])],    //
//     for i = n-2, ..., 0.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *U   Pointer to the first element of the unit upper triangular  //
//                 matrix.                                                    //
//     double *B   Pointer to the column vector, (n x 1) matrix, B.           //
//     double *x   Pointer to the column vector, (n x 1) matrix, x.           //
//     int     n   The number of rows or columns of the matrix U.             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N * (N + 1) / 2], B[N], x[N];                                 //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     Unit_Upper_Triangular_Solve_ut(A, B, x, n);                            //
//     printf(" The solution is \n");                                         //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void
Unit_Upper_Triangular_Solve_ut (double *U, double B[], double x[], int n)
{
  int i, k;

//         Solve the linear equation Ux = B for x, where U is a unit upper
//         triangular matrix stored in upper triangular form.

  for (k = n - 1, U += ((n * (n + 1)) >> 1) - 1; k >= 0; U -= n - --k)
    {
      x[k] = B[k];
      for (i = k + 1; i < n; i++)
	x[k] -= x[i] * *(U + i - k);
    }
}


////////////////////////////////////////////////////////////////////////////////
//  void Unit_Upper_Triangular_Inverse_ut(double *U,  int n)                  //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the inverse of the unit upper triangular       //
//      matrix U.  The diagonal is not addressed.                             //
//                                                                            //
//     The matrix U is stored in upper triangular form, i.e. the elements of  //
//     are stored sequentially starting with the first row which contains n   //
//     elements, followed by the second row which consists of n-1 elements,   //
//     and so on until the last row which consists of 1 element.              //
//                                                                            //
//     The algorithm follows:                                                 //
//        Let M be the inverse of U, then U M = I,                            //
//          M[i][j] = -( U[i][i+1] M[i+1][j] + ... + U[i][j] M[j][j] ),       //
//     for i = n-2, ... , 0,  j = n-1, ..., i+1.                              //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *U   On input, the pointer to the first element of the matrix   //
//                 whose upper triangular elements form the matrix which is   //
//                 to be inverted. On output, the upper triangular part is    //
//                 replaced by the inverse.  The diagonal elements are not    //
//                 modified.                                                  //
//     int     n   The number of rows and/or columns of the matrix U.         //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double U[N * (N + 1) / 2];                                             //
//                                                                            //
//     (your code to create the matrix U)                                     //
//     Unit_Upper_Triangular_Inverse_ut(U, N);                                //
//     printf(" The inverse is \n");                                          //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void
Unit_Upper_Triangular_Inverse_ut (double *U, int n)
{
  int i, j, k;
  double *p_i, *p_k;
  double sum;

//         Invert the superdiagonal part of the matrix U row by row where
//         the diagonal elements are assumed to be 1.0.

  for (i = n - 1, p_i = U + ((n * (n + 1)) >> 1) - 1; i >= 0; p_i -= n - --i)
    {
      for (j = n - 1; j > i; j--)
	{
	  sum = 0.0;
	  for (k = i + 1, p_k = p_i + n - i; k < j; p_k += n - k++)
	    {
	      sum += *(p_i + k - i) * *(p_k + j - k);
	    }
	  sum += *(p_i + j - i);
	  *(p_i + j - i) = -sum;
	}
    }
}

////////////////////////////////////////////////////////////////////////////////
// File: upper_triangular.c                                                   //
// Routines:                                                                  //
//    Upper_Triangular_Solve                                                  //
//    Upper_Triangular_Inverse                                                //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Upper_Triangular_Solve(double *U, double *B, double x[], int n)       //
//                                                                            //
//  Description:                                                              //
//     This routine solves the linear equation Ux = B, where U is an n x n    //
//     upper triangular matrix.  (The subdiagonal part of the matrix is       //
//     not addressed.)                                                        //
//     The algorithm follows:                                                 //
//                  x[n-1] = B[n-1]/U[n-1][n-1], and                          //
//     x[i] = [B[i] - (U[i][i+1] * x[i+1]  + ... + U[i][n-1] * x[n-1])]       //
//                                                                 / U[i][i], //
//     for i = n-2, ..., 0.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *U   Pointer to the first element of the upper triangular       //
//                 matrix.                                                    //
//     double *B   Pointer to the column vector, (n x 1) matrix, B.           //
//     double *x   Pointer to the column vector, (n x 1) matrix, x.           //
//     int     n   The number of rows or columns of the matrix U.             //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix U is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], x[N];                                            //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     err = Upper_Triangular_Solve(&A[0][0], B, x, n);                       //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else printf(" The solution is \n");                                    //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Upper_Triangular_Solve (double *U, double B[], double x[], int n)
{
  int i, k;

//         Solve the linear equation Ux = B for x, where U is an upper
//         triangular matrix.

  for (k = n - 1, U += n * (n - 1); k >= 0; U -= n, k--)
    {
      if (*(U + k) == 0.0)
	return -1;		// The matrix U is singular
      x[k] = B[k];
      for (i = k + 1; i < n; i++)
	x[k] -= x[i] * *(U + i);
      x[k] /= *(U + k);
    }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  int Upper_Triangular_Inverse(double *U,  int n)                           //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the inverse of the upper triangular matrix U.  //
//     The subdiagonal part of the matrix is not addressed.                   //
//     The algorithm follows:                                                 //
//        Let M be the inverse of U, then U M = I,                            //
//     M[n-1][n-1] = 1.0 / U[n-1][n-1] and                                    //
//     M[i][j] = -( U[i][i+1] M[i+1][j] + ... + U[i][j] M[j][j] ) / U[i][i],  //
//     for i = n-2, ... , 0,  j = n-1, ..., i+1.                              //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *U   On input, the pointer to the first element of the matrix   //
//                 whose upper triangular elements form the matrix which is   //
//                 to be inverted. On output, the upper triangular part is    //
//                 replaced by the inverse.  The subdiagonal elements are     //
//                 not modified.                                              //
//     int     n   The number of rows and/or columns of the matrix U.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix U is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double U[N][N];                                                        //
//                                                                            //
//     (your code to create the matrix U)                                     //
//     err = Upper_Triangular_Inverse(&U[0][0], N);                           //
//     if (err < 0) printf(" Matrix U is singular\n");                        //
//     else {                                                                 //
//        printf(" The inverse is \n");                                       //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Upper_Triangular_Inverse (double *U, int n)
{
  int i, j, k;
  double *p_i, *p_k;
  double sum;

//         Invert the diagonal elements of the upper triangular matrix U.

  for (k = 0, p_k = U; k < n; p_k += (n + 1), k++)
    {
      if (*p_k == 0.0)
	return -1;
      else
	*p_k = 1.0 / *p_k;
    }

//         Invert the remaining upper triangular matrix U.

  for (i = n - 2, p_i = U + n * (n - 2); i >= 0; p_i -= n, i--)
    {
      for (j = n - 1; j > i; j--)
	{
	  sum = 0.0;
	  for (k = i + 1, p_k = p_i + n; k <= j; p_k += n, k++)
	    {
	      sum += *(p_i + k) * *(p_k + j);
	    }
	  *(p_i + j) = -*(p_i + i) * sum;
	}
    }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// File: upper_triangular_ut.c                                                //
// Routines:                                                                  //
//    Upper_Triangular_Solve_ut                                               //
//    Upper_Triangular_Inverse_ut                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Upper_Triangular_Solve_ut(double *U, double *B, double x[], int n)    //
//                                                                            //
//  Description:                                                              //
//     This routine solves the linear equation Ux = B, where U is an n x n    //
//     upper triangular matrix.                                               //
//                                                                            //
//     The matrix U is stored in upper triangular form, i.e. the elements of  //
//     are stored sequentially starting with the first row which contains n   //
//     elements, followed by the second row which consists of n-1 elements,   //
//     and so on until the last row which consists of 1 element.              //
//                                                                            //
//     The algorithm follows:                                                 //
//                  x[n-1] = B[n-1]/U[n-1][n-1], and                          //
//     x[i] = [B[i] - (U[i][i+1] * x[i+1]  + ... + U[i][n-1] * x[n-1])]       //
//                                                                 / U[i][i], //
//     for i = n-2, ..., 0.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *U   Pointer to the first element of the upper triangular       //
//                 matrix.                                                    //
//     double *B   Pointer to the column vector, (n x 1) matrix, B.           //
//     double *x   Pointer to the column vector, (n x 1) matrix, x.           //
//     int     n   The number of rows or columns of the matrix U.             //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix U is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N * (N + 1) / 2], B[N], x[N];                                 //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     err = Upper_Triangular_Solve_ut(A, B, x, n);                           //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else printf(" The solution is \n");                                    //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Upper_Triangular_Solve_ut (double *U, double B[], double x[], int n)
{
  int i, k;

//         Solve the linear equation Ux = B for x, where U is an upper
//         triangular matrix.

  for (k = n - 1, U += ((n * (n + 1)) >> 1) - 1; k >= 0; U -= n - --k)
    {
      if (*U == 0.0)
	return -1;		// The matrix U is singular
      x[k] = B[k];
      for (i = k + 1; i < n; i++)
	x[k] -= x[i] * *(U + i - k);
      x[k] /= *U;
    }

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  int Upper_Triangular_Inverse_ut(double *U,  int n)                        //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the inverse of the upper triangular matrix U.  //
//                                                                            //
//     The matrix U is stored in upper triangular form, i.e. the elements of  //
//     are stored sequentially starting with the first row which contains n   //
//     elements, followed by the second row which consists of n-1 elements,   //
//     and so on until the last row which consists of 1 element.              //
//                                                                            //
//     The algorithm follows:                                                 //
//        Let M be the inverse of U, then U M = I,                            //
//     M[n-1][n-1] = 1.0 / U[n-1][n-1] and                                    //
//     M[i][j] = -( U[i][i+1] M[i+1][j] + ... + U[i][j] M[j][j] ) / U[i][i],  //
//     for i = n-2, ... , 0,  j = n-1, ..., i+1.                              //
//                                                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *U   On input, the pointer to the first element of the matrix   //
//                 whose upper triangular elements form the matrix which is   //
//                 to be inverted. On output, the upper triangular part is    //
//                 replaced by the inverse.  The subdiagonal elements are     //
//                 not modified.                                              //
//     int     n   The number of rows and/or columns of the matrix U.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix U is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double U[N * (N + 1) / 2];                                             //
//                                                                            //
//     (your code to create the matrix U)                                     //
//     err = Upper_Triangular_Inverse_ut(U, N);                               //
//     if (err < 0) printf(" Matrix U is singular\n");                        //
//     else {                                                                 //
//        printf(" The inverse is \n");                                       //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Upper_Triangular_Inverse_ut (double *U, int n)
{
  int i, j, k;
  double *p_i, *p_k;
  double sum;

//         Invert the diagonal elements of the upper triangular matrix U.

  for (k = 0, p_k = U; k < n; p_k += (n - k++))
    {
      if (*p_k == 0.0)
	return -1;
      else
	*p_k = 1.0 / *p_k;
    }

//         Invert the remaining upper triangular matrix U.

  for (i = n - 1, p_i = U + ((n * (n + 1)) >> 1) - 1; i >= 0; p_i -= n - --i)
    {
      for (j = n - 1; j > i; j--)
	{
	  sum = 0.0;
	  for (k = i + 1, p_k = p_i + n - i; k <= j; p_k += n - k++)
	    {
	      sum += *(p_i + k - i) * *(p_k + j - k);
	    }
	  *(p_i + j - i) = -*p_i * sum;
	}
    }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// File: doolittle.c                                                          //
// Routines:                                                                  //
//    Doolittle_LU_Decomposition                                              //
//    Doolittle_LU_Solve                                                      //
//                                                                            //
// Required Externally Defined Routines:                                      //
//    Unit_Lower_Triangular_Solve                                             //
//    Upper_Triangular_Solve                                                  //
////////////////////////////////////////////////////////////////////////////////

//                    Required Externally Defined Routines
void Unit_Lower_Triangular_Solve (double *L, double B[], double x[], int n);
int Upper_Triangular_Solve (double *U, double B[], double x[], int n);

////////////////////////////////////////////////////////////////////////////////
//  int Doolittle_LU_Decomposition(double *A, int n)                          //
//                                                                            //
//  Description:                                                              //
//     This routine uses Doolittle's method to decompose the n x n matrix A   //
//     into a unit lower triangular matrix L and an upper triangular matrix U //
//     such that A = LU.                                                      //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.                                                        //
//     Note!  In Doolittle's method the diagonal elements of L are 1 and are  //
//            not stored.                                                     //
//     Note!  The determinant of A is the product of the diagonal elements    //
//            of U.  (det A = det L * det U = det U).                         //
//     This routine is suitable for those classes of matrices which when      //
//     performing Gaussian elimination do not need to undergo partial         //
//     pivoting, e.g. positive definite symmetric matrices, diagonally        //
//     dominant band matrices, etc.                                           //
//     For the more general case in which partial pivoting is needed use      //
//                  Doolittle_LU_Decomposition_with_Pivoting.                 //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is Doolittle_LU_Solve      //
//     (see below).                                                           //
//                                                                            //
//     The Doolittle method is given by evaluating, in order, the following   //
//     pair of expressions for k = 0, ... , n-1:                              //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                 for j = k, k+1, ... , n-1                  //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                          / U[k][k]                                         //
//                                 for i = k+1, ... , n-1.                    //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
//       forms the lower triangular matrix.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A   Pointer to the first element of the matrix A[n][n].        //
//     int     n   The number of rows or columns of the matrix A.             //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N];                                                        //
//                                                                            //
//     (your code to intialize the matrix A)                                  //
//                                                                            //
//     err = Doolittle_LU_Decomposition(&A[0][0], N);                         //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else { printf(" The LU decomposition of A is \n");                     //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Doolittle_LU_Decomposition (double *A, int n)
{
  int i, j, k, p;
  double *p_k, *p_row, *p_col;

//         For each row and column, k = 0, ..., n-1,
//            find the upper triangular matrix elements for row k
//            and if the matrix is non-singular (nonzero diagonal element).
//            find the lower triangular matrix elements for column k.

  for (k = 0, p_k = A; k < n; p_k += n, k++)
    {
      for (j = k; j < n; j++)
	{
	  for (p = 0, p_col = A; p < k; p_col += n, p++)
	    *(p_k + j) -= *(p_k + p) * *(p_col + j);
	}
      if (*(p_k + k) == 0.0)
	return -1;
      for (i = k + 1, p_row = p_k + n; i < n; p_row += n, i++)
	{
	  for (p = 0, p_col = A; p < k; p_col += n, p++)
	    *(p_row + k) -= *(p_row + p) * *(p_col + k);
	  *(p_row + k) /= *(p_k + k);
	}
    }
  return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  int Doolittle_LU_Solve(double *LU, double *B, double *x,  int n)           //
//                                                                            //
//  Description:                                                              //
//     This routine uses Doolittle's method to solve the linear equation      //
//     Ax = B.  This routine is called after the matrix A has been decomposed //
//     into a product of a unit lower triangular matrix L and an upper        //
//     triangular matrix U without pivoting.  The argument LU is a pointer to //
//     the matrix the subdiagonal part of which is L and the superdiagonal    //
//     together with the diagonal part is U. (The diagonal part of L is 1 and //
//     is not stored.)   The matrix A = LU.                                   //
//     The solution proceeds by solving the linear equation Ly = B for y and  //
//     subsequently solving the linear equation Ux = y for x.                 //
//                                                                            //
//  Arguments:                                                                //
//     double *LU  Pointer to the first element of the matrix whose elements  //
//                 form the lower and upper triangular matrix factors of A.   //
//     double *B   Pointer to the column vector, (n x 1) matrix, B            //
//     double *x   Solution to the equation Ax = B.                           //
//     int     n   The number of rows or columns of the matrix LU.            //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], x[N];                                            //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     err = Doolittle_LU_Decomposition(&A[0][0], N);                         //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else {                                                                 //
//        err = Doolittle_LU_Solve(&A[0][0], B, x, n);                        //
//        if (err < 0) printf(" Matrix A is singular\n");                     //
//        else printf(" The solution is \n");                                 //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Doolittle_LU_Solve (double *LU, double B[], double x[], int n)
{

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix with an implied 1 along the diagonal.

  Unit_Lower_Triangular_Solve (LU, B, x, n);

//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Lx = B and U is an upper triangular matrix.

  return Upper_Triangular_Solve (LU, x, x, n);
}

////////////////////////////////////////////////////////////////////////////////
// File: doolittle_pivot.c                                                    //
// Routines:                                                                  //
//    Doolittle_LU_Decomposition_with_Pivoting                                //
//    Doolittle_LU_with_Pivoting_Solve                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Doolittle_LU_Decomposition_with_Pivoting(double *A, int pivot[],      //
//                                                                    int n)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Doolittle's method to decompose the n x n matrix A   //
//     into a unit lower triangular matrix L and an upper triangular matrix U //
//     such that A = LU.                                                      //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.  Note!  In Doolittle's method the diagonal elements of //
//     L are 1 and are not stored.                                            //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is                         //
//                      Doolittle_LU_with_Pivoting_Solve.                     //
//                                                                            //
//     The Doolittle method with partial pivoting is:  Determine the pivot    //
//     row and interchange the current row with the pivot row, then assuming  //
//     that row k is the current row, k = 0, ..., n - 1 evaluate in order the //
//     following pair of expressions                                          //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                 for j = k, k+1, ... , n-1                  //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                          / U[k][k]                                         //
//                                 for i = k+1, ... , n-1.                    //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
//       forms the lower triangular matrix.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A       Pointer to the first element of the matrix A[n][n].    //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     int     n       The number of rows or columns of the matrix A.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N];                                                        //
//     int    pivot[N];                                                       //
//                                                                            //
//     (your code to create matrix A)                                         //
//     err = Doolittle_LU_Decomposition_with_Pivoting(&A[0][0], pivot, N);    //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else { printf(" The LU decomposition of A is \n");                     //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //



int
Doolittle_LU_Decomposition_with_Pivoting (double *A, int pivot[], int n)
{
  int i, j, k;
  double *p_k, *p_row, *p_col;
  double max;


//         For each row and column, k = 0, ..., n-1,

  for (k = 0, p_k = A; k < n; p_k += n, k++)
    {

//            find the pivot row

      pivot[k] = k;
      max = fabs (*(p_k + k));
      for (j = k + 1, p_row = p_k + n; j < n; j++, p_row += n)
	{
	  if (max < fabs (*(p_row + k)))
	    {
	      max = fabs (*(p_row + k));
	      pivot[k] = j;
	      p_col = p_row;
	    }
	}

//     and if the pivot row differs from the current row, then
//     interchange the two rows.

      if (pivot[k] != k)
	for (j = 0; j < n; j++)
	  {
	    max = *(p_k + j);
	    *(p_k + j) = *(p_col + j);
	    *(p_col + j) = max;
	  }

//                and if the matrix is singular, return error


      if (*(p_k + k) == 0.0)
	return -1;

//      otherwise find the lower triangular matrix elements for column k.

      for (i = k + 1, p_row = p_k + n; i < n; p_row += n, i++)
	{
	  *(p_row + k) /= *(p_k + k);
	}

//            update remaining matrix

      for (i = k + 1, p_row = p_k + n; i < n; p_row += n, i++)
	for (j = k + 1; j < n; j++)
	  *(p_row + j) -= *(p_row + k) * *(p_k + j);

    }

  return 0;
}



////////////////////////////////////////////////////////////////////////////////
//  int Doolittle_LU_with_Pivoting_Solve(double *LU, double *B, int pivot[],  //
//                                       double *x,  int n)                   //
//                                                                            //
//  Description:                                                              //
//     This routine uses Doolittle's method to solve the linear equation      //
//     Ax = B.  This routine is called after the matrix A has been decomposed //
//     into a product of a unit lower triangular matrix L and an upper        //
//     triangular matrix U with pivoting.  The argument LU is a pointer to the//
//     matrix the subdiagonal part of which is L and the superdiagonal        //
//     together with the diagonal part is U. (The diagonal part of L is 1 and //
//     is not stored.)   The matrix A = LU.                                   //
//     The solution proceeds by solving the linear equation Ly = B for y and  //
//     subsequently solving the linear equation Ux = y for x.                 //
//                                                                            //
//  Arguments:                                                                //
//     double *LU      Pointer to the first element of the matrix whose       //
//                     elements form the lower and upper triangular matrix    //
//                     factors of A.                                          //
//     double *B       Pointer to the column vector, (n x 1) matrix, B.       //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     double *x       Solution to the equation Ax = B.                       //
//     int     n       The number of rows or columns of the matrix LU.        //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], x[N];                                            //
//     int    pivot[N];                                                       //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     err = Doolittle_LU_Decomposition_with_Pivoting(&A[0][0], pivot,  N);   //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else {                                                                 //
//        err = Doolittle_LU_with_Pivoting_Solve(&A[0][0], B, pivot, x, N);   //
//        if (err < 0) printf(" Matrix A is singular\n");                     //
//        else printf(" The solution is \n");                                 //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Doolittle_LU_with_Pivoting_Solve (double *A, double B[], int pivot[],
				  double x[], int n)
{
  int i, k;
  double *p_k;
  double dum;

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix with an implied 1 along the diagonal.

  for (k = 0, p_k = A; k < n; p_k += n, k++)
    {
      if (pivot[k] != k)
	{
	  dum = B[k];
	  B[k] = B[pivot[k]];
	  B[pivot[k]] = dum;
	}
      x[k] = B[k];
      for (i = 0; i < k; i++)
	x[k] -= x[i] * *(p_k + i);
    }

//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Lx = B and U is an upper triangular matrix.

  for (k = n - 1, p_k = A + n * (n - 1); k >= 0; k--, p_k -= n)
    {
      if (pivot[k] != k)
	{
	  dum = B[k];
	  B[k] = B[pivot[k]];
	  B[pivot[k]] = dum;
	}
      for (i = k + 1; i < n; i++)
	x[k] -= x[i] * *(p_k + i);
      if (*(p_k + k) == 0.0)
	return -1;
      x[k] /= *(p_k + k);
    }

  return 0;
}

int Gaussian_Elimination (double *A, double B[], int pivot[], double X[],
			  int n);
static void print_results (void);

double A[N][N], AC[N][N], B[N], BC[N];
double X[N];
int pivot[N];
////////////////////////////////////////////////////////////////////////////////
// File: gauss_elimination.c                                                  //
// Routines:                                                                  //
//    Gaussian_Elimination                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Gaussian_Elimination(double *A, int n, double *B)                     //
//                                                                            //
//     Solve the linear system of equations AX=B where A is an n x n matrix   //
//     B is an n-dimensional column vector (n x 1 matrix) for the             //
//     n-dimensional column vector (n x 1 matrix) X.                          //
//                                                                            //
//     This routine performs partial pivoting and the elements of A are       //
//     modified during computation.  The result X is returned in B.           //
//     If the matrix A is singular, the return value of the function call is  //
//     -1. If the solution was found, the function return value is 0.         //
//                                                                            //
//  Arguments:                                                                //
//     double *A      On input, the pointer to the first element of the       //
//                    matrix A[n][n].  On output, the matrix A is destroyed.  //
//     int     n      The number of rows and columns of the matrix A and the  //
//                    dimension of B.                                         //
//     double *B      On input, the pointer to the first element of the       //
//                    vector B[n].  On output, the vector B is replaced by the//
//                    vector X, the solution of AX = B.                       //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N];                                                  //
//                                                                            //
//     (your code to create the matrix A and vector B )                       //
//     err = Gaussian_Elimination((double*)A, NROWS, B);                      //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else { printf(" The Solution is: \n"); ...                             //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int
Gaussian_Elimination (double *A, double B[], int pivot[], double X[], int n)
{
  int row, i, j, pivot_row;
  double max, dum, *pa, *pA, *A_pivot_row;
  // for each variable find pivot row and perform forward substitution
  pa = A;
  for (row = 0; row < (n - 1); row++, pa += n)
    {
      //  find the pivot row
      A_pivot_row = pa;
      max = fabs (*(pa + row));
      pA = pa + n;
      pivot_row = row;
      for (i = row + 1; i < n; pA += n, i++)
	if ((dum = fabs (*(pA + row))) > max)
	  {
	    max = dum;
	    A_pivot_row = pA;
	    pivot_row = i;
	  }
      if (max == 0.0)
	return -1;		// the matrix A is singular
      // and if it differs from the current row, interchange the two rows.
      if (pivot_row != row)
	{
	  for (i = row; i < n; i++)
	    {
	      dum = *(pa + i);
	      *(pa + i) = *(A_pivot_row + i);
	      *(A_pivot_row + i) = dum;
	    }
	  dum = B[row];
	  B[row] = B[pivot_row];
	  B[pivot_row] = dum;
	}
      // Perform forward substitution
      for (i = row + 1; i < n; i++)
	{
	  pA = A + i * n;
	  dum = -*(pA + row) / *(pa + row);
	  *(pA + row) = 0.0;
	  for (j = row + 1; j < n; j++)
	    *(pA + j) += dum * *(pa + j);
	  B[i] += dum * B[row];
	}
    }
  // Perform backward substitution
  pa = A + (n - 1) * n;
  for (row = n - 1; row >= 0; pa -= n, row--)
    {
      if (*(pa + row) == 0.0)
	return -1;		// matrix is singular
      dum = 1.0 / *(pa + row);
      for (i = row + 1; i < n; i++)
	*(pa + i) *= dum;
      B[row] *= dum;
      for (i = 0, pA = A; i < row; pA += n, i++)
	{
	  dum = *(pA + row);
	  for (j = row + 1; j < n; j++)
	    *(pA + j) -= dum * *(pa + j);
	  B[i] -= dum * B[row];
	}
      pivot[i]++;
    }
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
// File: multiply_matrices.c                                                  //
// Routine(s):                                                                //
//    Multiply_Matrices                                                       //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Multiply_Matrices(double *C, double *A, int nrows, int ncols,        //
//                                                    double *B, int mcols)   //
//                                                                            //
//  Description:                                                              //
//     Post multiply the nrows x ncols matrix A by the ncols x mcols matrix   //
//     B to form the nrows x mcols matrix C, i.e. C = A B.                    //
//     The matrix C should be declared as double C[nrows][mcols] in the       //
//     calling routine.  The memory allocated to C should not include any     //
//     memory allocated to A or B.                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *C    Pointer to the first element of the matrix C.             //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of the matrices A and C.               //
//     int    ncols The number of columns of the matrices A and the           //
//                   number of rows of the matrix B.                          //
//     double *B    Pointer to the first element of the matrix B.             //
//     int    mcols The number of columns of the matrices B and C.            //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     #define NB                                                             //
//     double A[M][N],  B[N][NB], C[M][NB];                                   //
//                                                                            //
//     (your code to initialize the matrices A and B)                         //
//                                                                            //
//     Multiply_Matrices(&C[0][0], &A[0][0], M, N, &B[0][0], NB);             //
//     printf("The matrix C is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Multiply_Matrices(double *C, double *A, int nrows, int ncols, double *B, int mcols)
{
   double *pB;
   double *p_B;
   int i,j,k;
   for (i = 0; i < nrows; A += ncols, i++)
      for (p_B = B, j = 0; j < mcols; C++, p_B++, j++) {
         pB = p_B;
         *C = 0.0;
         for (k = 0; k < ncols; pB += mcols, k++)
            *C += *(A+k) * *pB;
      }
}
void
Input_data ()
{
  int i, j, n;
  n = N;
  int row, col;
  row = n;
  col = n;
  double A[n][n] = { {5.0, 7.0, 3.0}, {7.0, 11.0, 2.0}, {3.0, 2.0, 6.0} };
  double B[n] = { 2.0, 3.0, 5.0 };
  double array_A[n*n];
  printf ("Prog: LU_Solve.cpp - Solve the system AX = B using the LU decomposition method.\n\n\n");
  // Save matrices A and vector B, the routine Gaussian_Elimination() //
  // Destroys A and the Solution is returned in B.                    //
  for (i = 0; i < row; i++)
  {
	for (j = 0; j < col; j++)
	{
		array_A[i * row + j] = A[i][j];
	}
  }
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
	{
	  AC[i][j] = A[i][j];
	}
      BC[i] = B[i];
    }
  X[i] = 0;
  X[i] = B[i];
  pivot[i] = 0;
  Doolittle_LU_Decomposition (&A[0][0], N);
  Doolittle_LU_Solve (&A[0][0], &B[0], &X[0], N);
  print_results ();
  printf ("with Pivoting\n");
  Doolittle_LU_Decomposition_with_Pivoting (&A[0][0], &pivot[0], N);
  print_results ();
}

void
print_results (void)
{
  int i, j;
  printf ("******************** Solve Ax = B ********************\n\n");
  printf ("where A = \n");
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
	printf ("%6.3f   ", AC[i][j]);
      printf ("\n");
    }
  printf ("and B = \n");
  for (i = 0; i < N; i++)
    printf ("%6.3f   ", BC[i]);
  printf ("\n\n");
  printf ("The solution is X = \n");
  for (i = 0; i < N; i++)
    {
      printf ("%6.3f   ", X[i]);
      // Assigning values to the arrays X and pivot to fix warning message
    }
  printf ("\n\n");
  printf ("\n\n\n");
  return;
}

int
main ()
{
  Input_data ();
  return 0;
}
