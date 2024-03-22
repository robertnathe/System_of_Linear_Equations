#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;
template<typename T>
using matrix = std::vector<std::vector<T>>;

int Tridiagonal(const matrix<double>& A, double subdiagonal[], double diagonal[], double superdiagonal[], double array_A[], double B[], vector<double>& X, int num_rows, int num_cols);
int Tridiagonal_LU_Decomposition( double subdiagonal[], double diagonal[], double superdiagonal[], int num_cols);
int Tridiagonal_LU_Solve( double subdiagonal[], double diagonal[], double superdiagonal[], double B[], double x[], int num_cols);
int Jacobi(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter);
int Jacobi_Residual(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter);
int Gauss_Seidel(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter);
int SOR(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter);
int ReadTridiagonal(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols);
int ReadMatrixMM(matrix<double> & A, int& num_rows, int& num_cols); // For matrices 
int ReadVectorMM(vector<double>& B, int& num_cols); // For vectors
int WriteMatrixMM(const matrix<double>& A, int num_rows, int num_cols);
int WriteVectorMM(const std::vector<double>& X, int num_rows, int num_cols);
void PrintVector(const std::vector<double>& vector);
void PrintVector2D(const matrix<double>& A);
int print_results (double subdiagonal[], double diagonal[], double superdiagonal[], double array_A[], double B[], vector<double>& X, int num_rows, int num_cols);
std::vector<double> CreateZeroVector(int num_rows);

int Tridiagonal(const matrix<double>& A, double subdiagonal[], double diagonal[], double superdiagonal[], double array_A[], double B[], vector<double>& X, int num_rows, int num_cols)
{
  int i;
  double m;
  //forward elimination
  for (i=1;i<=num_cols-1;i++)
  {
    m = subdiagonal[i]/diagonal[i-1];
    diagonal[i] = diagonal[i] - (m * superdiagonal[i-1]);
    B[i] = B[i] - (m * B[i-1]);
  }
  //backward substitution
  X[num_cols]=B[num_cols]/diagonal[num_cols];
  for(i=num_cols-1;i>=0;i--)
  {
    X[i]=(B[i]-superdiagonal[i]*X[i+1])/diagonal[i];
  }
  print_results (&subdiagonal[0], &diagonal[0], &superdiagonal[0], &array_A[0], &B[0], X, num_rows, num_cols);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// File: tridiagonal.c                                                        //
// Contents:                                                                  //
//    Tridiagonal_LU_Decomposition                                            //
//    Tridiagonal_LU_Solve                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Tridiagonal_LU_Decomposition(double *subdiagonal, double *diagonal,       //
//                                           double *superdiagonal, int n )   //
//                                                                            //
//  Description:                                                              //
//     This routine decomposes a tridiagonal matrix A into the product of     //
//     a unit lower triangular (bidiagonal) matrix L and an upper triangular  //
//     (bidiagonal) matrix U, A = LU.                                         //
//     The tridiagonal matrix A is defined by the three vectors, subdiagonal, //
//     diagonal, and superdiagonal, where the i-th component of subdiagonal is//
//     subdiagonal[i] = A[i+1][i], for i = 0, ..., n - 2; the i-th component  //
//     of diagonal is diagonal[i] = A[i][i], for i = 0, ..., n - 1; and the   //
//     i-th component of superdiagonal is superdiagonal[i] = A[i][i+1], for   //
//     i = 0, ..., n - 2.                                                     //
//     The algorithm proceeds by decomposing the matrix A into the product    //
//     of a unit lower triangular (bidiagonal) matrix, stored in subdiagonal, //
//     and an upper triangular (bidiagonal) matrix, stored in diagonal and    //
//     and superdiagonal.                                                     //
//     After performing the LU decomposition for A, call Tridiagonal_LU_Solve //
//     to solve the equation Ax = B for x given B.                            //
//                                                                            //
//     This routine can fail if A[0][0] = 0 or if during the LU decomposition //
//     the diagonal element of U becomes 0.  This does not imply that the     //
//     matrix A is singular.  If A is positive definite or if A is diagonally //
//     dominant then the procedure should not fail.                           //
//                                                                            //
//  Arguments:                                                                //
//     double subdiagonal[]                                                   //
//        On input, subdiagonal[i] is the subdiagonal element A[i+1][i].      //
//        On output, subdiagonal[i] is the subdiagonal of the unit lower      //
//        triangular matrix L in the LU decomposition of A.                   //
//     double diagonal[]                                                      //
//        On input, diagonal[i] is the diagonal element A[i][i] of the matrix //
//        A.  On output, diagonal[i] is the diagonal of the upper triangular  //
//        matrix U in the LU decomposition of A.                              //
//     double superdiagonal[]                                                 //
//        On input, superdiagonal[i] is the superdiagonal element A[i][i+1] of//
//        the matrix A.  On output, superdiagonal[i] is the superdiagonal of  //
//        the upper triangular matrix U, which agrees with the input.         //
//     int     n   The number of rows and/or columns of the matrix A.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - A zero occurred on the diagonal of U.                     //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double subdiagonal[N], diagonal[N], superdiagonal[N];                  //
//                                                                            //
//     (your code to create subdiagonal, diagonal, and superdiagonal)         //
//     err = Tridiagonal_LU_Decomposition(subdiagonal, diagonal,              //
//                                                         superdiagonal, N); //
//     if (err < 0) printf(" Matrix A failed the LU decomposition\n");        //
//     else { printf(" The Solution is: \n"); ...                             //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Tridiagonal_LU_Decomposition( double subdiagonal[], double diagonal[],
                                                double superdiagonal[], int num_cols )
{
   int i;
   for (i = 0; i < (num_cols-1); i++) {
      if (diagonal[i] == 0.0) return -1;
      subdiagonal[i] /= diagonal[i];
      diagonal[i+1] -= subdiagonal[i] * superdiagonal[i];
   }
   if (diagonal[num_cols-1] == 0.0) return -1;
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
// int Tridiagonal_LU_Solve(double subdiagonal[], double diagonal[],          //
//                 double superdiagonal[],  double B[], double x[],  int n)   //
//                                                                            //
//  Description:                                                              //
//     This routine uses the LU decomposition from the routine above,         //
//     Tridiagonal_LU_Decomposition, to solve the linear equation Ax = B,     //
//     where A = LU, L is the unit lower triangular (bidiagonal) matrix with  //
//     subdiagonal subdiagonal[] and diagonal all 1's, and U is the upper     //
//     triangular (bidiagonal) matrix with diagonal diagonal[] and            //
//     superdiagonal superdiagonal[].                                         //
//     The solution proceeds by solving the linear equation Ly = B for y and  //
//     subsequently solving the linear equation Ux = y for x.                 //
//                                                                            //
//  Arguments:                                                                //
//     double subdiagonal[]                                                   //
//        The subdiagonal of the unit lower triangular matrix L in the LU     //
//        decomposition of A.                                                 //
//     double diagonal[]                                                      //
//        The diagonal of the upper triangular matrix U in the LU decomposi-  //
//        tion of A.                                                          //
//     double superdiagonal[]                                                 //
//        The superdiagonal of the upper triangular matrix U.                 //
//     double B[]                                                             //
//        Pointer to the column vector, (n x 1) matrix, B.                    //
//     double x[]                                                             //
//        Solution to the equation Ax = B.                                    //
//     int     n   The number of rows and/or columns of the matrix LU.        //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix U is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double subdiagonal[N], diagonal[N], superdiagonal[N];                  //
//     double B[N], x[N];                                                     //
//                                                                            //
//     (your code to create subdiagonal, diagonal, superdiagonal, and B)      //
//     err = Tridiagonal_LU_Decomposition(subdiagonal, diagonal,              //
//                                                          superdiagonal, N);//
//     if (err < 0) printf(" Matrix A is fails the LU decomposition\n");      //
//     else {                                                                 //
//        err = Tridiagonal_LU_Solve(subdiagona, diagonal, superdiagonal, B,  //
//                                                                      x, n);//
//        if (err < 0) printf(" Matrix A is singular\n");                     //
//        else printf(" The solution is \n");                                 //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Tridiagonal_LU_Solve( double subdiagonal[], double diagonal[],
                       double superdiagonal[], double B[], double x[], int num_cols)
{
   int i;
//         Check that all diagonal elements are nonzero.
//         If a diagonal element is zero then U is singular, so return
//         signalling an error.
   for (i = 0; i < num_cols; i++) if (diagonal[i] == 0.0) return -1;
//         Solve the linear equation Ly = B for y, where L is a lower
//         triangular matrix.
   x[0] = B[0];
   for (i = 1; i < num_cols; i++) x[i] = B[i] - subdiagonal[i-1] * x[i-1];
//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Ly = B and U is an upper triangular matrix.
   x[num_cols-1] /= diagonal[num_cols-1];
   for (i = num_cols-2; i >= 0; i--) {
      x[i] -= superdiagonal[i] * x[i+1];
      x[i] /= diagonal[i];
   }
   return 0;
}

int Jacobi(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter) {

 std::cout << std::endl << "Implementing the Jacobi method: " << endl;
 // Initialize previous iteration vector
 vector<double> XPrevious(num_rows, 0.0);
 // Loop for iterations
 for (int Iter = 1; Iter <= Max_Iter; ++Iter) {
   double distance_max = 0.0;
   // Update each element of X
   for (int i = 0; i < num_rows; ++i) {
     double sum = 0.0;
     for (int j = 0; j < num_cols; ++j) {
       if (i != j) {
         sum += A[i][j] * XPrevious[j];
       }
     }
     X[i] = (B[i] - sum) / A[i][i];
   }
   // Calculate distance metric and check for convergence
   for (int j = 0; j < num_cols; ++j) {
     double distance = std::abs(X[j] - XPrevious[j]);
     distance_max = std::max(distance_max, distance);
   }
   if (distance_max < tolerance) {
     std::cout << "The convergence distance is less than the tolerance of " << tolerance << "." << endl;
     std::cout << "The number of iterations is " << Iter << "." << endl;
     return 0; // Convergence achieved
   }
   // Update previous iteration vector
   XPrevious = X;
 }
 // Maximum iterations reached without convergence
 std::cout << "Maximum number of iterations reached before Jacobi could finish." << std::endl;
 return 1; // Not converged
}

int Jacobi_Residual(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter) {
  // Initialize iteration counter and residual
  int iter = 0;
  double residual = 1.0;
  // Loop until convergence or maximum iterations reached
  while (iter < Max_Iter && residual > tolerance) {
    // Update each element of X using Jacobi iteration
    residual = 0.0;
//    #pragma omp parallel for reduction(+:residual)
    for (int i = 0; i < num_rows; ++i) {
      double new_x = 0.0;
      for (int j = 0; j < num_cols; ++j) {
        if (i != j) {
          new_x += A[i][j] * X[j];
        }
      }
      new_x = (B[i] - new_x) / A[i][i];
      residual += std::pow(X[i] - new_x, 2.0);
      X[i] = new_x;
    }
    iter++;
  }
  // Return number of iterations or -1 if not converged
  return (residual <= tolerance) ? iter : -1;
}

int Gauss_Seidel (const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter)
{
 std::cout << endl << "Implementing the Gauss Seidel method: " << endl; 
 bool convergence_check {false};
 int i {0}, j {0}, Iter {1};
 double elem {0.0}, distance {0.0}, distance_max {0.0};
 double temporary_variable_one {0.0}, temporary_variable_two {0.0};
 double distance_squared {0.0};
 double distance_new {0.0};
 double norm_variable_squared {0.0};
 double norm_variable {0.0};
 double tolerance_variable {0.0};
 vector<double> XPrevious;
 for (int j=0; j<num_cols; j++) {
   XPrevious.push_back(elem);
 }
 std::cout << std::setprecision(5);
 while(convergence_check != true) {   
   distance_max = 0.0;   
   for (i = 0; i < num_rows; i++)
   {
     X[i] = (B[i] / A[i][i]);
     for (j = 0; j < num_rows; j++)
     {
       if (j == i) // Without the diagonals
         continue;
       X[i] = X[i] - ((A[i][j] / A[i][i]) * XPrevious[j]);
     }
     distance = abs( X[i] - XPrevious[i] );
    distance_max = std::max(distance, distance_max);
    temporary_variable_one = abs( X[i] - XPrevious[i] );
     distance_squared += temporary_variable_one * temporary_variable_one;
     temporary_variable_two = abs( X[i] - XPrevious[i] );
     norm_variable_squared += temporary_variable_two;
     XPrevious[i] = X[i];
   }
   distance_new = sqrt(distance_squared);
   norm_variable = sqrt(norm_variable_squared);
   tolerance_variable = distance_new / norm_variable;
   distance_squared = 0.0;
   norm_variable_squared = 0.0; 
   if (norm_variable > 100.0)
   {
     break;
   }
   if ( (tolerance_variable) < tolerance) { 
//   if ( (distance_max) < tolerance) {
     std::cout << "The convergence distance is less than the " << tolerance << "." << endl;
    std::cout << "The number of iterations is " << Iter << "." << endl;
    convergence_check = true;
     break;
   }
   if (Iter >= Max_Iter) {
    std::cout << "The distance_max is the following: " << distance_max << endl;
    std::cout << "The tolerance is the following: " << tolerance << endl;
    std::cout << "The number of iterations is " << Iter << "." << endl;
     std::cout << "Maximum number of iterations reached before Jacobi_Iterative_Solve could finish." << std::endl;
     convergence_check = true;
     break;
   }   
   Iter++;
 }
 convergence_check = false;
 XPrevious.clear();
 XPrevious.shrink_to_fit();
 return 0;
}

int SOR (const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter)
{ 
 std::cout << endl << "Implementing the SOR method: " << endl;
 bool convergence_check {false};
 int i {0}, j {0}, Iter {1};
 double elem {1.0}, omega {1.25}, distance {0.0}, distance_max {0.0};
 double temporary_variable_one {0.0}, temporary_variable_two {0.0};
 double distance_squared {0.0};
 double distance_new {0.0};
 double norm_variable_squared {0.0};
 double norm_variable {0.0};
 double tolerance_variable {0.0}; 
 vector<double> XPrevious;
// Change omega to compare the Gauss-Seidel and SOR methods.
// std::cout << "Please Enter the omega factor for the SOR method." << endl;
// std::cin >> omega;
 for (j=0; j<num_cols; j++) {
   XPrevious.push_back(elem);
 }
 std::cout << std::setprecision(5); 
 while(convergence_check != true) {   
   distance_max = 0.0;   
   for (i = 0; i < num_rows; i++)
   {
     X[i] = (B[i] / A[i][i]);
     for (j = 0; j < num_rows; j++)
     {
       if (j == i) // Without the diagonals
         continue;
       X[i] = X[i] - ((A[i][j] / A[i][i]) * XPrevious[j]);
     }
     // SOR factor omega calculated on the next line
     X[i] = omega*X[i] +((1-omega)*XPrevious[i]);
     distance = abs( X[i] - XPrevious[i] );
    distance_max = std::max(distance, distance_max);
    temporary_variable_one = abs( X[i] - XPrevious[i] );
     distance_squared += temporary_variable_one * temporary_variable_one;
     temporary_variable_two = abs( X[i] - XPrevious[i] );
     norm_variable_squared += temporary_variable_two;
     XPrevious[i] = X[i];
   }
   distance_new = sqrt(distance_squared);
   norm_variable = sqrt(norm_variable_squared);
   tolerance_variable = distance_new / norm_variable;
   distance_squared = 0.0;
   norm_variable_squared = 0.0; 
   if (norm_variable > 100.0)
   {
     break;
   }
   if ( (tolerance_variable) < tolerance) { 
//   if ( (distance_max) < tolerance) {
     std::cout << "The convergence distance is less than the " << tolerance << "." << endl;
    std::cout << "The number of iterations is " << Iter << "." << endl;
    convergence_check = true;
     break;
   }
   if (Iter >= Max_Iter) {
    std::cout << "The distance_max is the following: " << distance_max << endl;
    std::cout << "The tolerance is the following: " << tolerance << endl;
    std::cout << "The number of iterations is " << Iter << "." << endl;
     std::cout << "Maximum number of iterations reached before Jacobi_Iterative_Solve could finish." << std::endl;
     convergence_check = true;
     break;
   }   
   Iter++;
 }
 convergence_check = false;
 XPrevious.clear();
 XPrevious.shrink_to_fit();
 return 0;
}

int ReadTridiagonal(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols)
{ 
  double *superdiagonal;
  superdiagonal = new double[num_cols];
//  superdiagonal[num_cols] = 0.0;
  double *diagonal;
  diagonal = new double[num_cols];
  double *subdiagonal;
  subdiagonal = new double[num_cols];
//  subdiagonal[0] = 0.0;
  double *array_A;
  array_A = new double[num_rows*num_cols];
  for (int i = 0; i < num_rows; i++)
  {
    for (int j = 0; j < num_cols; j++)
	{
	  array_A[i * num_cols + j] = A[i][j];
	  if (i == j)
	  {
	    diagonal[i] = array_A[i * num_cols + j];
      }
      else if  (i == j-1)
      {
	    superdiagonal[i] = array_A[i*num_cols+j];
	  }
	  else if (i == j+1)
	  {
	    subdiagonal[i] = array_A[i*num_cols+j];
	  }
	}
  }
  Tridiagonal (A, &subdiagonal[0], &diagonal[0], &superdiagonal[0], &array_A[0], &B[0], X, num_rows, num_cols);
  return 0;
}

int ReadMatrixMM(matrix<double>& A, int& num_rows, int& num_cols) {

 int number_of_entries_A {0};
 int i_index {0}, j_index {0};
 double elem {0.0};
 FILE* myFile = fopen("A.dat", "r");
 if (myFile == NULL) {
   std::cerr << "Error Reading File" << endl;
   exit(0);
 }
 // Skip header comments
 fscanf(myFile, "%*[^\n]\n"); // Read and discard header line
 // Read dimensions
 if (fscanf(myFile, "%u %u %u\n", &num_rows, &num_cols, &number_of_entries_A) != 3) {
   std::cerr << "Error reading matrix dimensions from A.dat" << endl;
   fclose(myFile);
   return -1;
 }
 // Resize A to accommodate num_num_rows and num_num_num_cols
 A.resize(num_rows);
 for (int i = 0; i < num_rows; ++i) {
   A[i].resize(num_cols);
 }
 // Read non-zero elements by row and column indices
 for (int i = 0; i < number_of_entries_A; ++i) {
   if (fscanf(myFile, "%u %u %lf\n", &i_index, &j_index, &elem) != 3) {
     std::cerr << "Error reading matrix entries from A.dat" << endl;
     fclose(myFile);
     return -1;
   }
   i_index--; // Adjust for zero-based indexing
   j_index--;
   A[i_index][j_index] = elem;
 }
 fclose(myFile);
 return 0;
}

// ReadVectorMM implementation for vectors
int ReadVectorMM(vector<double>& B, int& num_cols) {
    (void) num_cols; 
   FILE *myFile2;
   myFile2 = fopen ("B.dat", "r");
   int dim_B_Array[3];
   int i_index {0}, j_index {0};
   double value {0.0};
   while (myFile2 == NULL)
   {
    std::cout << "Error Reading File" << endl;
     exit (0);
   } 
   fscanf (myFile2, "%*s %*s %*s %*s %*s");
   for (int i = 0; i < 3; i++)
   {
     fscanf (myFile2, "%u,", &dim_B_Array[i]);
   }
   for (int i = 0; i < dim_B_Array[1]; i++)
     B.push_back(0.0);
   for (int i = 0; i < dim_B_Array[1]; i++)
   {
     fscanf (myFile2, "%u,", &i_index);
     i_index--;
     fscanf (myFile2, "%u,", &j_index);
     j_index--;
     fscanf (myFile2, "%lf,", &value);
     if (value != 0.0) 
     {
       B[i] = value;
     }
   }
   fclose (myFile2);
 return 0;
}

int WriteMatrixMM(const matrix<double>& A, int num_rows, int num_cols) {
 // Open file in writing mode
 int CountNonZeroEntries {0};
 FILE* fptr = fopen("A_out.txt", "w");
 if (fptr == NULL) {
   std::cerr << "Error opening file for writing matrix A." << std::endl;
   return -1; // Indicate error
 }
 for (int i = 0; i < num_rows; ++i) {
   for (int j = 0; j < num_cols; ++j) {
    if (A[i][j] != 0.0)
    {
      CountNonZeroEntries++;
    }
   }
 }
 // Write Matrix Market header information for a sparse matrix
 fprintf(fptr, "%%MatrixMarket_Input_matrix_A.dat matrix coordinate pattern general\n");
// fprintf(fptr, "%d %d %d\n", num_rows, num_cols, num_rows*num_cols); // Count all entries
 fprintf(fptr, "%d %d %d\n", num_rows, num_cols, CountNonZeroEntries); // Count non-zero entries
 // Write only non-zero elements of A
 // I changed this section of code.
 for (int i = 0; i < num_rows; ++i) {
   for (int j = 0; j < num_cols; ++j) {
     if (A[i][j] != 0.0) 
     { // Check for non-zero value
       fprintf(fptr, "%u %u %lf\n", i + 1, j + 1, A[i][j]);
     }
   }
 }
 // Close the file
 fclose(fptr);
 return 0; // Indicate success
}

int WriteVectorMM(const std::vector<double>& X, int num_rows, int num_cols) {
 // Open file in writing mode
 FILE* fptr = fopen("X.dat", "w");
 if (fptr == NULL) {
   std::cerr << "Error opening file for writing X vector." << std::endl;
   return -1; // Indicate error
 }
 // Write Matrix Market header information
 fprintf(fptr, "%%MatrixMarket_Input_matrix_X.dat matrix coordinate pattern general\n");
 fprintf(fptr, "%d %d %d\n", 1, num_cols, num_cols); // All entries are assumed non-zero
 std::cout << "%%MatrixMarket_Input_matrix_X.dat matrix coordinate pattern general\n";
 // Write each element of X
 for (int i = 0; i < num_rows; ++i) {
   fprintf(fptr, "%u %u %lf\n", 1, i+1, X[i]); // Row index always 1 for a vector
   std::cout << 1 << "   " << i+1 << "   "<< X[i] << std::endl;
 }
 // Close the file
 fclose(fptr);
 return 0; // Indicate success
}

void PrintVector(const std::vector<double>& vector) {
 std::cout << "Displaying the vector: " << endl;
 std::cout << std::setprecision(5);
 for (const double& value : vector) {
   std::cout << value << " ";
 }
 std::cout << std::endl;
}

void PrintVector2D(const matrix<double>& A) {
 std::cout << "Displaying the 2D vector:" << endl;
 for (const auto& row : A) {
   for (const double& value : row) {
     std::cout << value << " ";
   }
   std::cout << std::endl;
 }
}

int print_results (double subdiagonal[], double diagonal[], double superdiagonal[], double array_A[], double B[], vector<double>& X, int num_rows, int num_cols)
{
  int i, j;
  std::cout << "The Tridiagonal output is the following: " << std::endl;
  printf ("******************** Solve Ax = B ********************\n\n");
  printf ("where A = \n");
  for (i = 0; i < num_rows; i++)
    {
      for (j = 0; j < num_cols; j++)
	{
	  printf ("%6.5f   ", array_A[i * num_cols + j]);
	}
      printf ("\n");
    }
  printf ("\n");
  printf ("and B = \n");
  for (i = 0; i < num_cols; i++)
    printf ("%6.5f   ", B[i]);
  printf ("\n\n");
  FILE *myFile3;
  myFile3 = fopen("X.dat","w+");
  if (myFile3 == NULL)
  {
    printf("Error writing to file.\n");
    exit(0);
  }
  fprintf(myFile3,"%%MatrixMarket_Output_vector_C.dat matrix coordinate pattern general\n");
  fprintf (myFile3,"%d %d %d\n", 1, num_rows, num_cols);
  printf ("The solution is X = \n");
  for (i = 0; i < num_cols; i++)
  {
    fprintf(myFile3, "%d %d %lf\n", 1, i+1, X[i]);
    printf ("%6.5f    ", X[i]);
  }
  fclose(myFile3);
  printf ("\n\n");
  return 0;
}

std::vector<double> CreateZeroVector(int num_rows) {
 return std::vector<double>(num_rows);
}

int main ()
{
   int num_rows {0}, num_cols {0};
   vector <double> B;
   ReadVectorMM(B, num_cols);
   matrix<double> A;
   ReadMatrixMM(A,num_rows,num_cols);         
   std::cout << std::setprecision(5);
   int Max_Iter {1000};
   double tolerance {0.00001};
   double omega {1.25};
   (void) omega;
   vector<double> X = CreateZeroVector(num_rows);  
   ReadTridiagonal(A, B, X, num_rows, num_cols);
   for (int j = 0; j < num_cols; j++)
   {
     X[j] = 0.0; 
   }
   ReadVectorMM(B, num_cols);  
   std::cout << "The Jacobi output is the following: " << std::endl;
   Jacobi (A, B, X, num_rows, num_cols, tolerance, Max_Iter);
   PrintVector(X);
   Gauss_Seidel (A, B, X, num_rows, num_cols, tolerance, Max_Iter);
   PrintVector(X);
//   WriteVectorMM(X,num_rows,num_cols);
//   WriteMatrixMM(A,num_rows,num_cols);
   SOR (A, B, X, num_rows, num_cols, tolerance, Max_Iter);
   PrintVector(X);
   A.clear();
   A.shrink_to_fit();
   B.clear();
   B.shrink_to_fit();
   X.clear();
   X.shrink_to_fit();
   return 0;
}
