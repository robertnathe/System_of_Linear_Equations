#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstddef> // for size_t

using namespace std;
template<typename T>
using matrix = std::vector<std::vector<T>>;

bool convergedUsingL1Norm(const std::vector<double>& X, const std::vector<double>& XPrevious, double tolerance, int size);
int Tridiagonal(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols);
int Tridiagonal_LU_Decomposition( double lowerDiagonal[], double mainDiagonal[], double upperDiagonal[], int num_cols);
int Tridiagonal_LU_Solve( double lowerDiagonal[], double mainDiagonal[], double upperDiagonal[], double B[], double x[], int num_cols);
int Jacobi_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter);
int Jacobi_L2_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter);
int Gauss_Seidel_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter);
int SOR_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter);
int ReadMatrixMarketMatrix(matrix<double> & A, int& num_rows, int& num_cols); // For matrices 
int ReadMatrixMarketVector(vector<double>& B, int& num_cols); // For vectors
int Write1DArrayToMatrixMarketFile(const double B[], int num_rows); 
int Write2DArrayToMatrixMarketFile(const double array_A[], int num_rows, int num_cols);
int WriteMatrixMarketMatrix(const matrix<double>& A, int num_rows, int num_cols);
int WriteMatrixMarketVector(const std::vector<double>& X, int num_rows, int num_cols);
int WriteTridiagonalToVectorMarketFile(const matrix<double>& A, double B[], vector<double>& X, int num_rows, int num_cols);
void PrintVector(const std::vector<double>& vector);
void PrintMatrix(const matrix<double>& A);
std::vector<double> CreateVectorFilledWithValue(int num_rows);

// Helper function for convergence check
bool convergedUsingL1Norm(const std::vector<double>& X, const std::vector<double>& XPrevious, double tolerance, int size) {
  double distance_max = 0.0;
  for (int i = 0; i < size; ++i) {
    distance_max = std::max(distance_max, std::abs(X[i] - XPrevious[i]));
  }
  return distance_max < tolerance;
}

int Tridiagonal(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols)
{ 
  // Allocate memory for diagonals (consider using smart pointers)
  vector<double> lowerDiagonal(num_cols);
  vector<double> mainDiagonal(num_cols);
  vector<double> upperDiagonal(num_cols);
  // Extract diagonals from matrix A
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      if (i == j) {
        mainDiagonal[i] = A[i][j]; // Use matrix accessor for clarity
      } else if (i == j - 1) {
        upperDiagonal[i] = A[i][j];
      } else if (i == j + 1) {
        lowerDiagonal[i] = A[i][j];
      }
    }
  }
  // Forward elimination
  for (int i = 1; i < num_cols; ++i) {
    double factor = lowerDiagonal[i] / mainDiagonal[i - 1];
    mainDiagonal[i] -= factor * upperDiagonal[i - 1];
    B[i] -= factor * B[i - 1];
  }
  // Backward substitution
  X[num_cols - 1] = B[num_cols - 1] / mainDiagonal[num_cols - 1];
  for (int i = num_cols - 2; i >= 0; --i) {
    X[i] = (B[i] - upperDiagonal[i] * X[i + 1]) / mainDiagonal[i];
  }
  WriteTridiagonalToVectorMarketFile(A,&B[0],X,num_rows,num_cols); 
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
int Tridiagonal_LU_Decomposition( double lowerDiagonal[], double mainDiagonal[],
                                                double upperDiagonal[], int num_cols )
{
   int i;
   for (i = 0; i < (num_cols-1); i++) {
      if (mainDiagonal[i] == 0.0) return -1;
      lowerDiagonal[i] /= mainDiagonal[i];
      mainDiagonal[i+1] -= lowerDiagonal[i] * upperDiagonal[i];
   }
   if (mainDiagonal[num_cols-1] == 0.0) return -1;
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
int Tridiagonal_LU_Solve( double lowerDiagonal[], double mainDiagonal[],
                       double upperDiagonal[], double B[], double x[], int num_cols)
{
   int i;
//         Check that all diagonal elements are nonzero.
//         If a diagonal element is zero then U is singular, so return
//         signalling an error.
   for (i = 0; i < num_cols; i++) if (mainDiagonal[i] == 0.0) return -1;
//         Solve the linear equation Ly = B for y, where L is a lower
//         triangular matrix.
   x[0] = B[0];
   for (i = 1; i < num_cols; i++) x[i] = B[i] - lowerDiagonal[i-1] * x[i-1];
//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Ly = B and U is an upper triangular matrix.
   x[num_cols-1] /= mainDiagonal[num_cols-1];
   for (i = num_cols-2; i >= 0; i--) {
      x[i] -= upperDiagonal[i] * x[i+1];
      x[i] /= mainDiagonal[i];
   }
   return 0;
}

int Jacobi_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter) {
 std::cout << std::endl << "Implementing the Jacobi_L1_Norm method: " << std::endl;
 // Initialize previous iteration vector
 vector<double> XPrevious(num_rows, 0.0);
 // Loop for iterations
 for (int Iter = 1; Iter <= Max_Iter; ++Iter) {
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
   // Check for convergence using L1 norm
   if (convergedUsingL1Norm(X, XPrevious, tolerance, num_cols)) {
     std::cout << "The convergence distance (L1 norm) is less than the tolerance of " << tolerance << "." << endl;
     std::cout << "The number of iterations is " << Iter << "." << endl;
     return 0;  // Convergence achieved
   } 
   // Update previous iteration vector
   XPrevious = X;
 }
 // Maximum iterations reached without convergence
 std::cout << "Maximum number of iterations reached before Jacobi_L1_Distance could finish." << std::endl;
 return 1; // Not converged
}

int Jacobi_L2_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter) {
  std::cout << std::endl << "Implementing the Jacobi_L2_Norm method: " << std::endl;
  // Input validation
  if (num_rows != num_cols || num_rows <= 0) {
    std::cerr << "Error: Invalid matrix dimensions. Matrix must be square and have positive dimensions." << std::endl;
    return -1;
  }
  if (tolerance <= 0.0) {
    std::cerr << "Error: Tolerance must be a positive value." << std::endl;
    return -1;
  }
  if (Max_Iter <= 0) {
    std::cerr << "Error: Maximum iterations must be a positive value." << std::endl;
    return -1;
  }
  // Initialize iteration counter and residual
  int iter = 0;
  double residual = 1.0;
  // Loop until convergence or maximum iterations reached
  while (iter < Max_Iter && residual > tolerance) {
    // Update each element of X using Jacobi_Residual iteration
    residual = 0.0;
    // OpenMP parallel section with reduction
    #pragma omp parallel for reduction(+:residual)
    for (size_t i = 0; i < num_rows; ++i) {
      double updated_value = 0.0;
      for (size_t j = 0; j < num_cols; ++j) {
        if (i != j) {
          updated_value += A[i][j] * X[j];
        }
      }
      updated_value = (B[i] - updated_value) / A[i][i];
      // Check for division by zero (singular matrix)
      //if (std::abs(A[i][i]) < std::numeric_limits<double>::epsilon()) {
      //  std::cerr << "Error: Singular matrix encountered during Jacobi iteration." << std::endl;
      //  return -1;
      //}
      residual += std::pow(X[i] - updated_value, 2.0);
      X[i] = updated_value;
    }
    iter++;
  }
  std::cout << "The convergence distance is less than the tolerance of " << tolerance << "." << std::endl;
  std::cout << "The number of iterations is " << iter << "." << endl;
  // Return number of iterations or -1 if not converged
  return (residual <= tolerance) ? iter : -1;
}

int Gauss_Seidel_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter)
{
 std::cout << endl << "Implementing the Gauss_Seidel_L1_Norm method: " << std::endl; 
 bool convergence_check {false};
 int i {0}, Iter {1};
 double elem {0.0};
 double distance_max {0.0};
 double distance_squared {0.0};
 double norm_variable_squared {0.0};
 vector<double> XPrevious;
 for (int j=0; j<num_cols; j++) {
   XPrevious.push_back(elem);
 }
 std::cout << std::setprecision(5);
 while(convergence_check != true) {   
   distance_max = 0.0;   
   for (i = 0; i < num_rows; i++)
   {
	 // Solve for the i-th element of X (assuming diagonal dominance for stability)
     double updated_value = B[i] / A[i][i];
     for (int j = 0; j < num_rows; ++j) {
       if (j != i) {
         updated_value -= (A[i][j] / A[i][i]) * XPrevious[j];
       }
     }
     X[i] = updated_value;  
     // Update solution and calculate distance
     double distance = std::abs(X[i] - XPrevious[i] );
     distance_max = std::max(distance_max, distance);
     distance_squared = distance*distance;
     // Update norm variables (combine calculations)
     norm_variable_squared += distance*distance;
     XPrevious[i] = X[i];
   }
   double distance_new = sqrt(distance_squared);
   double norm_variable = sqrt(norm_variable_squared);
   double tolerance_variable = distance_new / norm_variable;
   distance_squared = 0.0;
   norm_variable_squared = 0.0;    
   // Check convergence conditions and break the loop if met
  if (norm_variable > 100.0) {
    // Handle large norm scenario
    break;
  } else if (tolerance_variable < tolerance) {
    std::cout << "The convergence distance is less than the tolerance of " << tolerance << "." << endl;
    std::cout << "The number of iterations is " << Iter << "." << endl;
    convergence_check = true;
    break;
  } else if (Iter >= Max_Iter) {
    std::cout << "The distance_max is the following: " << distance_max << endl;
    std::cout << "The tolerance is the following: " << tolerance << endl;
    std::cout << "The number of iterations is " << Iter << "." << endl;
    std::cout << "Maximum number of iterations reached before Jacobi_L1_Distance_Iterative_Solve could finish." << std::endl;
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

int SOR_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, int num_rows, int num_cols, double tolerance, int Max_Iter)
{ 
 std::cout << endl << "Implementing the SOR_L1_Norm method: " << endl;
 bool convergence_check {false};
 int i {0}, j {0}, Iter {1};
 double elem {1.0}, omega {1.25};
 double distance_squared {0.0};
 double norm_variable_squared {0.0}; 
 vector<double> XPrevious;
// Change omega to compare the Gauss-Seidel and SOR methods.
// std::cout << "Please Enter the omega factor for the SOR method." << endl;
// std::cin >> omega;
 for (j=0; j<num_cols; j++) {
   XPrevious.push_back(elem);
 }
 std::cout << std::setprecision(5); 
 while(convergence_check != true) {   
   double distance_max = 0.0;   
   for (i = 0; i < num_rows; i++)
   {
	 // Solve for the i-th element of X (assuming diagonal dominance for stability)
     double updated_value = B[i] / A[i][i];
     for (int j = 0; j < num_rows; ++j) {
       if (j != i) {
         updated_value -= (A[i][j] / A[i][i]) * XPrevious[j];
       }
     }
     X[i] = updated_value;  
     // SOR factor omega calculated on the next line  
     X[i] = omega * X[i] + (1 - omega) * XPrevious[i];
     // Update solution and calculate distance
     double distance = std::abs(X[i] - XPrevious[i]);
     distance_max = std::max(distance_max, distance);
     // Update norm variables (combine calculations)
     distance_squared = distance*distance;
     norm_variable_squared += distance_squared;
     XPrevious[i] = X[i];
   }
   double distance_new = sqrt(distance_squared);
   double norm_variable = sqrt(norm_variable_squared);
   double tolerance_variable = distance_new / norm_variable;
   distance_squared = 0.0;
   norm_variable_squared = 0.0; 
  // Check convergence conditions and break the loop if met
  if (norm_variable > 100.0) {
   // Handle large norm scenario
   break;
  } else if (tolerance_variable < tolerance) {
    std::cout << "The convergence distance is less than the tolerance of " << tolerance << "." << endl;
    std::cout << "The number of iterations is " << Iter << "." << endl;
    convergence_check = true;
    break;
  } else if (Iter >= Max_Iter) {
    std::cout << "The distance_max is the following: " << distance_max << endl;
    std::cout << "The tolerance is the following: " << tolerance << endl;
    std::cout << "The number of iterations is " << Iter << "." << endl;
    std::cout << "Maximum number of iterations reached before Jacobi_L1_Distance_Iterative_Solve could finish." << std::endl;
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

int ReadMatrixMarketMatrix(matrix<double>& A, int& num_rows, int& num_cols) {

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
 if (fscanf(myFile, "%d %d %d\n", &num_rows, &num_cols, &number_of_entries_A) != 3) {
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
   if (fscanf(myFile, "%d %d %lf\n", &i_index, &j_index, &elem) != 3) {
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

// ReadMatrixMarketVector implementation for vectors
int ReadMatrixMarketVector(vector<double>& B, int& num_cols) {
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
     fscanf (myFile2, "%d,", &dim_B_Array[i]);
   }
   for (int i = 0; i < dim_B_Array[1]; i++)
     B.push_back(0.0);
   for (int i = 0; i < dim_B_Array[1]; i++)
   {
     fscanf (myFile2, "%d,", &i_index);
     i_index--;
     fscanf (myFile2, "%d,", &j_index);
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

int Write1DArrayToMatrixMarketFile(const double B[], int num_rows) {  
  //Use C++ streams for safer file handling
  ofstream outfile("B_out.dat");
  if (!outfile.is_open()) {
    cerr << "Error opening file for writing: " << "B_out.dat" << endl;
    return 1;
  }
  printf ("B =\n");
  // Write header information (assuming general coordinate pattern)
  outfile << "%%MatrixMarket_Output_vector_B.dat matrix coordinate pattern general\n";
  outfile << num_rows << " 1 " << num_rows << endl; // Adjust for 1D array
  // Write each element with row and column indices (starting from 1)
  for (int i = 0; i < num_rows; i++) {
    outfile << i + 1 << " " << 1 << " " << B[i] << endl;
    printf ("%6.5lf    ", B[i]);
  }
  std::cout << std::endl;
  outfile.close();
  return 0;
}

int Write2DArrayToMatrixMarketFile(const double array_A[], int num_rows, int num_cols) {
  // Use C++ streams for safer file handling
  ofstream outfile("A_out.dat");
  if (!outfile.is_open()) {
    cerr << "Error opening file for writing: A_out.dat" << endl;
    return 1;
  }
  // Write header information (assuming general coordinate pattern)
  outfile << "%%MatrixMarket_Output_vector_A.dat matrix coordinate pattern general\n";
  outfile << num_rows << " " << num_cols << " " << num_rows * num_cols << endl;
  printf ("A =\n");
  // Write each element with row and column indices (starting from 1)
  for (int i = 0; i < num_rows; i++) {
    for (int j = 0; j < num_cols; j++) {
      outfile << i + 1 << " " << j + 1 << " " << array_A[i * num_cols + j] << endl;
      printf ("%6.5lf    ", array_A[i * num_cols + j]);
    }
    std::cout << std::endl;
  }
  outfile.close();
  return 0;
}

int WriteMatrixMarketMatrix(const matrix<double>& A, int num_rows, int num_cols) {
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
 for (int i = 0; i < num_rows; ++i) {
   for (int j = 0; j < num_cols; ++j) {
     if (A[i][j] != 0.0) 
     { // Check for non-zero value
       fprintf(fptr, "%d %d %lf\n", i + 1, j + 1, A[i][j]);
     }
   }
 }
 // Close the file
 fclose(fptr);
 return 0; // Indicate success
}

int WriteMatrixMarketVector(const std::vector<double>& X, int num_rows, int num_cols) {
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
 for (size_t i = 0; i < num_rows; ++i) {
   if (X[i] != 0.0) 
   { // Check for non-zero value
     fprintf(fptr, "%d %zu %lf\n", 1, i+1, X[i]); // Row index always 1 for a vector
     std::cout << 1 << "   " << i+1 << "   "<< X[i] << std::endl;
   }
 }
 // Close the file
 fclose(fptr);
 return 0; // Indicate success
}

int WriteTridiagonalToVectorMarketFile(const matrix<double>& A, double B[], vector<double>& X, int num_rows, int num_cols)
{
  std::cout << endl << "Implementing the Tridiagonal method: " << std::endl;
  printf ("A =\n");
  PrintMatrix(A);  
  Write1DArrayToMatrixMarketFile(&B[0], num_rows);   
  printf ("X =\n");
  PrintVector(X);
  WriteMatrixMarketVector(X, num_rows, num_cols);
  return 0;
}

void PrintVector(const std::vector<double>& vector) {
// std::cout << "Displaying the vector: " << endl;
 std::cout << std::setprecision(5);
 for (const double& value : vector) {
   std::cout << value << " ";
 }
 std::cout << std::endl;
}

void PrintMatrix(const matrix<double>& A) {
// std::cout << "Displaying the 2D vector:" << endl;
 for (const auto& row : A) {
   for (const double& value : row) {
     std::cout << value << " ";
   }
   std::cout << std::endl;
 }
}

std::vector<double> CreateVectorFilledWithValue(int num_rows) {
  return std::vector<double>(num_rows, 0.0);
}

int main ()
{
   int num_rows {0};
   int num_cols {0};
   vector <double> B;
   ReadMatrixMarketVector(B, num_cols);
   matrix<double> A;
   ReadMatrixMarketMatrix(A,num_rows,num_cols);         
   std::cout << std::setprecision(5);
   int Max_Iter {1000};
   double tolerance {0.00001};
   double omega {1.25};
   (void) omega;
   std::vector<double> X = CreateVectorFilledWithValue(num_rows);
   // Using time point and system_clock
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();
   Tridiagonal(A, B, X, num_rows, num_cols);
   for (size_t j = 0; j < num_cols; j++) {
     X[j] = 0.0;
   }
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end - start;
   std::time_t end_time = std::chrono::system_clock::to_time_t(end);
   std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
   ReadMatrixMarketVector(B, num_cols);  
   // Using time point and system_clock
   start = std::chrono::system_clock::now();
   Jacobi_L1_Norm(A, B, X, num_rows, num_cols, tolerance, Max_Iter);
   end = std::chrono::system_clock::now();
   elapsed_seconds = end - start;
   std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n"; 
   printf ("X =\n");
   PrintVector(X);
   // Using time point and system_clock
   start = std::chrono::system_clock::now();           
   Jacobi_L2_Norm(A, B, X, num_rows, num_cols, tolerance, Max_Iter);
   end = std::chrono::system_clock::now();
   elapsed_seconds = end - start;
   std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";    
   printf ("X =\n");
   PrintVector(X);
   // Using time point and system_clock
   start = std::chrono::system_clock::now();  
   Gauss_Seidel_L1_Norm(A, B, X, num_rows, num_cols, tolerance, Max_Iter);
   end = std::chrono::system_clock::now();
   elapsed_seconds = end - start;
   std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";   
   printf ("X =\n");
   PrintVector(X);
   // Using time point and system_clock
   start = std::chrono::system_clock::now();  
   SOR_L1_Norm(A, B, X, num_rows, num_cols, tolerance, Max_Iter);
   end = std::chrono::system_clock::now();
   elapsed_seconds = end - start;
   std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";   
   printf ("X =\n");
   PrintVector(X);
   //   WriteMatrixMarketVector(X,num_rows,num_cols);
   //   WriteMatrixMarketMatrix(A,num_rows,num_cols);
   A.clear();
   A.shrink_to_fit();
   B.clear();
   B.shrink_to_fit();
   X.clear();
   X.shrink_to_fit();
   return 0;
}
