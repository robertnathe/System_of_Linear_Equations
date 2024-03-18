#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std; 
#define a(i, j) A[i][j] // Assuming A is a matrix (use appropriate access method)
template<typename T>
using matrix = std::vector<std::vector<T>>;

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
using namespace std;

#define a(i, j) A[i][j] // Assuming A is a matrix (use appropriate access method)

template<typename T>
using matrix = std::vector<std::vector<T>>;

int SOR(const matrix<double>& A, vector<double>& B, vector<double>& X, unsigned int num_rows, unsigned int num_cols, double tolerance, unsigned int Max_Iter);
int Gauss_Seidel(const matrix<double>& A, vector<double>& B, vector<double>& X, unsigned int num_rows, unsigned int num_cols, double tolerance, unsigned int Max_Iter);
int Jacobi(const matrix<double>& A, vector<double>& B, vector<double>& X, unsigned int num_rows, unsigned int num_cols, double tolerance, unsigned int Max_Iter);
int ReadMatrixMM(matrix<double> & A, unsigned int& num_rows, unsigned int& num_cols); // For matrices 
int ReadVectorMM(vector<double>& B, unsigned int& num_cols); // For vectors
int WriteMatrixMM(const matrix<double>& A, unsigned int num_rows, unsigned int num_cols);
int WriteVectorMM(const std::vector<double>& X, unsigned int num_rows, unsigned int num_cols);
void PrintVector(const std::vector<double>& vector);
void PrintVector2D(const matrix<double>& A);
std::vector<double> CreateZeroVector(unsigned int num_rows);

template<typename T>
using matrix = std::vector< std::vector<T> >;

int SOR (const matrix<double>& A, vector<double>& B, vector<double>& X, unsigned int num_rows, unsigned int num_cols, double tolerance, unsigned int Max_Iter)
{ 
 std::cout << endl << "Implementing the SOR method: " << endl;
 bool convergence_check {false};
 unsigned int i {0}, j {0}, Iter {1};
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

int Gauss_Seidel (const matrix<double>& A, vector<double>& B, vector<double>& X, unsigned int num_rows, unsigned int num_cols, double tolerance, unsigned int Max_Iter)
{
 std::cout << endl << "Implementing the Gauss Seidel method: " << endl; 
 bool convergence_check {false};
 unsigned int i {0}, j {0}, Iter {1};
 double elem {0.0}, distance {0.0}, distance_max {0.0};
 double temporary_variable_one {0.0}, temporary_variable_two {0.0};
 double distance_squared {0.0};
 double distance_new {0.0};
 double norm_variable_squared {0.0};
 double norm_variable {0.0};
 double tolerance_variable {0.0};
 vector<double> XPrevious;
 for (unsigned int j=0; j<num_cols; j++) {
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

int Jacobi(const matrix<double>& A, vector<double>& B, vector<double>& X, unsigned int num_rows, unsigned int num_cols, double tolerance, unsigned int Max_Iter) {

 std::cout << std::endl << "Implementing the Jacobi method: " << endl;
 // Initialize previous iteration vector
 vector<double> XPrevious(num_rows, 0.0);
 // Loop for iterations
 for (unsigned int Iter = 1; Iter <= Max_Iter; ++Iter) {
   double distance_max = 0.0;
   // Update each element of X
   for (unsigned int i = 0; i < num_rows; ++i) {
     double sum = 0.0;
     for (unsigned int j = 0; j < num_cols; ++j) {
       if (i != j) {
         sum += A[i][j] * XPrevious[j];
       }
     }
     X[i] = (B[i] - sum) / A[i][i];
   }
   // Calculate distance metric and check for convergence
   for (unsigned int j = 0; j < num_cols; ++j) {
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

int ReadMatrixMM(matrix<double>& A, unsigned int& num_rows, unsigned int& num_cols) {

 unsigned int number_of_entries_A {0};
 unsigned int i_index {0}, j_index {0};
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
 for (unsigned int i = 0; i < num_rows; ++i) {
   A[i].resize(num_cols);
 }
 // Read non-zero elements by row and column indices
 for (unsigned int i = 0; i < number_of_entries_A; ++i) {
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
int ReadVectorMM(vector<double>& B, unsigned int& num_cols) {
    (void) num_cols; 
   FILE *myFile2;
   myFile2 = fopen ("B.dat", "r");
   unsigned int dim_B_Array[3];
   unsigned int i_index {0}, j_index {0};
   double value {0.0};
   while (myFile2 == NULL)
   {
    std::cout << "Error Reading File" << endl;
     exit (0);
   } 
   fscanf (myFile2, "%*s %*s %*s %*s %*s");
   for (unsigned int i = 0; i < 3; i++)
   {
     fscanf (myFile2, "%u,", &dim_B_Array[i]);
   }
   for (unsigned int i = 0; i < dim_B_Array[1]; i++)
     B.push_back(0.0);
   for (unsigned int i = 0; i < dim_B_Array[1]; i++)
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

int WriteMatrixMM(const matrix<double>& A, unsigned int num_rows, unsigned int num_cols) {
 // Open file in writing mode
 unsigned int CountNonZeroEntries {0};
 FILE* fptr = fopen("A_out.txt", "w");
 if (fptr == NULL) {
   std::cerr << "Error opening file for writing matrix A." << std::endl;
   return -1; // Indicate error
 }
 for (unsigned int i = 0; i < num_rows; ++i) {
   for (unsigned int j = 0; j < num_cols; ++j) {
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
 for (unsigned int i = 0; i < num_rows; ++i) {
   for (unsigned int j = 0; j < num_cols; ++j) {
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

int WriteVectorMM(const std::vector<double>& X, unsigned int num_rows, unsigned int num_cols) {
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
 for (unsigned int i = 0; i < num_rows; ++i) {
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

std::vector<double> CreateZeroVector(unsigned int num_rows) {
 return std::vector<double>(num_rows);
}

int main ()
{
   unsigned int num_rows {0}, num_cols {0};
   vector <double> B;
   ReadVectorMM(B, num_cols);
   matrix<double> A;
   ReadMatrixMM(A,num_rows,num_cols);         
   std::cout << std::setprecision(5);
   unsigned int Max_Iter {1000};
   double tolerance {0.00001};
   vector<double> X = CreateZeroVector(num_rows);
   Jacobi (A, B, X, num_rows, num_cols, tolerance, Max_Iter);
//   WriteVectorMM(X,num_rows,num_cols);
   PrintVector2D(A);
   PrintVector(B);
   PrintVector(X);
   Gauss_Seidel (A, B, X, num_rows, num_cols, tolerance, Max_Iter);
//   WriteVectorMM(X,num_rows,num_cols);
//   WriteMatrixMM(A,num_rows,num_cols);
   PrintVector2D(A);
   PrintVector(B);
   PrintVector(X);
   SOR (A, B, X, num_rows, num_cols, tolerance, Max_Iter);
   PrintVector2D(A);
   PrintVector(B);
   PrintVector(X);
   A.clear();
   A.shrink_to_fit();
   B.clear();
   B.shrink_to_fit();
   X.clear();
   X.shrink_to_fit();
   return 0;
}
