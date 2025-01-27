#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <sstream>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::ostringstream;
using std::tuple;
using std::get;

// Type definition for 2D Matrix
template<typename T>
using Matrix = std::vector<std::vector<T>>;

inline double& a(Matrix<double>& mat, size_t i, size_t j) {
    return mat[i][j];
}

// Function Prototypes
void PrintVector2D(const Matrix<double>& TempMatrix);
void PrintVector(const double TempVector[], unsigned int nrows);
int Output_Results(double X[], double B[], unsigned int nrows, unsigned int ncols);
int ReadMatrixMarketMatrix(Matrix<double>& A, size_t& num_rows, size_t& num_cols);
int ReadMatrixMarketVector(vector<double>& B, size_t& num_rows);

// Other function prototypes...
int Output_Data(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter);
int Write1DArrayToMatrixMarketFile(const double B[], size_t num_rows);
int Write2DArrayToMatrixMarketFile(const double array_A[], size_t num_rows, size_t num_cols);
int WriteMatrixMarketMatrix(const Matrix<double>& A, size_t num_rows, size_t num_cols);
int WriteMatrixMarketVector(const vector<double>& X, size_t num_cols);
int WriteTridiagonalToVectorMarketFile(double B[], vector<double>& X, size_t num_rows);
void PrintExecutionTime(const std::chrono::time_point<std::chrono::system_clock>& start,
                        const std::chrono::time_point<std::chrono::system_clock>& end);
void PrintVector(const string& label, const vector<double>& vec);
void PrintVector(const string& label, const double vec[], unsigned int n);
void PrintVector(const string& label, const double vec[], size_t size);
void PrintMatrix(const string& label, const Matrix<double>& A);
void initializeTestMatrices(Matrix<double>& A, Matrix<double>& AMat, Matrix<double>& BMat, Matrix<double>& CMat, size_t nrows, size_t ncols);
void initializeTestMatricesWithOffset(Matrix<double>& AMat, Matrix<double>& BMat, size_t nrows, size_t ncols);
void Identity_Matrix_ut(Matrix<double>& A, unsigned int n);
double Sum_over_Rows(const Matrix<double>& A, double X[], unsigned int nrows);
double Sum_over_Column(const Matrix<double>& A, unsigned int col);
double Sum_over_Row(const Matrix<double>& A, unsigned int row);
double Trace_of_Matrix(const Matrix<double>& A);
int Testing_General_Purpose_Parameters(Matrix<double>& A, size_t nrows, size_t ncols);
void Copy_Matrix(vector<double> A[], vector<double> B[], unsigned int nrows, unsigned int ncols);
void Copy_Vector(double *d, double *s, unsigned int ncols);
void Get_Row(double X[], vector <double> A[], unsigned int kth_row, unsigned int ncols);
void Get_Column(double X[], vector <double> A[], unsigned int kth_column);
void Set_Row(double X[], vector <double> A[], unsigned int kth_row, unsigned int ncols);
void Set_Column(double X[], vector <double> A[], unsigned int kth_column, unsigned int nrows);
void Get_Diagonal(double X[], vector <double> A[], unsigned int nrows, unsigned int ncols);
void Set_Diagonal(double X[], vector <double> A[], unsigned int nrows, unsigned int ncols);
void Set_Diagonal_to_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols);
void Get_SubMatrix(vector<double> A[], vector<double> S[], unsigned int mrows, unsigned int mcols, unsigned int row, unsigned int col);
void Set_SubMatrix(vector<double> A[], unsigned int mrows, unsigned int mcols, unsigned int row, unsigned int col, double Scalar);
void Fill_Matrix_with_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols);
int Join_Rows(vector<double> A[], vector<double> BMatrix[], vector <double> CMatrix[],  unsigned int nrows, unsigned int ncols, unsigned int mcols);
int Join_Columns(vector<double> A[], vector<double> BMatrix[], vector <double> CMatrix[],  unsigned int nrows, unsigned int ncols, unsigned int mcols);
void Interchange_Rows(Matrix<double>& A, unsigned int row1, unsigned int row2, unsigned int ncols);
void Interchange_Columns(vector<double> A[], unsigned int col1, unsigned int col2, unsigned int nrows);
int Vec_Swap(double MyVector1[], double MyVector2[], unsigned int Vec_Size);
int Test_Vec_Swap(double MyVector1[], double MyVector2[], unsigned int Vec_Size);

// Main function
int main() {
    double scalarValue {10.0}; // Example scalar value
    (void) scalarValue;
    size_t num_rows = 0, num_cols = 0;
    Matrix<double> A;
    vector<double> B, X;
    // Read Matrix A
    if (ReadMatrixMarketMatrix(A, num_rows, num_cols) != 0) {
        std::cerr << "Failed to read Matrix A." << std::endl;
        return -1;
    }
    // Read vector B
    size_t vector_rows = 0;
    if (ReadMatrixMarketVector(B, vector_rows) != 0) {
        std::cerr << "Failed to read vector B." << std::endl;
        return -1;
    }
    if (num_rows != vector_rows) {
        std::cerr << "Error: Matrix and vector dimensions do not match." << std::endl;
        return -1;
    }
    // Initialize vector X
    X.resize(num_rows, 0.0);
    // Begin implementing your code here.
    
    
    
    // End implementing your code here.
    // Print matrices and vectors
    //WriteMatrixMarketMatrix(A, num_rows, num_cols);
	//WriteMatrixMarketVector(B, num_cols);
    PrintMatrix("A", A);
    PrintVector("B", B);
    PrintVector("X", X);
    // Cleanup
    A.clear();
    B.clear();
    X.clear();
    return 0;
}

// Function Definitions...

void PrintVector2D(const Matrix<double>& TempMatrix) {
    std::cout << "Displaying the 2D vector:" << std::endl;
    std::cout << std::setprecision(5);
    for (const auto& row : TempMatrix) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

void Identity_Matrix_ut(Matrix<double>& A, unsigned int n) {
    A.resize(n, vector<double>(n, 0.0));
    for (unsigned int i=0; i<n; i++) A[i][i] = 1.0;
}

// Sum of all elements in each row
double Sum_over_Rows(const Matrix<double>& A, double X[], unsigned int nrows) {
    for (unsigned int i=0; i<nrows; i++) {
        X[i] = 0.0;
        for (unsigned int j=0; j<A[i].size(); j++) 
            X[i] += A[i][j];
    }
    return 0;
}

// Sum of elements in a specific column
double Sum_over_Column(const Matrix<double>& A, unsigned int col) {
    double sum = 0.0;
    for (unsigned int i=0; i<A.size(); i++)
        if (col < A[i].size()) sum += A[i][col];
    return sum;
}

// Sum of elements in a specific row
double Sum_over_Row(const Matrix<double>& A, unsigned int row) {
    double sum = 0.0;
    if (row < A.size()) 
        for (double val : A[row]) sum += val;
    return sum;
}

double Trace_of_Matrix(const Matrix<double>& A) {
    double trace = 0.0;
    for (unsigned int i=0; i<A.size(); i++)
        if (i < A[i].size()) trace += A[i][i];
    return trace;
}

int Testing_General_Purpose_Parameters(Matrix<double>& A, size_t nrows, size_t ncols) {
    Matrix<double> A_Copy(nrows, std::vector<double>(ncols));
    for (unsigned int i = 0; i < nrows; ++i) {
        for (unsigned int j = 0; j < ncols; ++j) {
            A_Copy[i][j] = A[i][j];
        }
    }
    return 0;
}

void Copy_Matrix(vector<double> A[], vector<double> B[], unsigned int nrows, unsigned int ncols)
{
  unsigned int i,j;
  for (i = 0; i < nrows; i ++ )
  {
    for (j = 0; j < ncols; j ++ )
    {
      B[i][j] = A[i][j];
    }   
  }
}

 extern void Copy_Vector(double *d, double *s, unsigned int ncols)
{
   memcpy(d, s, sizeof(double) * ncols);
}

void Get_Row(double X[], vector <double> A[], unsigned int kth_row, unsigned int ncols)
{
  unsigned int j;
  for (j = 0; j < ncols; j++)
    X[j] = A[kth_row][j];
}

void Get_Column(double X[], vector <double> A[], unsigned int kth_column, unsigned int nrows)
{
  unsigned int i;
  for (i = 0; i < nrows; i++)
    X[i] = A[i][kth_column];
}

void Set_Row(double X[], vector <double> A[], unsigned int kth_row, unsigned int ncols)
{
  unsigned int j;
  for (j = 0; j < ncols; j++)
    A[kth_row][j] = X[j];
}

void Set_Column(double X[], vector <double> A[], unsigned int kth_column, unsigned int nrows)
{
  unsigned int i;
  for (i = 0; i < nrows; i++)
    A[i][kth_column] = X[i];
}

void Get_Diagonal(double X[], vector <double> A[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i, j, n;
   n = (nrows < ncols) ? nrows: ncols;
   for (i=0; i < n; i++)
     for (j=0; j < n; j++)
       X[i] = A[i][i];
}

void Set_Diagonal(double X[], vector <double> A[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i, j, n;
   if (nrows < ncols) n = nrows; else n = ncols;
   for (i=0; i < n; i++)
     for (j=0; j < n; j++)
       A[i][i] = X[i];
}

void Set_Diagonal_to_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols)
{
   unsigned int i, n  = (nrows < ncols) ? nrows : ncols;
   for (i = 0; i < n; i++)
     A[i][i] = Scalar;
}

void Get_SubMatrix(vector<double> A[], vector<double> S[], unsigned int mrows, unsigned int mcols, unsigned int row, unsigned int col)
{
  unsigned int i {0}, j {0};
  if (row == 0)
  {
    for(i = 0; i < (mrows-row); i++)
    {
      for (j=0; j < (mcols-col); j++)
      {
        S[i][j] = A[row+i][col+j];
      }
    }
  }
  if (row != 0)
  {  
    for (i = 0; i <= (mrows-row); i++) 
    {
      for (j = 0; j <= (mcols-col); j++)
      {
        S[i][j] = A[row+i][col+j];
      }
    }
  }
}

void Set_SubMatrix(vector<double> A[], unsigned int mrows, unsigned int mcols, unsigned int row, unsigned int col, double Scalar)
{
  unsigned int i {0}, j {0};
  if (row == 0)
  {
    for(i = 0; i < (mrows-row); i++)
    {
      for (j=0; j < (mcols-col); j++)
      {
		A[row+i][col+j] = Scalar;
      }
    }
  }
  if (row != 0)
  {  
    for (i = 0; i <= (mrows-row); i++) 
    {
      for (j = 0; j <= (mcols-col); j++)
      {
		  A[row+i][col+j] = Scalar;
      }
    }
  }
}

void Fill_Matrix_with_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols)
{
   unsigned int i,j;
   for (i = 0; i < nrows; i++)
      for ( j = 0; j < ncols; j++) A[i][j] = Scalar;
}

int Join_Rows(vector<double> A[], vector<double> BMatrix[], vector <double> CMatrix[],  unsigned int nrows, unsigned int ncols, unsigned int mcols)
{
  unsigned int i, j;
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols+mcols; j++)
      if (j < ncols) { 
        CMatrix[i][j] = A[i][j];
      }
      else {
	    CMatrix[i][j] = BMatrix[i][j-ncols];	  
	  }
  }
  return 0;	
}

int Join_Columns(vector<double> A[], vector<double> BMatrix[], vector <double> CMatrix[],  unsigned int nrows, unsigned int ncols, unsigned int mcols)
{
  unsigned int i, j;
  for (i = 0; i < (nrows+mcols); i++) {
    for (j = 0; j < ncols; j++)
      if (i < nrows) { 
        CMatrix[i][j] = A[i][j];
      }
      else {
	    CMatrix[i][j] = BMatrix[i-nrows][j];	  
	  }
  }
  return 0;	
}

void Interchange_Rows(Matrix<double>& A, unsigned int row1, unsigned int row2, unsigned int ncols) {
    unsigned int j;
    double temp;
    for (j = 0; j < ncols; j++) {
        temp = a(A, row1, j);
        a(A, row1, j) = a(A, row2, j);
        a(A, row2, j) = temp;
    }
}

void Interchange_Columns(vector<double> A[], unsigned int col1, unsigned int col2, unsigned int nrows)
{
   unsigned int i;
   double temp;
   for (i = 0; i < nrows; i++)
   {
     temp = A[i][col1];
     A[i][col1] = A[i][col2];
     A[i][col2] = temp;  
   }
}

int Vec_Swap(double MyVector1[], double MyVector2[], unsigned int Vec_Size)
{
  unsigned int i;
  double TempVar {0.0};
  for (i = 0; i < Vec_Size; i++)
  {
    TempVar = MyVector2[i];
	MyVector2[i] = MyVector1[i];
	MyVector1[i] = TempVar;
  }
  return 0;
}

int Test_Vec_Swap(double MyVector1[], double MyVector2[], unsigned int Vec_Size)
{
  Vec_Swap(&MyVector1[0], &MyVector2[0], Vec_Size); 
  return 0;	
}

void PrintVector(const double TempVector[], unsigned int nrows) {
    std::cout << "Displaying the vector: " << std::endl;
    std::cout << std::setprecision(5);
    for (unsigned int i = 0; i < nrows; i++)
        std::cout << TempVector[i] << "   "; // Changed spacing for readability
    std::cout << std::endl;
}

int Output_Results(double X[], double B[], unsigned int nrows, unsigned int ncols)
{
  unsigned int i;
  std::cout << std::setprecision (5) << std::endl;
  std::cout << "Displaying Output_Results" << std::endl;
  std::cout << "The vector X is the following: " << std::endl;
  PrintVector(&X[0], nrows);
  std::cout << "The vector B is the following: " << std::endl;
  PrintVector(&B[0], ncols);
  // Create a new file named "C.dat"
  std::ofstream outputFile("C.dat");
  if (outputFile.is_open()) {
    // Write some text into the file
    outputFile << "%%MatrixMarket_Output_vector_C.dat Matrix coordinate pattern general" << std::endl;
    outputFile << 1 << " " << nrows << " " << nrows;
    outputFile << std::endl;
    for (i = 0; i < ncols; i++)
      outputFile <<  1 << " " << i+1 << " " << X[i] << std::endl;
    // Close the file
    outputFile.close();
  } else {
    std::cout << "Error writing to file." << std::endl;
  }
  return 0;
}

int Output_Data(Matrix<double>& A, double X[], double B[], unsigned int nrows, unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter) {
    (void) max_iter;
    (void) tolerance;
    (void) eigenvalue;
    std::cout << std::setprecision(7) << std::endl;
    std::cout << "******************** Solve Ax = B ********************" << std::endl;
    std::cout << "The Matrix A is the following: " << std::endl;
    PrintVector2D(A);
    std::cout << "The vector X is the following: " << std::endl;
    PrintVector(X, ncols);
    std::cout << "The vector B is the following: " << std::endl;
    PrintVector(B, ncols);
    std::ofstream outputFile("C.dat");
    if (outputFile.is_open()) {
        outputFile << "%%MatrixMarket_Output_vector_C.dat Matrix coordinate pattern general" << std::endl;
        outputFile << 1 << " " << ncols << " " << nrows << std::endl; // Corrected output dimensions
        for (unsigned int i = 0; i < ncols; i++) { // Declared i here
            outputFile << 1 << " " << i + 1 << " " << X[i] << std::endl;
        }
        outputFile.close();
    } else {
        std::cout << "Error writing to file." << std::endl;
    }
    return 0;
}

int Write1DArrayToMatrixMarketFile(const double B[], size_t num_rows) {  
    ofstream outfile("B_out.dat");
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: B_new_out.dat" << std::endl;
        return 1;
    }
    outfile << "%%MatrixMarket Matrix coordinate real general\n";
    outfile << 1 << " " << num_rows << " " << num_rows << "\n";
    outfile << std::fixed << std::setprecision(6);
    //std::cout << "B:\n";
    for (size_t i = 0; i < num_rows; ++i) {
        outfile << 1 << " " << (i + 1) << " " << B[i] << "\n";
        //printf("%7.6lf    ", B[i]);
    }
    //std::cout << "\n";
    outfile.close();
    return 0;
}

int Write2DArrayToMatrixMarketFile(const double array_A[], size_t num_rows, size_t num_cols) {
    ofstream outfile("A_out.dat");
    if (!outfile.is_open()) {
        std::cerr << "Error opening file for writing: A_out.dat" << std::endl;
        return 1;
    }
    outfile << "%%%%MatrixMarket Matrix coordinate real general\n"
            << num_rows << " " << num_cols << " " << num_rows * num_cols << "\n";    
    ostringstream file_buffer;
    file_buffer << std::scientific << std::setprecision(6);
    std::cout << "A:\n";
    for (size_t i = 0; i < num_rows; ++i) {
        ostringstream console_line;
        console_line << std::fixed << std::setprecision(5);
        for (size_t j = 0; j < num_cols; ++j) {
            const double val = array_A[i * num_cols + j];
            file_buffer << i+1 << " " << j+1 << " " << val << "\n";   
            console_line << std::setw(6) << val;
            if (j < num_cols-1) console_line << "    ";
        }
        console_line << "\n";
        std::cout << console_line.str();
    }
    outfile << file_buffer.str();
    outfile.close();
    return 0;
}

int WriteMatrixMarketMatrix(const Matrix<double>& A, size_t num_rows, size_t num_cols) {
    FILE* fptr = fopen("A_out.dat", "w");
    if (!fptr) return -1;
    vector<tuple<size_t, size_t, double>> entries;
    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_cols; ++j) {
            if (A[i][j] != 0.0) {
                entries.emplace_back(i+1, j+1, A[i][j]);
            }
        }
    }
    fprintf(fptr, "%%%%MatrixMarket Matrix coordinate real general\n");
    fprintf(fptr, "%zu %zu %zu\n", num_rows, num_cols, entries.size());
    for (const auto& entry : entries) {
        fprintf(fptr, "%zu %zu %.15g\n", 
                get<0>(entry), 
                get<1>(entry), 
                get<2>(entry));
    }
    fclose(fptr);
    return 0;
}

int WriteMatrixMarketVector(const vector<double>& X, size_t num_cols) {
    FILE* fptr = fopen("X_out.dat", "w");
    if (!fptr) return -1;
    fprintf(fptr, "%%%%MatrixMarket Matrix coordinate real general\n");
    fprintf(fptr, "1 %zu %zu\n", num_cols, X.size());
    for (size_t i = 0; i < X.size(); ++i) {
        fprintf(fptr, "1 %zu %.15g\n", i + 1, X[i]);
    }
    fclose(fptr);
    return 0;
}

int WriteTridiagonalToVectorMarketFile(double B[], vector<double>& X, size_t num_rows) {
    Write1DArrayToMatrixMarketFile(B, num_rows);
    WriteMatrixMarketVector(X, X.size());
    return 0;
}

int ReadMatrixMarketMatrix(Matrix<double>& A, size_t& num_rows, size_t& num_cols) {
    ifstream file("A.dat");
    if (!file.is_open()) return -1;
    string line;
    while (getline(file, line) && line[0] == '%');
    size_t num_entries;
    istringstream header_stream(line);
    if (!(header_stream >> num_rows >> num_cols >> num_entries)) return -1;
    A.resize(num_rows, vector<double>(num_cols, 0.0));
    for (size_t i = 0; i < num_entries; ++i) {
        size_t row, col;
        double value;
        if (!(file >> row >> col >> value)) return -1;
        row--; col--;
        A[row][col] = value;
    }
    return 0;
}

int ReadMatrixMarketVector(vector<double>& B, size_t& num_rows) {
    ifstream file("B.dat");
    if (!file.is_open()) {
        std::cerr << "Error opening vector file B.dat" << std::endl;
        return -1;
    }
    string line;
    while (getline(file, line) && line[0] == '%'); // Skip comments
    size_t num_entries;
    istringstream header_stream(line);
    size_t declared_rows, declared_cols;
    if (!(header_stream >> declared_rows >> declared_cols >> num_entries)) return -1;
    // Check if it's a single-row matrix
    if (declared_rows == 1 && declared_cols > 1) {
        num_rows = declared_cols; // Treat as a vector
        B.resize(num_rows);
        for (size_t i = 0; i < num_entries; ++i) {
            size_t row, col;
            double value;
            if (!(file >> row >> col >> value)) return -1;
            B[col - 1] = value; // Fill vector
        }
    } else {
        std::cerr << "Error: Expected a single row for vector." << std::endl;
        return -1;
    }   
    return 0;
}


void PrintExecutionTime(const std::chrono::time_point<std::chrono::system_clock>& start,
                        const std::chrono::time_point<std::chrono::system_clock>& end) {
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << "s\n";
}

void PrintVector(const string& label, const vector<double>& vec) {
    std::cout << label << ":\n[";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << std::fixed << std::setprecision(5) << vec[i];
        if (i != vec.size()-1) std::cout << ", ";
    }
    std::cout << "]\n" << std::endl;
}

void PrintVector(const string& label, const double vec[], unsigned int n) {
    std::cout << label << ":\n[";
    for (unsigned int i = 0; i < n; ++i) {
        std::cout << std::fixed << std::setprecision(5) << vec[i];
        if (i != n - 1) std::cout << ", ";
    }
    std::cout << "]\n" << std::endl;
}

void PrintVector(const string& label, const double vec[], size_t size) {
    std::cout << label << ": [";
    for (size_t i = 0; i < size; ++i) {
        std::cout << std::fixed << std::setprecision(5) << vec[i];
        if (i < size - 1) std::cout << ", ";
    }
    std::cout << "]\n";
}

void PrintMatrix(const string& label, const Matrix<double>& A) {
    std::cout << label << ":\n";
    for (const auto& row : A) {
        for (double value : row) {
            std::cout << std::fixed << std::setprecision(5) << std::setw(10) << value;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void initializeTestMatrices(Matrix<double>& A, Matrix<double>& AMat, Matrix<double>& BMat, Matrix<double>& CMat, size_t nrows, size_t ncols) {
    A.assign(nrows, std::vector<double>(ncols, 0.0));
    AMat.assign(nrows, std::vector<double>(ncols));
    BMat.assign(nrows, std::vector<double>(ncols));
    CMat.assign(nrows, std::vector<double>(ncols));
    for (size_t i = 0; i < nrows; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            AMat[i][j] = static_cast<double>(i * nrows + j);
            BMat[i][j] = static_cast<double>(i * nrows + j);
        }
    }
}

void initializeTestMatricesWithOffset(Matrix<double>& AMat, Matrix<double>& BMat, size_t nrows, size_t ncols) {
    AMat.assign(nrows, std::vector<double>(ncols));
    BMat.assign(nrows, std::vector<double>(ncols));
    for (size_t i = 0; i < nrows; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            AMat[i][j] = static_cast<double>(i) * static_cast<double>(nrows) + static_cast<double>(j) + 0.5;
            BMat[i][j] = static_cast<double>(i) * static_cast<double>(nrows) + static_cast<double>(j) - 0.5;
        }
    }
}

// Matrix-scalar multiplication
void Multiply_Matrix_by_Scalar(Matrix<double>& A, double Scalar) {
    for (auto& row : A)
        for (double& val : row)
            val *= Scalar;
}

void Multiply_Matrix_by_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols) 
{
   unsigned int i, j;
   for (i = 0; i < nrows; i++) {
     for (j = 0; j < ncols; j++) { 
		 A[i][j] *= Scalar;
	 }
   }
}

// Matrix-scalar division
void Divide_Matrix_by_Scalar(Matrix<double>& A, double Scalar) {
    if (Scalar == 0) return;
    for (auto& row : A)
        for (double& val : row)
            val /= Scalar;
}

void Divide_Matrix_by_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols) 
{
   unsigned int i, j;
   double z = 1.0 / Scalar;
   for (i = 0; i < nrows; i++) {
     for (j = 0; j < ncols; j++) { 
		 A[i][j] *= z;
	 }
   }
}

int Zero_Matrix(Matrix<double>& A) {
    for (auto& row : A) {
        fill(row.begin(), row.end(), 0.0);
    }
    return 0;
}

int Add_Matrices(vector<double> CMat[], vector<double> AMat[], vector<double> BMat[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i, j;
   for (i = 0; i < nrows; i++) {
     for (j = 0; j < ncols; j++) { 
		 CMat[i][j] = AMat[i][j] + BMat[i][j];
	 }
   }
   return 0;
}

int Add_Matrices(Matrix<double>& CMat, const Matrix<double>& AMat, const Matrix<double>& BMat) {
    if (AMat.size() != BMat.size() || AMat[0].size() != BMat[0].size()) return -1;
    CMat.resize(AMat.size(), vector<double>(AMat[0].size()));
    for (unsigned int i=0; i<AMat.size(); i++)
        for (unsigned int j=0; j<AMat[i].size(); j++)
            CMat[i][j] = AMat[i][j] + BMat[i][j];
    return 0;
}

int Subtract_Matrices(vector<double> CMat[], vector<double> AMat[], vector<double> BMat[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i, j;
   for (i = 0; i < nrows; i++) {
     for (j = 0; j < ncols; j++) { 
		 CMat[i][j] = AMat[i][j] - BMat[i][j];
	 }
   }
   return 0;
}

void Multiply_Matrices(vector<double> CMat[], vector <double> AMat[], unsigned int nrows, vector <double> BMat[], unsigned int mcols) 
{
   unsigned int i,j,k;
   for(i = 0; i < nrows; ++i)
     for(j = 0; j < mcols; ++j)
       for(k = 0; k < nrows; ++k)
       {
         CMat[i][j] += AMat[i][k] * BMat[k][j];
       }
}

extern int Multiply_Matrix_by_Vector(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i,j;
   double sum = 0.0;
   for (i = 0; i < nrows; i++)
   {
     sum = 0.0;
     for (j = 0; j < ncols; j++)
     {
       sum += A[i][j] * X[j];
       B[i] = sum;
     }
   }
   return 0;
}

int Multiply_Vector_by_Vector(vector<double> A[], unsigned int nrows, unsigned int ncols, double MyVector1[], double MyVector2[])
{
   unsigned int i,j;
   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) A[i][j] = MyVector1[i] * MyVector2[j];
   }
   return 0;
}

extern int array_Multiply_Matrix_by_Vector(double array_A[], double X[], double B[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i,j;
   double sum = 0.0;
   for (i = 0; i < nrows; i++)
   {
     sum = 0.0;
     for (j = 0; j < ncols; j++)
     {
       sum += array_A[i*ncols+j] * X[j];
       B[i] = sum;
     }
   }
   return 0;
}

// Row transformation: row1 = row1 + Scalar*row2
void Row_Transformation(Matrix<double>& A, double Scalar, unsigned int row1, unsigned int row2) {
    if (row1 >= A.size() || row2 >= A.size()) return;
    for (unsigned int j=0; j<A[row1].size(); j++)
        A[row1][j] += Scalar * A[row2][j];
}

// Column transformation: col1 = col1 + Scalar*col2
void Column_Transformation(Matrix<double>& A, double Scalar, unsigned int col1, unsigned int col2) {
    for (unsigned int i=0; i<A.size(); i++) {
        if (col1 < A[i].size() && col2 < A[i].size())
            A[i][col1] += Scalar * A[i][col2];
    }
}

void Row_Transformation(vector<double> A[], double Scalar, unsigned int row1, unsigned int row2, unsigned int ncols)
{
  unsigned int j;
  for(j = 0; j < ncols; ++j)
	A[row2][j] = A[row2][j] + Scalar * A[row1][j];
}

void Column_Transformation(vector<double> A[], double Scalar, unsigned int col1, unsigned int col2, unsigned int nrows)
{
  unsigned int i;
  for(i = 0; i < nrows; ++i)
	A[i][col2] = A[i][col2] + Scalar * A[i][col1];
}

void Multiply_Column_by_Scalar(Matrix<double>& A, double Scalar, size_t col)
{
  unsigned int i;
  for(i = 0; i < col; ++i)
	A[i][col] *= Scalar;
}

void Multiply_Row_by_Scalar(vector<double> A[], double Scalar, unsigned int row, unsigned int ncols)
{
  unsigned int j;
  for(j = 0; j < ncols; ++j)
	A[row][j] *= Scalar;
}

int Multiply_Vector_by_Scalar(double MyVector1[], double Scalar, unsigned int nrows)
{
  unsigned int i;
  for (i=0; i < nrows; i++)
    MyVector1[i] *= Scalar;
  return 0;
}

extern void Divide_Vector_by_Scalar(double X[], double TempScalar, unsigned int ncols)
{
   double z = 1.0 / TempScalar;

   for (; ncols > 0; ncols--) *X++ *= z;
}

int Add_Vectors(double MyVector3[], double MyVector2[], double MyVector1[], unsigned int nrows) 
{
   unsigned int i;
   for (i = 0; i < nrows; i++) MyVector3[i] = MyVector2[i] + MyVector1[i];
   return 0;
}

int Subtract_Vectors(double MyVector3[], double MyVector2[], double MyVector1[], unsigned int nrows) 
{
   unsigned int i;
   for (i = 0; i < nrows; i++) MyVector3[i] = MyVector2[i] - MyVector1[i];
   return 0;
}

int Zero_Vector(double MyVector1[], unsigned int nrows)
{
  unsigned int i;
  for (i = 0; i < nrows; i++)
    MyVector1[i] = 0.0;
  return 0;
}

int Canonical_Basis_Vector(double* MyVector1, unsigned int nrows, unsigned int i) {
    unsigned int k;
    for (k = 0; k < nrows; k++) MyVector1[k] = 0.0;
    if (i < nrows) MyVector1[i] = 1.0; // Remove the i >= 0 check
    return 0;
}

int Full_Pivoting_Solution(vector<vector<double>>& A, vector<double>& B) {
    if (A.empty() || A[0].empty() || B.empty()) {
        std::cerr << "Error: Input Matrix/vector is empty.\n";
        return 1;
    }
    const size_t n = A.size(); // Use n for both rows and cols since it's square
    if (n != A[0].size() || n != B.size()) {
        std::cerr << "Error: Incompatible Matrix/vector dimensions.\n";
        return 1;
    }
    // No need to check for rectangularity again, squareness check covers it.
    vector<size_t> col_indices(n); // Keep track of column swaps
    for (size_t i = 0; i < n; ++i) {
        col_indices[i] = i;
    }
    for (size_t k = 0; k < n; ++k) {
        size_t pivot_row = k;
        size_t pivot_col = k;
        double max_val = 0.0; // Initialize to zero for correct comparison
        // Find the pivot efficiently
        for (size_t i = k; i < n; ++i) {
            for (size_t j = k; j < n; ++j) {
                double abs_val = std::fabs(A[i][j]); // Use std::fabs
                if (abs_val > max_val) {
                    max_val = abs_val;
                    pivot_row = i;
                    pivot_col = j;
                }
            }
        }
        if (max_val == 0.0) { // Check for singular Matrix
            std::cerr << "Error: Matrix is singular (no unique solution).\n";
            return 1;
        }
        // Row swap (optimized)
        if (pivot_row != k) {
            std::swap(A[k], A[pivot_row]);
            std::swap(B[k], B[pivot_row]);
        }
        // Column swap (optimized using col_indices)
        if (pivot_col != k) {
            for (size_t i = 0; i < n; ++i) {
                std::swap(A[i][k], A[i][pivot_col]);
            }
            std::swap(col_indices[k], col_indices[pivot_col]); // Update indices
        }
    }
    PrintVector(&B[0], static_cast<unsigned int>(n));
    PrintVector2D(A);
    return 0;
}

int Rook_Pivoting_Solution(Matrix<double>& A, std::vector<double>& B, size_t nrows, size_t ncols) {
    if (nrows == 0 || ncols == 0 || A.size() != nrows || B.size() != nrows) {
        return -1; // Handle invalid input
    }
    for (unsigned int k = 0; k < std::min(nrows, ncols); ++k) {
        // Find the maximum element in column k (row pivot)
        unsigned int pivot_row = k;
        double max_row_pivot = std::abs(A[k][k]); // Initialize with the diagonal element
        for (unsigned int i = k + 1; i < nrows; ++i) { // Start from k+1
            const double current = std::abs(A[i][k]);
            if (current > max_row_pivot) {
                max_row_pivot = current;
                pivot_row = i;
            }
        }
        // Find the maximum element in row k (column pivot)
        unsigned int pivot_col = k;
        double max_col_pivot = std::abs(A[k][k]); // Initialize with the diagonal element
        for (unsigned int j = k + 1; j < ncols; ++j) { // Start from k+1
            const double current = std::abs(A[k][j]);
            if (current > max_col_pivot) {
                max_col_pivot = current;
                pivot_col = j;
            }
        }
        if (max_row_pivot >= max_col_pivot) {
            // Swap rows k and pivot_row if necessary
            if (pivot_row != k) {
                std::swap_ranges(A[k].begin(), A[k].end(), A[pivot_row].begin()); // Use std::swap_ranges
                std::swap(B[k], B[pivot_row]);
            }
        } else {
            // Swap columns k and pivot_col if necessary
            if (pivot_col != k) {
                for (unsigned int i = 0; i < nrows; ++i) {
                    std::swap(A[i][k], A[i][pivot_col]);
                }
            }
        }
    }
    PrintVector(&B[0], static_cast<unsigned int>(nrows));
    PrintVector2D(A);
    return 0;
}

int Row_Pivoting_Solution(Matrix<double>& A, double B[], size_t nrows) {
    for (unsigned int k = 0; k < nrows; ++k) {
        // Find the best pivot using std::max_element and a lambda
        unsigned int p = k;
        double maxPivotRow = std::abs(A[k][k]);
        for (unsigned int i = k + 1; i < nrows; ++i) {
            if (std::abs(A[i][k]) > maxPivotRow) {
                maxPivotRow = std::abs(A[i][k]);
                p = i;
            }
        }
        // Only swap if necessary
        if (p != k) {
            std::swap(A[p], A[k]); // Swap entire rows efficiently
            std::swap(B[p], B[k]);
        }
        PrintVector(&B[0], static_cast<unsigned int>(nrows));
        PrintVector2D(A);
    }
    return 0;
}

int Column_Pivoting_Solution(Matrix<double>& A, double B[], size_t ncols) {
    if (ncols == 0 || A.empty() || A[0].size() != ncols) {
        return -1; // Handle invalid input
    }
    for (unsigned int k = 0; k < ncols; ++k) {
        // Find the best pivot in the current column (starting from row k)
        unsigned int p = k;
        double maxPivot = std::abs(A[k][k]);
        for (unsigned int i = k + 1; i < ncols; ++i) { // Start from k+1
            if (double absVal = std::abs(A[i][k]); absVal > maxPivot) {
                maxPivot = absVal;
                p = i;
            }
        }
        // Swap rows k and p (only if necessary)
        if (p != k) {
            std::swap(A[k], A[p]); // Efficient row swap
            std::swap(B[k], B[p]);
        }
    }
    PrintVector(B, static_cast<unsigned int>(ncols));
    PrintVector2D(A);
    return 0;
}

double Sum_over_Row(vector<double> A[], unsigned int ncols, unsigned int row) 
{
   double sum_row = 0.0;
   unsigned int j;
   for (j = 0; j < ncols; j++)
     sum_row += A[row][j];
   return sum_row;
}

double Sum_over_Rows(vector<double> A[], double X[], unsigned int nrows, unsigned int ncols) 
{
  unsigned int i, j;
  double sum_row;
  for (i = 0; i < nrows; i++) {
	sum_row = 0.0;
    for (j = 0; j < ncols; j++) {
      sum_row += A[i][j];
    }
    X[i] = sum_row;
  }
  return X[0];
}

double Sum_over_Column(vector<double> A[], unsigned int nrows, unsigned int col) 
{
  double sum_column = 0.0;
  unsigned int i;
  for (i = 0; i < nrows; i++)
    sum_column += A[i][col];
  return sum_column;
}

double Sum_over_Columns(vector<double> A[], double X[], unsigned int nrows, unsigned int ncols) 
{
  unsigned int i, j;
  double sum_column;
  for (j = 0; j < ncols; j++) {
	sum_column = 0.0;
    for (i = 0; i < nrows; i++) {
      sum_column += A[i][j];
    }
    X[j] = sum_column;
  }
  return X[0];
}

extern double Inner_Product(double U[], double V[], unsigned int nrows)
{
   double Inner_Product = 0.0;
   for (nrows--; nrows > 0; nrows--) Inner_Product +=  U[nrows] * V[nrows];
   return Inner_Product;
}

extern double Vector_Max_Norm(double TempVector[], unsigned int ncols)
{
   double norm = 0.0;
   double TempScalar;
   unsigned int i;

   for (i = 0; i < ncols; i++) if (norm < ( TempScalar = fabs( TempVector[i] ) ) ) norm = TempScalar;

   return norm;
}

extern double Vector_L2_Norm(double v[], unsigned int ncols)
{
   double norm = 0.0;
   unsigned int i;
   for (i = 0; i < ncols; i++) norm +=  v[i] * v[i];
   return sqrt(norm);
}

void Identity_Matrix_ut(vector<double> A[], unsigned int n) 
{
   unsigned int i, j;
   for (i = 0; i < n; ) {
      A[i][i] = 1.0;
      for (j = ++i; j < n; j++) A[i][j] = 0.0;
   }
}

void Identity_Matrix_lt(vector<double> A[], unsigned int n) 
{
   unsigned int i, j;
   for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) A[i][j] = 0.0;
      A[i][j] = 1.0;
   }
}

int Identity_Matrix(vector<double> A[], unsigned int nrows)
{
   unsigned int i,j;
   for (i = 0; i < nrows; i++) {
      for (j = 0; j < nrows; j++) {
		  if (i != j)
		    A[i][j] = 0.0;
		  else
		    A[i][j] = 1.0;
	  }
   }
   return 0; 
}

double Trace_of_Matrix(vector<double> A[], unsigned int nrows) 
{
   double trace = 0.0;
   unsigned int i;
   for (i = 0; i < nrows; i++)
     trace += A[i][i];
   return trace;
}

extern int TransposeMatrix(vector<double> A[], unsigned int nrows, unsigned int ncols)
{
  unsigned int i,j;
  double elem {0.0};
  Matrix<double> transpose;
  // Inserting elements into vector
  for (i = 0; i < nrows; i++) {
    // Vector to store column elements
    vector<double> v2;
    for (j = 0; j < ncols; j++) {
	  elem = A[i][j];
      v2.push_back(elem);
      if (j == (ncols)) {
	    v2.push_back(elem);
	  }
    }
    // to create the 2D vector
    transpose.push_back(v2);
  }
  //calculating transpose
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {	
	  transpose[j][i] = A[i][j];
    }
  }
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols; j++) {
	  A[i][j] = transpose[i][j];
    }
  }
  transpose.clear();
  transpose.shrink_to_fit();	
  return 0;
}

extern int TransposeSquareMatrix(vector<double> A[], unsigned int nrows, unsigned int ncols)
{
  unsigned int i,j;
  double elem {0.0};
  Matrix<double> transpose;
    // Inserting elements into vector
    for (i = 0; i < nrows; i++) {
      // Vector to store column elements
      vector<double> v2;
        for (j = 0; j < ncols; j++) {
		  elem = A[i][j];
          v2.push_back(elem);
            if (j == (ncols)) {
				v2.push_back(elem);
			}
        }
        // Pushing back above 1D vector
        // to create the 2D vector
        transpose.push_back(v2);
    }
    //calculating transpose
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {	
			  transpose[j][i] = A[i][j];
        }
    }
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
			A[i][j] = transpose[i][j];
        }
    }
    transpose.clear();
    transpose.shrink_to_fit();	
  return 0;
}

extern double Bilinear_Function(double u[], vector<double> A[], double v[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i,j;
   double product = 0.0;
   double sum;
   for (i = 0; i < nrows; i++)
   {
	  sum = 0.0;
      for (j = 0; j < ncols; j++)
        sum += A[i][j] * v[j];
      product += sum * u[i];
   }
   return product;
}

int Subtract_Scalar_from_Diagonal(Matrix<double>& Matrix, double scalar) {
    const size_t numRows = Matrix.size();
    const size_t numCols = (numRows > 0) ? Matrix[0].size() : 0;
    const size_t size = std::min(numRows, numCols);
    for (size_t i = 0; i < size; ++i) {
        Matrix[i][i] -= scalar;
    }
    return 0;
}

extern int Add_Scalar_to_Diagonal(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols)
{
  unsigned int i, n;
  n = ( nrows < ncols ) ? nrows : ncols;
  for (i = 0; i < n; i++)
    A[i][i] = A[i][i] + Scalar;
  return 0;
}
