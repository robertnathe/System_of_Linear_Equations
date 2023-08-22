#ifndef ArithmeticRoutines_CPP
#define ArithmeticRoutines_CPP

#include <string.h>                                 // required for memcpy()
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cfloat>
#include <math.h>		// required for fabs() and for sqrt()
using namespace std;

namespace {
	
#define a(i,j) a[(i)*nrows+(j)]

template<typename T>
using matrix = std::vector< std::vector<T> >;

/*
Iterate over vector of vectors and for each of the
nested vector print its contents
*/

void PrintVector2D(const vector<double> TempMatrix[], unsigned int nrows)
{
  std::cout << "Displaying the 2D vector:" << endl;
/*
Iterate over vector of vectors and for each of the
nested vector print its contents
*/
// Displaying the 2D vector
  std::cout << std::setprecision(5);
  for (unsigned int i = 0; i < nrows; i++) {
    for (
      auto it = TempMatrix[i].begin();
        it != TempMatrix[i].end(); it++)
        cout << *it << " ";
      cout << endl;
    }
}
void PrintVector(const double TempVector[], unsigned int nrows){
  std::cout << "Displaying the vector: " << endl;
  std::cout << std::setprecision(5);
  for (unsigned int i = 0; i < nrows; i++)
    std::cout << TempVector[i] << "   ";
  std::cout << endl;
}

int Testing_Arithmetic_Matrix_Parameters(vector<double> A[], double X[], double B[], vector <double> AMat[], vector<double> BMat[], vector<double> CMat[], double MyVector1[0], double MyVector2[0], double MyVector3[0], double Scalar, unsigned int nrows, unsigned int ncols)
{
//  Zero_Matrix(&A[0], nrows, ncols);
//  Identity_Matrix_lt(&A[0], nrows); 
//  PrintVector2D(&A[0], nrows);
//  Zero_Matrix(&A[0], nrows, ncols);
//  Identity_Matrix_ut(&A[0], nrows); 
//  PrintVector2D(&A[0], nrows);
//  Zero_Matrix(&A[0], nrows, ncols);	
//  PrintVector2D(&A[0], nrows);
//  Identity_Matrix(&A[0], nrows);
//  PrintVector2D(&A[0], nrows);	
//  Sum_over_Columns(&A[0], &X[0], nrows, ncols);
//  PrintVector(&X[0], nrows); 	
//  Sum_over_Rows(&A[0], &X[0], nrows, ncols);
//  PrintVector(&X[0], nrows); 	
//  double sum_column {0.0};
//  unsigned int col {0};
//  PrintVector2D(&A[0], nrows);
//  sum_column = Sum_over_Column(&A[0], nrows, col);
//  std::cout << "The Sum_over_Column is " << sum_column << endl;
//  double sum_row {0.0};
//  unsigned int row {2};
//  PrintVector2D(&A[0], nrows);
//  sum_row = Sum_over_Row(&A[0], ncols, row);
//  std::cout << "The Sum_over_Row is " << sum_row << endl;
//  double trace {0.0};
//  PrintVector2D(&A[0], nrows);
//  trace = Trace_of_Matrix(&A[0], nrows);
//  std::cout << "The trace of the matrix is " << trace << endl;
//  PrintVector2D(&A[0], nrows);
//  unsigned int row1 {1}, row2 {2};
//  Row_Transformation(&A[0], Scalar, row1, row2, ncols);
//  PrintVector2D(&A[0], nrows);
//  unsigned int col1 {2}, col2 {0};
//  Column_Transformation(&A[0], Scalar, col1, col2, nrows);
//  PrintVector2D(&A[0], nrows);
//  unsigned int row1 {0}, row2 {1};
//  PrintVector2D(&A[0], nrows);
//  Row_Transformation(&A[0], Scalar, row1, row2, ncols);
//  PrintVector2D(&A[0], nrows);   
//  unsigned int col1 {0}, col2 {1};	
//  PrintVector2D(&A[0], nrows);
//  Column_Transformation(&A[0], Scalar, col1, col2, nrows, ncols);
//  PrintVector2D(&A[0], nrows);
//  PrintVector2D(&A[0], nrows);
//  unsigned int row = 0;
//  Multiply_Row_by_Scalar(&A[0], Scalar, row, ncols);
//  PrintVector2D(&A[0], nrows);
//  PrintVector2D(&A[0], nrows);
//  unsigned int col = 2;
//  Multiply_Column_by_Scalar(&A[0], Scalar, col, nrows);
//  PrintVector2D(&A[0], nrows); 
//  unsigned int i,j;
//  unsigned int mcols = ncols;
//  for (i = 0; i < nrows; i++) {
//    for (j = 0; j < ncols; j++) { 
//	  CMat[i][j] = 0.0;
//	  AMat[i][j] = i*nrows+j;
//	  BMat[i][j] = i*nrows+j;
//    }	 
//  }
//  PrintVector2D(&AMat[0], nrows);
//  PrintVector2D(&BMat[0], nrows);
//  Multiply_Matrices(&CMat[0], &AMat[0], nrows, ncols, &BMat[0], mcols);
//  PrintVector2D(&CMat[0], nrows);
//  PrintVector2D(&A[0], nrows);
//  Multiply_Matrix_by_Scalar(&A[0], Scalar, nrows, ncols);
//  PrintVector2D(&A[0], nrows);
//  Divide_Matrix_by_Scalar(&A[0], Scalar, nrows, ncols);
//  PrintVector2D(&A[0], nrows);
//  unsigned int i, j;
//  for (i = 0; i < nrows; i++) {
//    for (j = 0; j < ncols; j++) { 
//	  AMat[i][j] = i*nrows+j+0.5;
//	  BMat[i][j] = i*nrows+j-0.5;
//	}
//  }
//  PrintVector2D(&AMat[0], nrows);
//  PrintVector2D(&BMat[0], nrows);
//  Add_Matrices(&CMat[0], &AMat[0], &BMat[0], nrows, ncols);	
//  PrintVector2D(&CMat[0], nrows);
//  PrintVector2D(&AMat[0], nrows);
//  PrintVector2D(&BMat[0], nrows);
//  Subtract_Matrices(&CMat[0], &AMat[0], &BMat[0], nrows, ncols);	
//  PrintVector2D(&CMat[0], nrows);
//	Multiply_Vector_by_Vector(&A[0], nrows, ncols, &MyVector1[0], &MyVector2[0]);
//	PrintVector2D(&A[0], nrows);
//	Zero_Vector(&MyVector1[0], nrows);
//	PrintVector(&MyVector1[0], nrows);
//	unsigned int i {1};
//	Canonical_Basis_Vector(&MyVector1[0], i, nrows);
//	PrintVector(&MyVector1[0], nrows);
//	Add_Vectors(&MyVector3[0], &MyVector2[0], &MyVector1[0], nrows);
//	PrintVector(&MyVector3[0], nrows);
//	Subtract_Vectors(&MyVector3[0], &MyVector2[0], &MyVector1[0], nrows);
//    PrintVector(&MyVector3[0], nrows);
//	Multiply_Vector_by_Scalar(&MyVector1[0], Scalar, nrows);
//	PrintVector(&MyVector1[0], nrows);
//	std::cout << "The vector X is the following: " << endl;
//	PrintVector(&X[0], nrows);
//	std::cout << "The vector B is the following: " << endl;
//	PrintVector(&B[0], nrows);
//	product = Bilinear_Function(&X[0], &A[0], &B[0], nrows, ncols);
//	std::cout << "The product X'AB is the following: " << endl;
//	std::cout << setprecision(7) << product << endl;
//    PrintVector2D(&A[0], nrows);
//    std::cout << "The scalar is " << Scalar << endl;
//     Subtract scalar from diagonal
//    Subtract_Scalar_from_Diagonal(&A[0], Scalar, nrows, ncols);
//    std::cout << "The matrix A is " << endl;
//    PrintVector2D(&A[0], nrows);
//	std::cout << setprecision(7) << product << endl;
//    PrintVector2D(&A[0], nrows);
//    std::cout << "The scalar is " << Scalar << endl;
//    Add scalar to diagonal(&A[0], Scalar, nrows, ncols);
//    std::cout << "The matrix A is " << endl;
//    PrintVector2D(&A[0], nrows);
//    Copy vector
//    vector <double> vd(X);
//    cout << "Old vector elements are : ";
//    for (i=0; i<X.size(); i++)
//        cout << X[i] << " ";
//    cout << endl;
//    cout << "New vector elements are : ";
//    for (i=0; i<vd.size(); i++)
//        cout << vd[i] << " ";
//    cout<< endl;
//	std::cout << "The computed inner product is the following: " << endl;
//	inner_product = std::inner_product(X.begin(), X.end(), B.begin(), 0);
//	std::cout << setprecision(7) << inner_product << endl;
//    std::cout << "The computed norm is the following: " << endl;
//    norm = Vector_Max_Norm(&B[0], ncols);
//    std::cout << setprecision(7) << norm << endl;
//	double TempScalar {1.0};
//    if ( TempScalar != 0.0)  Divide_Vector_by_Scalar(&X[0], TempScalar,ncols);
//    std:: cout <<   "The vector X is : " << endl;
//    for (unsigned int i = 0; i < ncols; i++)
//      std::cout << X[i] << "   ";
//    std::cout << endl;
//    std::cout << "The computed norm is the following: " << endl;
//    norm = Vector_L2_Norm(&B[0], ncols);
//    std::cout << setprecision(7) << norm << endl;
//    Multiply_Matrix_by_Vector(&A[0], &X[0], &B[0], nrows, ncols);
//    Output_Data (&A[0], &X[0], &B[0], nrows, ncols);
//     Input &X[0] insted of input &X0[0] This is for inverse power method Parameter_Solution.
//    Parameter_Solution (&A[0], &X[0], &B[0], nrows, ncols, eigenvalue, tolerance, max_iter);
//    Output_Data (&A[0], &X[0], &B[0], nrows, ncols, eigenvalue, tolerance, max_iter);
//    std::cout << "The vector P before LUPDecompose is the following: " << endl;
//	for (i = 0; i < ncols; i++){
//	  std::cout << P[i] << endl;
//	}
//	std::cout << endl;
//    PrintVector2D(&A[0], nrows);
//    RowPivot(&A[0], &X[0], &B[0], &P[0], ncols, tolerance);
//    PrintVector2D(&A[0], nrows);
//    Output_Data (&A[0], &X[0], &B[0], nrows, ncols, eigenvalue, tolerance, max_iter);
//    std::cout << "The vector P after LUPDecompose is the following: " << endl;
//	for (i = 0; i < ncols; i++){
//	  std::cout << P[i] << endl;
//	}
//	std::cout << endl;
//    PrintVector(&X[0], nrows);
//    std::cout << "The X vector should be the following: -0.436 0.430 5.12" << endl;

//	Testing_Function (&A[0], &X[0], &B[0], nrows, ncols);
//    Output_Results (&X[0], &B[0], nrows, ncols, tolerance, max_iter);

	std::cout << "The matrix A is the following: " << endl;
	PrintVector2D(&A[0], nrows);
    return 0;
}

class Input {
  public:
  unsigned int Input_Begin(void)
  {
    unsigned int i,j, max_iter {1000};
    (void) max_iter;
    double tolerance {0.00001};
    (void) tolerance;
    double eigenvalue {0.0};
    (void) eigenvalue;
    double Scalar {2.5};
    (void) Scalar;
    FILE *myFile;
    myFile = fopen ("A.dat", "r");
    unsigned int dimArray[3];
    unsigned int nrows, number_of_entries_A;
    unsigned int ncols;
    unsigned int i_index, j_index;
    double value, elem {0.0};
    while (myFile == NULL)
    {
	  std::cout << "Error Reading File" << endl;
      exit (0);
    }
    fscanf (myFile, "%*s %*s %*s %*s %*s");
    for (i = 0; i < 3; i++)
    {
    fscanf (myFile, "%u,", &dimArray[i]);
    }
    nrows = dimArray[0];
    ncols = dimArray[1];
    number_of_entries_A = dimArray[2];
    vector <double> array_A;
    for (i = 0; i < number_of_entries_A; i++)
    {
      fscanf (myFile, "%u,", &i_index);
      i_index--;
      fscanf (myFile, "%u,", &j_index);
      j_index--;
      fscanf (myFile, "%lf,", &value);
//    Change program to use the single index array_A
      array_A.push_back(value);
    }
    fclose (myFile);
    FILE *myFile2;
    myFile2 = fopen ("B.dat", "r");
    unsigned int dim_B_Array[3];
    unsigned int number_of_entries_B;
    unsigned int col_B;
    (void) col_B;
    while (myFile2 == NULL)
    {
	  std::cout << "Error Reading File" << endl;
      exit (0);
    }
    fscanf (myFile2, "%*s %*s %*s %*s %*s");
    for (i = 0; i < 3; i++)
    {
      fscanf (myFile2, "%u,", &dim_B_Array[i]);
    }
    col_B = dim_B_Array[1];
    number_of_entries_B = dim_B_Array[2];
    vector <double> B;
    vector <double> X;
    vector <double> U;
    vector <double> V;
    vector <double> W;
    vector <int> P;
    for (i = 0; i < number_of_entries_B; i++)
    {
      fscanf (myFile2, "%u,", &i_index);
      i_index--;
      fscanf (myFile2, "%u,", &j_index);
      j_index--;
      fscanf (myFile2, "%lf,", &value);
      B.push_back(value);
      X.push_back(0.0);
      U.push_back(0.0);
      V.push_back(0.0);
      W.push_back(0.0);
      P.push_back(0);
    }
    fclose (myFile2);
    // Initializing the vector of vectors
//    vector<vector<double> > A;
    matrix<double> A, AMat, BMat, CMat;
    // Inserting elements into vector
    for (i = 0; i < nrows; i++) {
      // Vector to store column elements
      vector<double> Row1, Row2, Row3, Row4;
        for (j = 0; j < ncols; j++) {
		  elem = array_A[i*nrows+j];
          Row1.push_back(elem);
          Row2.push_back(elem);
          Row3.push_back(elem);
          Row4.push_back(elem);
            if (j == (ncols)) {
				Row1.push_back(elem);
				Row2.push_back(elem);
				Row3.push_back(elem);
				Row4.push_back(elem);				
			}
        }
        // Pushing back above 1D vector
        // to create the 2D vector
        A.push_back(Row1);
        AMat.push_back(Row2);
        BMat.push_back(Row3);
        CMat.push_back(Row4);
    }
    matrix<double> Augmented;
    // Elements to insert in column
    // Inserting elements into vector
    for (i = 0; i < nrows; i++) {
    // Vector to store column elements
      vector<double> Augment;
      for (j = 0; j < ncols+1; j++) {
        Augment.push_back(elem);
      }
      // Pushing back above 1D vector
      // to create the 2D vector
      Augmented.push_back(Augment);
    }
    for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < ncols+1; j++)
	  {
	    if (j != (ncols)) {
	      Augmented[i][j] = A[i][j];
	    }
	    if (j == (ncols)) {
	      Augmented[i][ncols] = B[i];
	    }
      }
    }
    
//    Full_Pivoting_Solution(&A[0], &X[0], &B[0], nrows, ncols);
    
//    Testing_General_Purpose_Parameters (&A[0], &X[0], &B[0], nrows, ncols);
//    Testing_Arithmetic_Matrix_Parameters (&A[0], &X[0], &B[0], nrows, ncols); 
//    Rook_Pivoting_Solution(&A[0], &X[0], &B[0], nrows, ncols);
//    Row_Pivoting_Solution(&A[0], &X[0], &B[0], nrows, ncols);
//    Column_Pivoting_Solution(&A[0], &X[0], &B[0], nrows, ncols);
    array_A.clear();
    array_A.shrink_to_fit();
    B.clear();
    B.shrink_to_fit();
    X.clear();
    X.shrink_to_fit();
    A.clear();
    A.shrink_to_fit();
    AMat.clear();
    AMat.shrink_to_fit();
    BMat.clear();
    BMat.shrink_to_fit();
    CMat.clear();
    CMat.shrink_to_fit();
    Augmented.clear();
    Augmented.shrink_to_fit();
    U.clear();
    U.shrink_to_fit();
    V.clear();
    V.shrink_to_fit();
    W.clear();
    W.shrink_to_fit();
    P.clear();
    P.shrink_to_fit();
//    std::cin.clear(); // reset any error flags
//    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//    ignore any characters in the input buffer until we find an enter character
//    std::cin.get(); // get one more char from the user
    return 0;
  }
};

} // close namespace

////////////////////////////////////////////////////////////////////////////////
// File: add_matrices.c                                                       //
// Routine(s):                                                                //
//    Add_Matrices                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Add_Matrices(double *C, double *A, double *B, int nrows, int ncols)  //
//                                                                            //
//  Description:                                                              //
//     This routine computes C = A + B where C is an nrows x ncols real matrix//
//     and A and B are given nrows x ncols real matrices.                     //
//                                                                            //
//     All matrices should be declared as " double X[nrows][ncols] " in the   //
//     calling routine, where X = A, B, C.                                    //
//                                                                            //
//  Arguments:                                                                //
//     double *C    Address of the first element of the matrix C.             //
//     double *A    Address of the first element of the matrix A.             //
//     double *B    Address of the first element of the matrix B.             //
//     int    nrows The number of rows of matrices A, B, and C.               //
//     int    ncols The number of columns of the matrices A, B, and C.        //
//                                                                            //
//  Return Values:                                                            //
//      void                                                                  //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  B[M][N], C[M][N];                                     //
//                                                                            //
//     (your code to initialize the matrices A and B)                         //
//                                                                            //
//     Add_Matrices((double *) C, &A[0][0], &B[0][0], M, N);                  //
//     printf("The matrix C = A + B is \n"); ...                              //
////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
// File: subtract_matrices.c                                                  //
// Routine(s):                                                                //
//    Subtract_Matrices                                                       //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Subtract_Matrices(double *C, double *A, double *B, int nrows,        //
//                                                               int ncols)   //
//                                                                            //
//  Description:                                                              //
//     This routine computes C = A - B where C is an nrows x ncols real matrix//
//     and A and B are given nrows x ncols real matrices.                     //
//                                                                            //
//     All matrices should be declared as " double X[nrows][ncols] " in the   //
//     calling routine, where X = A, B, C.                                    //
//                                                                            //
//  Arguments:                                                                //
//     double *C    Pointer to the first element of the matrix C.             //
//     double *A    Pointer to the first element of the matrix A.             //
//     double *B    Pointer to the first element of the matrix B.             //
//     int    nrows The number of rows of matrices A and B.                   //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  B[M][N], C[M][N];                                     //
//                                                                            //
//     (your code to initialize the matrices A and B)                         //
//                                                                            //
//     Subtract_Matrices(&C[0][0], &A[0][0], &B[0][0], M, N);                 //
//     printf("The matrix C = A - B  is \n"); ...                             //
////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
// File: sub_scalar_from_diagonal.c                                           //
// Routine(s):                                                                //
//    Subtract_Scalar_from_Diagonal                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Subtract_Scalar_from_Diagonal(double *A, double x, int nrows,        //
//                                                               int ncols)   //
//                                                                            //
//  Description:                                                              //
//     Subtract the scalar x from the each element of the diagonal of the     //
//     matrix A.                                                              //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     double x     Scalar to be subtracted from each diagonal element of the //
//                  matrix A.                                                 //
//     int    nrows The number of rows of matrix A.                           //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  x;                                                    //
//                                                                            //
//     (your code to initialize the matrix A and scalar x)                    //
//                                                                            //
//     Subtract_Scalar_from_Diagonal(&A[0][0], x, M, N);                      //
//     printf("The matrix A is \n"); ... }                                    //
////////////////////////////////////////////////////////////////////////////////
extern int Subtract_Scalar_from_Diagonal(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols)
{
   unsigned int i, n;
   n = ( nrows < ncols ) ? nrows : ncols;
   for (i = 0; i < n; i++)
     A[i][i] = A[i][i] - Scalar;
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

////////////////////////////////////////////////////////////////////////////////
// File: mul_matrix_by_scalar.c                                               //
// Routine(s):                                                                //
//    Multiply_Matrix_by_Scalar                                               //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Multiply_Matrix_by_Scalar(double *A, double x, int nrows, int ncols) //
//                                                                            //
//  Description:                                                              //
//     Multiply each element of the matrix A by the scalar x.                 //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     double x     Scalar to multipy each element of the matrix A.           //
//     int    nrows The number of rows of matrix A.                           //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  x;                                                    //
//                                                                            //
//     (your code to initialize the matrix A and scalar x)                    //
//                                                                            //
//     Multiply_Matrix_by_Scalar(&A[0][0], x, M, N);                          //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Multiply_Matrix_by_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols) 
{
   unsigned int i, j;
   for (i = 0; i < nrows; i++) {
     for (j = 0; j < ncols; j++) { 
		 A[i][j] *= Scalar;
	 }
   }
}
////////////////////////////////////////////////////////////////////////////////
// File: div_matrix_by_scalar.c                                               //
// Routine(s):                                                                //
//    Divide_Matrix_by_Scalar                                                 //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Divide_Matrix_by_Scalar(double *A, double x, int nrows, int ncols)   //
//                                                                            //
//  Description:                                                              //
//     Divide each element of the matrix A by the scalar x.                   //
//     i.e.       A[i][j] <- A[i][j] / x for all i,j.                         //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     double x     The non-zero scalar used to divide each element of the    //
//                  matrix A.                                                 //
//     int    nrows The number of rows of the matrix A.                       //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  x;                                                    //
//                                                                            //
//     (your code to initialize the matrix A and scalar x)                    //
//                                                                            //
//     if ( x != 0.0)  Divide_Matrix_by_Scalar(&A[0][0], x, M, N);            //
//      printf("The matrix A is \n"); ...                                     //
////////////////////////////////////////////////////////////////////////////////
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
void Multiply_Matrices(vector<double> CMat[], vector <double> AMat[], unsigned int nrows, unsigned int ncols, vector <double> BMat[], unsigned int mcols) 
{
   unsigned int i,j,k;
   for(i = 0; i < nrows; ++i)
     for(j = 0; j < mcols; ++j)
       for(k = 0; k < nrows; ++k)
       {
         CMat[i][j] += AMat[i][k] * BMat[k][j];
       }
}
////////////////////////////////////////////////////////////////////////////////
// File: multiply_matrix_by_vector.c                                          //
// Routine(s):                                                                //
//    Multiply_Matrix_by_Vector                                               //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Multiply_Matrix_by_Vector(double u[], double *A, int nrows,          //
//                                                   int ncols, double v[])   //
//                                                                            //
//  Description:                                                              //
//     Post multiply the nrows x ncols matrix A by the column vector v        //
//     to form the column vector u, i.e. u = A v.                             //
//     The matrix A should be declared as "double A[nrows][ncols]" in the     //
//     calling routine.  The vector v declared as "double v[ncols]" and       //
//     the vector u declared as "double u[nrows]" in the calling routine.     //
//                                                                            //
//  Arguments:                                                                //
//     double *u    Pointer to the first element of the vector u.             //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of the matrix A and the number of      //
//                  components of the column vector u.                        //
//     int    ncols The number of columns of the matrices A and the           //
//                  number of components of the column vector v.              //
//     double *v    Pointer to the first element of the vector v.             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  u[M], v[N];                                           //
//                                                                            //
//     (your code to initialize the matrix A and column vector v)             //
//                                                                            //
//     Multiply_Matrix_by_Vector(u, &A[0][0], M, N, v);                       //
//     printf("The vector u is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
// File: row_transformation.c                                                 //
// Routine(s):                                                                //
//    Row_Transformation                                                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Row_Transformation(double *A, double x, int row1, int row2,          //
//                                                                 int ncols) //
//                                                                            //
//  Description:                                                              //
//     Multiply the row 'row1' by x and add to 'row2' of the  nrows x ncols   //
//     matrix A, i.e. A[row2][j] <- A[row2][j] + x * A[row1][j], for all j.   //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     double x     Scalar used to multiply each element of row1 of A.        //
//     int    row1  The row of A which is multiplied by x and then added to   //
//                  row2 of A.                                                //
//     int    row2  The row of A which is changed to the sum of its old value //
//                  and x times the corresponding element of row1.            //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N], x;                                                     //
//     int i, j;                                                              //
//                                                                            //
//     (your code to create the matrix A, the scalar x, the row number i and  //
//      row number j)                                                         //
//                                                                            //
//     if ( ( i < M ) && ( j < M) ) {                                         //
//        Row_Transformation(&A[0][0], i, j, N);                              //
//         printf("The matrix A is \n"); ...                                  //
//     } else printf("Illegal row numbers.\n");                               //
////////////////////////////////////////////////////////////////////////////////
void Row_Transformation(vector<double> A[], double Scalar, unsigned int row1, unsigned int row2, unsigned int ncols)
{
  unsigned int j;
  for(j = 0; j < ncols; ++j)
	A[row2][j] = A[row2][j] + Scalar * A[row1][j];
}
////////////////////////////////////////////////////////////////////////////////
// File: col_transformation.c                                                 //
// Routine(s):                                                                //
//    Column_Transformation                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Column_Transformation(double *A, double x, int col1, int col2,       //
//                                                      int nrows, int ncols) //
//                                                                            //
//  Description:                                                              //
//     Multiply the column 'col1' by x and add to 'col2' of the               //
//     nrows x ncols matrix A, i.e.                                           //
//           A[i][col2] <- A[i][col2] + x * A[i][col1], for all i.            //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     double x     Scalar used to multiply each element of column col1 of A. //
//     int    col1  The column of A which is multiplied by x and then added   //
//                  to col2 of A.                                             //
//     int    col2  The column of A which is changed to the sum of its old    //
//                  value and x times the corresponding element of col1.      //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N], x;                                                     //
//     int i, j;                                                              //
//                                                                            //
//     (your code to create the matrix A, the scalar x, the column number i   //
//      and column number j)                                                  //
//                                                                            //
//     if ( ( i < N ) && ( j < N) ) {                                         //
//        Column_Transformation(&A[0][0], x, i, j, M, N);                     //
//         printf("The matrix A is \n"); ...                                  //
//     } else printf("Illegal column numbers.\n");                            //
////////////////////////////////////////////////////////////////////////////////
void Column_Transformation(vector<double> A[], double Scalar, unsigned int col1, unsigned int col2, unsigned int nrows)
{
  unsigned int i;
  for(i = 0; i < nrows; ++i)
	A[i][col2] = A[i][col2] + Scalar * A[i][col1];
}
////////////////////////////////////////////////////////////////////////////////
// File: mul_matrix_row_by_scalar.c                                           //
// Routine(s):                                                                //
//    Multiply_Row_by_Scalar                                                  //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Multiply_Row_by_Scalar(double *A, double x, int row, int ncols)      //
//                                                                            //
//  Description:                                                              //
//     Multiply the row 'row' by x of the  nrows x ncols matrix A, i.e.       //
//                  A[row][j] <- x * A[row][j], for all j.                    //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     double x     Scalar used to multiply each element of row 'row' of A.   //
//     int    row   The row of A which is multiplied by x.                    //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N], x;                                                     //
//     int i;                                                                 //
//                                                                            //
//     (your code to create the matrix A, the scalar x, the row number i )    //
//                                                                            //
//     if ( i < M  ) {                                                        //
//        Multiply_Row_by_Scalar(&A[0][0], x, i, N);                          //
//         printf("The matrix A is \n"); ...                                  //
//     } else printf("Illegal row number.\n");                                //
////////////////////////////////////////////////////////////////////////////////
void Multiply_Row_by_Scalar(vector<double> A[], double Scalar, unsigned int row, unsigned int ncols)
{
  unsigned int j;
  for(j = 0; j < ncols; ++j)
	A[row][j] *= Scalar;
}
////////////////////////////////////////////////////////////////////////////////
// File: mul_matrix_col_by_scalar.c                                           //
// Routine(s):                                                                //
//    Multiply_Column_by_Scalar                                               //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Multiply_Column_by_Scalar(double *A, double x, int col, int nrows,   //
//                                                                 int ncols) //
//                                                                            //
//  Description:                                                              //
//     Multiply the column 'col' by x of the nrows x ncols matrix A, i.e.     //
//               A[i][col] <- x * A[i][col], for all i.                       //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     double x     Scalar used to multiply each element of column col of A.  //
//     int    col   The column of A which is multiplied by x.                 //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N], x;                                                     //
//     int i;                                                                 //
//                                                                            //
//     (your code to create the matrix A, the scalar x, the column number i)  //
//                                                                            //
//     if ( i < N )  {                                                        //
//        Multiply_Column_by_Scalar(&A[0][0], x, i, M, N);                    //
//         printf("The matrix A is \n"); ...                                  //
//     } else printf("Illegal column number.\n");                             //
////////////////////////////////////////////////////////////////////////////////
void Multiply_Column_by_Scalar(vector<double> A[], double Scalar, unsigned int col, unsigned int nrows) 
{
  unsigned int i;
  for(i = 0; i < nrows; ++i)
	A[i][col] *= Scalar;
}

////////////////////////////////////////////////////////////////////////////////
// File: sum_over_row.c                                                       //
// Routine(s):                                                                //
//    Sum_over_Row                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Sum_over_Row(double *A, int ncols, int row)                        //
//                                                                            //
//  Description:                                                              //
//     Sums over the row 'row' of the nrows x ncols matrix A.                 //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    ncols The number of columns of the matrix A.                    //
//     int    row   The row number which is summed.                           //
//                                                                            //
//  Return Values:                                                            //
//     double sum   The sum of the row 'row' of the matrix A.                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  sum;                                                  //
//                                                                            //
//     (your code to initialize the matrix A and the row number row)          //
//                                                                            //
//     if ( row < M )                                                         //
//        sum = Sum_over_Row(&A[0][0], N, row);                               //
//     printf("The row sum is \n"); ...                                       //
////////////////////////////////////////////////////////////////////////////////
double Sum_over_Row(vector<double> A[], unsigned int ncols, unsigned int row) 
{
   double sum_row = 0.0;
   unsigned int j;
   for (j = 0; j < ncols; j++)
     sum_row += A[row][j];
   return sum_row;
}
////////////////////////////////////////////////////////////////////////////////
// File: sum_over_rows.c                                                      //
// Routine(s):                                                                //
//    Sum_over_Rows                                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Sum_over_Rows(double v[], double *A, int nrows, int ncols)           //
//                                                                            //
//  Description:                                                              //
//     The ith component of the vector v[] is set to the sum of the ith row   //
//     of the matrix A.  This is the same as postmultiplying A by the vector  //
//     all of whose components are 1.                                         //
//                                                                            //
//  Arguments:                                                                //
//     double v[]   Pointer to the first element of the vector v the ith      //
//                  component of which is the sum over the ith row of A.      //
//     double *A    Pointer to the first element of the matrix A.             //
//                  matrix A.                                                 //
//     int    nrows The number of rows of matrix A or dimension of v[].       //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[M];                                                 //
//                                                                            //
//     (your code to initialize the matrix A)                                 //
//                                                                            //
//     Sum_over_Rows(v, &A[0][0], M, N);                                      //
//     printf("The row sums v are \n"); ... }                                 //
////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
// File: sum_over_col.c                                                       //
// Routine(s):                                                                //
//    Sum_over_Column                                                         //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Sum_over_Column(double *A, int nrows, int ncols, int col)          //
//                                                                            //
//  Description:                                                              //
//     Sum over column 'col' of the nrows x ncols matrix A.                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of matrix A.                           //
//     int    ncols The number of columns of the matrix A.                    //
//     int    col   The column which is summed.                               //
//                  Note that the sum is over A[*][col]                       //
//                                                                            //
//  Return Values:                                                            //
//     double sum   The sum of the column 'col'.                              //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  sum;                                                  //
//     int col;                                                               //
//                                                                            //
//     (your code to initialize the matrix A)                                 //
//                                                                            //
//     if ( (col > 0) && (col < N) )                                          //
//        sum = Sum_over_Column(&A[0][0], M, N, col);                         //
//     printf("The column sum is \n"); ...                                    //
////////////////////////////////////////////////////////////////////////////////
double Sum_over_Column(vector<double> A[], unsigned int nrows, unsigned int col) 
{
  double sum_column = 0.0;
  unsigned int i;
  for (i = 0; i < nrows; i++)
    sum_column += A[i][col];
  return sum_column;
}
////////////////////////////////////////////////////////////////////////////////
// File: sum_over_cols.c                                                      //
// Routine(s):                                                                //
//    Sum_over_Columns                                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Sum_over_Columns(double v[], double *A, int nrows, int ncols)        //
//                                                                            //
//  Description:                                                              //
//     The ith component of the vector v[] is set to the sum of the ith       //
//     column of the matrix A.  This is the same as premultiplying A by the   //
//     vector all of whose components are 1.                                  //
//                                                                            //
//  Arguments:                                                                //
//     double v[]   Pointer to the first element of the vector v the ith      //
//                  component of which is the sum over the ith column of A.   //
//     double *A    Pointer to the first element of the matrix A.             //
//                  matrix A.                                                 //
//     int    nrows The number of rows of matrix A.                           //
//     int    ncols The number of columns of the matrix A or the dimension of //
//                  the vector v[].                                           //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[N];                                                 //
//                                                                            //
//     (your code to initialize the matrix A)                                 //
//                                                                            //
//     Sum_over_Columns(v, &A[0][0], M, N);                                   //
//     printf("The column sums are v \n"); ...                                //
////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
// File: trace_of_matrix.c                                                    //
// Routine(s):                                                                //
//    Trace_of_Matrix                                                         //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Trace_of_Matrix(double *A, int n)                                  //
//                                                                            //
//  Description:                                                              //
//     The trace of a square matrix is the sum of the diagonal elements of    //
//     that matrix.  This routines calculates the trace of the square n x n   //
//     matrix A.                                                              //
//                                                                            //
//  Arguments:                                                                //
//     double *A      Pointer to the first element of the matrix A.           //
//     int    n       The number of rows and columns of the square matrix A.  //
//                                                                            //
//  Return Values:                                                            //
//     double trace;                                                          //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], trace;                                                 //
//                                                                            //
//     (your code to initialize the matrix A)                                 //
//                                                                            //
//     trace = Trace_of_Matrix(&A[0][0], N);                                  //
//     printf("The trace of A is \n"); ...                                    //
////////////////////////////////////////////////////////////////////////////////
double Trace_of_Matrix(vector<double> A[], unsigned int nrows) 
{
   double trace = 0.0;
   unsigned int i;
   for (i = 0; i < nrows; i++)
     trace += A[i][i];
   return trace;
}
////////////////////////////////////////////////////////////////////////////////
// File: zero_matrix.c                                                        //
// Routine(s):                                                                //
//    Zero_Matrix                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Zero_Matrix(double *A, int nrows, int ncols)                         //
//                                                                            //
//  Description:                                                              //
//     Set the nrows x ncols matrix A equal to the zero matrix, i.e.          //
//     A[i][j] = 0 for all i, j.                                              //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of the matrix A.                       //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[N][M];                                                        //
//                                                                            //
//     Zero_Matrix(&A[0][0], N, M);                                           //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
int Zero_Matrix(vector<double> A[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i,j;
   for (i = 0; i < nrows; i++) {
	   for (j = 0; j < ncols; j++)
     A[i][j] = 0.0;
   }
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
// File: identity_matrix.c                                                    //
// Routine(s):                                                                //
//    Identity_Matrix                                                         //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Identity_Matrix(double *A, int n)                                    //
//                                                                            //
//  Description:                                                              //
//     Set the square n x n matrix A equal to the identity matrix, i.e.       //
//     A[i][j] = 0 if i != j and A[i][i] = 1.                                 //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    n     The number of rows and columns of the matrix A.           //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N];                                                        //
//                                                                            //
//     Identity_Matrix(&A[0][0], N);                                          //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
// File: identity_matrix_ut.c                                                 //
// Routines:                                                                  //
//    Identity_Matrix_ut                                                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Identity_Matrix_ut( double *A, int n)                                //
//                                                                            //
//  Description:                                                              //
//     Set the square symmetric matrix A stored in upper triangular form      //
//     equal to the identity matrix.                                          //
//                                                                            //
//  Arguments:                                                                //
//     double *A                                                              //
//        On output A is set to the identity matrix where A is stored in      //
//        upper triangular form.  I.e  A[0] = 1, A[1] = 0, ,,,, A[n-1] = 0,   //
//        A[n] = 1, A[n+1] = 0,..., A[2n-2] = 0,..., A[n(n+1)/2] = 1.         //
//     int    n                                                               //
//        The dimension A, i.e. A is an n x n symmetric matrix.  It should    //
//        be declared as an array dimensioned at least a large as             //
//        n * (n + 1) / 2 in the calling routine.                             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[(N * (N + 1)) >>1];                                           //
//                                                                            //
//     Identity_Matrix_ut(A, n);                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void Identity_Matrix_ut(vector<double> A[], unsigned int n) 
{
   unsigned int i, j;
   for (i = 0; i < n; ) {
      A[i][i] = 1.0;
      for (j = ++i; j < n; j++) A[i][j] = 0.0;
   }
}
////////////////////////////////////////////////////////////////////////////////
// File: identity_matrix_lt.c                                                 //
// Routines:                                                                  //
//    Identity_Matrix_lt                                                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Identity_Matrix_lt( double *A, int n)                                //
//                                                                            //
//  Description:                                                              //
//     Set the square symmetric matrix A stored in lower triangular form      //
//     equal to the identity matrix.                                          //
//                                                                            //
//  Arguments:                                                                //
//     double *A                                                              //
//        On output A is set to the identity matrix where A is stored in      //
//        lower triangular form.  I.e  A[0] = 1, A[1] = 0, A[2] = 1,          //
//        A[3] = 0, A[4] = 0, A[5] = 1, ...,  A[n(n+1)/2] = 1.                //
//     int    n                                                               //
//        The dimension A, i.e. A is an n x n symmetric matrix.  It should    //
//        be declared as an array dimensioned at least a large as             //
//        n * (n + 1) / 2 in the calling routine.                             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[(N * (N + 1)) >>1];                                           //
//                                                                            //
//     Identity_Matrix_lt(A, n);                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
void Identity_Matrix_lt(vector<double> A[], unsigned int n) 
{
   unsigned int i, j;
   for (i = 0; i < n; i++) {
      for (j = 0; j < i; j++) A[i][j] = 0.0;
      A[i][j] = 1.0;
   }
}

////////////////////////////////////////////////////////////////////////////////
// File: bilinear_function.c                                                  //
// Routine(s):                                                                //
//    Bilinear_Function                                                       //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Bilinear_Function(double u[], double *A, double v[],               //
//                                                     int nrows, int ncols)  //
//                                                                            //
//  Description:                                                              //
//     Pre-multiply the nrows x ncols matrix A by the row vector u and post-  //
//     multiply by the column vector v, to form the scalar u'Av.              //
//     The matrix A should be declared as "double A[nrows][ncols]" in the     //
//     calling routine.  The vector u declared as "double u[nrows]" and       //
//     the vector v declared as "double v[ncols]" in the calling routine.     //
//                                                                            //
//  Arguments:                                                                //
//     double *u    Pointer to the first element of the vector u.             //
//     double *A    Pointer to the first element of the matrix A.             //
//     double *v    Pointer to the first element of the vector v.             //
//     int    nrows The number of rows of the matrix A and the number of      //
//                  components of the row vector u.                           //
//     int    ncols The number of columns of the matrices A and the           //
//                  number of components of the column vector v.              //
//                                                                            //
//  Return Values:                                                            //
//     Product of u' A v.                                                     //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  u[M], v[N];                                           //
//     double product;                                                        //
//                                                                            //
//     (your code to initialize the matrix A, row vector u and                //
//                                                         column vector v)   //
//                                                                            //
//     product = Bilinear_Function(u, &A[0][0], v, M, N);                     //
//     printf("The product u'Av is \n"); ...                                  //
////////////////////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////////////////////
// File: add_vectors.c                                                        //
// Routine(s):                                                                //
//    Add_Vectors                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Add_Vectors(double *w, double *u, double *v, int n)                  //
//                                                                            //
//  Description:                                                              //
//     Add vectors u and v to form the vector w, i.e. w = u + v, where        //
//     w[i] = u[i] + v[i].  All vectors u,v,w should be declared as           //
//     " double [n] " in the calling routine.                                 //
//                                                                            //
//  Arguments:                                                                //
//     double w[]   Resultant vector w = u + v.                               //
//     double u[]   A summand.                                                //
//     double v[]   The other summand.                                        //
//     int    n     The number of components of u, v, and w.                  //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double u[N],  v[N], w[N];                                              //
//                                                                            //
//     (your code to initialize the vector u and v)                           //
//                                                                            //
//     Add_Vectors(w, u, v, N);                                               //
//     printf("The vector w = u + v is \n"); ...                              //
////////////////////////////////////////////////////////////////////////////////
int Add_Vectors(double MyVector3[], double MyVector2[], double MyVector1[], unsigned int nrows) 
{
   unsigned int i;
   for (i = 0; i < nrows; i++) MyVector3[i] = MyVector2[i] + MyVector1[i];
   return 0;
}
////////////////////////////////////////////////////////////////////////////////
// File: subtract_vectors.c                                                   //
// Routine(s):                                                                //
//    Subtract_Vectors                                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Subtract_Vectors(double *w, double *u, double *v, int n)             //
//                                                                            //
//  Description:                                                              //
//     Subtract the vector v from the vector u to form the vector w, i.e.     //
//     w[i] = u[i] - v[i].  All vectors should be declared as "double [n]"    //
//     in the calling routine.                                                //
//                                                                            //
//  Arguments:                                                                //
//     double w[]   Resultant vector w = u - v.                               //
//     double u[]   The minuend.                                              //
//     double v[]   The subtrahend.                                           //
//     int    n     The number of components of u, v, and w.                  //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double u[N],  v[N], w[N];                                              //
//                                                                            //
//     (your code to initialize the vector u and v)                           //
//                                                                            //
//     Subtract_Vectors(w, u, v, N);                                          //
//     printf("The vector w = u - v is \n"); ...                              //
////////////////////////////////////////////////////////////////////////////////
int Subtract_Vectors(double MyVector3[], double MyVector2[], double MyVector1[], unsigned int nrows) 
{
   unsigned int i;
   for (i = 0; i < nrows; i++) MyVector3[i] = MyVector2[i] - MyVector1[i];
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
// File: mul_vector_by_scalar.c                                               //
// Routine(s):                                                                //
//    Multiply_Vector_by_Scalar                                               //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Multiply_Vector_by_Scalar(double *v, double x, int n)                //
//                                                                            //
//  Description:                                                              //
//     Multiply the vector v by the scalar x, i.e. multiply each component    //
//     of the vector v by the scalar x, v[i] <- v[i] * x for all i.           //
//                                                                            //
//  Arguments:                                                                //
//     double *v    Pointer to the first element of the vector v.             //
//     double x     Scalar which multiplies each element of the vector v.     //
//     int    n     The number of components of the vector v.                 //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N],  x;                                                       //
//                                                                            //
//     (your code to initialize the vector v and scalar x)                    //
//                                                                            //
//     Multiply_Vector_by_Scalar(v, x, N);                                    //
//     printf("The vector v is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
int Multiply_Vector_by_Scalar(double MyVector1[], double Scalar, unsigned int nrows)
{
  unsigned int i;
  for (i=0; i < nrows; i++)
    MyVector1[i] *= Scalar;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// File: div_vector_by_scalar.c                                               //
// Routine(s):                                                                //
//    Divide_Vector_by_Scalar                                                 //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Divide_Vector_by_Scalar(double *v, double x, int n)                  //
//                                                                            //
//  Description:                                                              //
//     Divide the vector v by the non-zero scalar x, i.e. divide each         //
//     component of the vector v by the scalar x, v[i] <- v[i] / x for all i. //
//                                                                            //
//  Arguments:                                                                //
//     double *v    Pointer to the first element of the vector v.             //
//     double x     Scalar which divides each element of the vector v.        //
//     int    n     The number of components of the vector v.                 //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N],  x;                                                       //
//                                                                            //
//     (your code to initialize the vector v and scalar x)                    //
//                                                                            //
//     if ( x != 0.0)  Divide_Vector_by_Scalar(v, x,N);                       //
//      printf("The vector v is \n"); ...                                     //
////////////////////////////////////////////////////////////////////////////////
extern void Divide_Vector_by_Scalar(double X[], double TempScalar, unsigned int ncols)
{
   double z = 1.0 / TempScalar;

   for (; ncols > 0; ncols--) *X++ *= z;
}

////////////////////////////////////////////////////////////////////////////////
// File: multiply_vector_by_vector.c                                          //
// Routine(s):                                                                //
//    Multiply_Vector_by_Vector                                               //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Multiply_Vector_by_Vector(double *A, int nrows, int ncols,           //
//                                                  double u[], double v[])   //
//                                                                            //
//  Description:                                                              //
//     Post multiply the column vector u by the row vector v to form the      //
//     nrows x ncols matrix A.  Here u is an nrows x 1 column vector and v    //
//     is an 1 x ncols row vector and  A[i][j] = u[i]v[j].                    //
//     The matrix A should be declared as "double A[nrows][ncols]" in the     //
//     calling routine.  The vector v declared as "double v[ncols]" and       //
//     the vector u declared as "double u[nrows]" in the calling routine.     //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of the matrix A and the number of      //
//                  components of the column vector u.                        //
//     int    ncols The number of columns of the matrices A and the           //
//                  number of components of the row vector v.                 //
//     double *u    Pointer to the first element of the vector u.             //
//     double *v    Pointer to the first element of the vector v.             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  u[M], v[N];                                           //
//                                                                            //
//     (your code to initialize the column vector u and row vector v)         //
//                                                                            //
//     Multiply_Vector_by_Vector(&A[0][0], M, N, u, v);                       //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
int Multiply_Vector_by_Vector(vector<double> A[], unsigned int nrows, unsigned int ncols, double MyVector1[], double MyVector2[])
{
   unsigned int i,j;
   for (i = 0; i < nrows; i++) {
      for (j = 0; j < ncols; j++) A[i][j] = MyVector1[i] * MyVector2[j];
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
// File: zero_vector.c                                                        //
// Routine(s):                                                                //
//    Zero_Vector                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Zero_Vector(double *A, int n)                                        //
//                                                                            //
//  Description:                                                              //
//     Set the vector A equal to the zero vector, i.e. A[i] = 0 for all i.    //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the vector A.             //
//     int    n     The number of components of the vector A.                 //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N];                                                           //
//                                                                            //
//     Zero_Vector(A, N);                                                     //
//     printf("The vector A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
int Zero_Vector(double MyVector1[], unsigned int nrows)
{
  unsigned int i;
  for (i = 0; i < nrows; i++)
    MyVector1[i] = 0.0;
  return 0;
}
////////////////////////////////////////////////////////////////////////////////
// File: canonical_basis_vector.c                                             //
// Routine(s):                                                                //
//    Canonical_Basis_Vector                                                  //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Canonical_Basis_Vector(double v[], int i, int n)                     //
//                                                                            //
//  Description:                                                              //
//     Set the n dimension vector v to the i_th canonical basis vector, i.e.  //
//     v[j] = 0 if j != i and v[i] = 1, j = 0,...,n-1, i = 0,...,n-1.         //
//     If n <= 1, v[0] is set to 1 regardless of i. If n > 1 and i >= n,      //
//     or i < 0, then v[i] is not set to 1, i.e. the vector v is set to the   //
//     zero vector.                                                           //
//                                                                            //
//  Arguments:                                                                //
//     double v[]   Vector of dimension n, upon return each component is 0.0  //
//                  except v[i] = 1.0.  If n <= 1, v[0] is set to 1.0.        //
//                  If i >= n or i < 0 then v[i] is not set.                  //
//     int    i     The component of v[] set to 1.0, other components are 0.0 //
//     int    n     The number of components of the vector v.                 //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N];                                                           //
//     int i;                                                                 //
//                                                                            //
//     Canonical_Basis_Vector(v, i, N);                                       //
//     printf("The ith canonical basis vector is \n"); ...                    //
////////////////////////////////////////////////////////////////////////////////
int Canonical_Basis_Vector(double MyVector1[], unsigned int i, unsigned int nrows)
{
   unsigned int j;
   if ( nrows <= 1) { MyVector1[0] = 1.0; return 1; }
   for (j = 0; j < nrows; j++) MyVector1[j] = 0.0;
   if ( (i >= 0) && (i < nrows) ) MyVector1[i] = 1.0;
   return 0;
}

extern int TransposeMatrix(vector<double> A[], unsigned int nrows, unsigned int ncols)
{
  unsigned int i,j;
  double elem {0.0};
  matrix<double> transpose;
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
//  PrintVector2D(&A[0], ncols);
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
//  PrintVector2D(&A[0], ncols);
  transpose.clear();
  transpose.shrink_to_fit();	
  return 0;
}

extern int TransposeSquareMatrix(vector<double> A[], unsigned int nrows, unsigned int ncols)
{
  unsigned int i,j;
  double elem {0.0};
  matrix<double> transpose;
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
//    PrintVector2D(&A[0], ncols);
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
//  PrintVector2D(&A[0], ncols);
    transpose.clear();
    transpose.shrink_to_fit();	
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// File: inner_product.c                                                      //
// Routines:                                                                  //
//    Inner_Product                                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Inner_Product(double u[], double v[], int n)                       //
//                                                                            //
//  Description:                                                              //
//     An inner product on a real vector space V is a positive definite       //
//     symmetric bilinear form < , > : V x V -> R.                            //
//     Given a Real vector space V with a basis e[i] for which the inner      //
//     product of the basis vectors <e[i],e[j]> = delta(i,j), delta(i,j) being//
//     the Kronecker delta function, the inner product of two vectors in V    //
//     is the sum of the component-wise products.  If dim V = 3,              //
//     u = u[0] e[0] + ... + u[n-1] e[n-1] and                                //
//     v = v[0] e[0] + ... + v[n-1] e[n-1] are vectors in V, then             //
//     <u,v> = <u[0]e[0]+...+u[n-1]e[n-1], v[0]e[0]+...+u[n-1]e[n-1]>         //
//          = u[0]v[0] <e[0],e[0]> + ... + u[0]v[n-1] <e[n-1],e[n-1]>         //
//           + ... +                                                          //
//           u[n-1]v[0] <e[n-1,e[0]> + ... + u[n-1]v[n-1] <e[n-1],e[n-1]>     //
//          =  u[0]v[0] + ... + u[n-1]v[n-1]                                  //
//                                                                            //
//     The arguments u and v should be declared as double u[N] and            //
//     double v[N] where N >= n in the calling program.                       //
//                                                                            //
//  Arguments:                                                                //
//     double u[]  Pointer to the first element of the vector u.              //
//     double v[]  Pointer to the first element of the vector v.              //
//     int     n   The number of components of the vectors u and v.           //
//                                                                            //
//  Return Values:                                                            //
//     Inner Product of u and v.                                              //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double u[N], v[N], inner_product;                                      //
//                                                                            //
//     (your code to intialize the vectors u and v)                           //
//     inner_product = Inner_Product(u,v,N);                                  //
//     printf(" <u,v> = %12.6f\n", inner_product);                            //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
extern double Inner_Product(double U[], double V[], unsigned int nrows)
{
   double Inner_Product = 0.0;
   for (nrows--; nrows > 0; nrows--) Inner_Product +=  U[nrows] * V[nrows];
   return Inner_Product;
}

////////////////////////////////////////////////////////////////////////////////
// File: vector_max_norm.c                                                    //
// Routines:                                                                  //
//    Vector_Max_Norm                                                         //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Vector_Max_Norm(double v[], int n)                                 //
//                                                                            //
//  Description:                                                              //
//     Given an n-d real vector space V with a basis consisting of vectors    //
//     with unit norm, the max norm on the vector space V is the maximum of   //
//     the absolute values of the components of a vector v with respect to    //
//     that basis i.e. for v = v[0]e[0] + ... + v[n-1]e[n-1], where e[0], ...,//
//     e[n-1] are the basis vectors for which                                 //
//     || e[0] || = ... = || e[n-1] || = 1, then                              //
//                   || v || = Max( |v[0]|, ..., |v[n-1]| ).                  //
//                                                                            //
//  Arguments:                                                                //
//     double v[]  Pointer to the first element of the vector v[n].           //
//     int     n   The number of elements of the vector v[].                  //
//                                                                            //
//  Return Values:                                                            //
//     max norm of v.                                                         //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N], norm;                                                     //
//                                                                            //
//     (your code to initialize the vector v)                                 //
//     norm = Vector_Max_Norm(v, N);                                          //
//     printf(" || v || = %12.6f\n", norm);                                   //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
extern double Vector_Max_Norm(double TempVector[], unsigned int ncols)
{
   double norm = 0.0;
   double TempScalar;
   unsigned int i;

   for (i = 0; i < ncols; i++) if (norm < ( TempScalar = fabs( TempVector[i] ) ) ) norm = TempScalar;

   return norm;
}

////////////////////////////////////////////////////////////////////////////////
// File: vector_l2_norm.c                                                     //
// Routines:                                                                  //
//    Vector_L2_Norm                                                          //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  double Vector_L2_Norm(double v[], int n)                                  //
//                                                                            //
//  Description:                                                              //
//     Given an n-d real vector space V with a basis consisting of vectors    //
//     with unit norm, the l2 norm on the vector space V is the sqrt of the   //
//     sum of squares of the components of a vector v with respect to that    //
//     basis i.e.                                                             //
//     for v = v[0]e[0] + ... + v[n-1]e[n-1], where e[0], ..., e[n-1] are     //
//     the basis vectors for which || e[0] || = ... = || e[n-1] || = 1,       //
//     then                                                                   //
//               || v || = sqrt( |v[0]|^2  + ... + |v[n-1]|^2 ).              //
//                                                                            //
//  Arguments:                                                                //
//     double v[]  Pointer to the first element of the vector v[n].           //
//     int     n   The number of elements of the vector v[].                  //
//                                                                            //
//  Return Values:                                                            //
//     l2 norm of v.                                                          //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N], norm;                                                     //
//                                                                            //
//     (your code to intialize the vector v)                                  //
//     norm = Vector_L2_Norm(v, N);                                           //
//     printf(" || v || = %12.6f\n", norm);                                   //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
extern double Vector_L2_Norm(double v[], unsigned int ncols)
{
   double norm = 0.0;
   unsigned int i;
   for (i = 0; i < ncols; i++) norm +=  v[i] * v[i];
   return sqrt(norm);
}

#endif
