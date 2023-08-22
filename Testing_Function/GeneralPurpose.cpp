#ifndef MAKE_GPP_FUNC
#define MAKE_GPP_FUNC

#include <string.h>                                 // required for memcpy()
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cfloat>
#include <math.h>		// required for fabs() and for sqrt()

using namespace std;

namespace {
	
#include "GeneralPurpose.h"

#define a(i,j) a[(i)*nrows+(j)]

template<typename T>
using matrix = std::vector< std::vector<T> >;

/*
Iterate over vector of vectors and for each of the
nested vector print its contents
*/
static void PrintVector2D(const vector<double> TempMatrix[], unsigned int nrows)
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
static void PrintVector(const double TempVector[], unsigned int nrows){
  std::cout << "Displaying the vector: " << endl;
  std::cout << std::setprecision(5);
  for (unsigned int i = 0; i < nrows; i++)
    std::cout << TempVector[i] << "   ";
  std::cout << endl;
}

int Testing_General_Purpose_Parameters(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols)
{
  unsigned int i,j;
  double elem {0.0};
  double Scalar {2.5};
  (void) Scalar;
  matrix<double> A_Copy;
    // Inserting elements into vector
    for (i = 0; i < nrows; i++) {
      // Vector to store column elements
      vector<double> Row2;
        for (j = 0; j < ncols; j++) {
		  elem = A[i][j];
          Row2.push_back(elem);
            if (j == (ncols)) {
				Row2.push_back(elem);
			}
        }
        // Pushing back above 1D vector
        // to create the 2D vector
        A_Copy.push_back(Row2);
    }
        
//	std::cout << "The matrix A is the following: " << endl;
//	PrintVector2D(&A[0], nrows);
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
//  std::cout << "Displaying the vector: " << endl;
//  std::cout << std::setprecision(5);
//  for (unsigned int i = 0; i < nrows; i++)
//    std::cout << X[i] << "   ";
//  std::cout << endl;

//	TransposeSquareMatrix (&A[0], nrows, ncols);
//	TransposeMatrix (&A[0], nrows, ncols);

//    Output_Results (&X[0], &B[0], nrows, ncols);

//  Copy_Matrix(&A[0], &A_Copy[0], nrows, ncols);
//  PrintVector2D(&A_Copy[0], ncols);
//  std::cout << "Displaying the 2D vector:" << endl;
/*
Iterate over vector of vectors and for each of the
nested vector print its contents
*/
// Displaying the 2D vector
//  std::cout << std::setprecision(5);
//  for (unsigned int i = 0; i < nrows; i++) {
//    for ( 
//    auto it = A[i].begin();
//      it != A[i].end(); it++)
//      cout << *it << " ";
//      cout << endl;
//  }

    vector <double> MyVector1;
    vector <double> MyVector2;
    unsigned int Vec_Size {0};
    for (i=0; i < nrows; i++)
    {
      MyVector1.push_back(100.0);
      MyVector2.push_back(200.0);
      Vec_Size++;
	}
    unsigned int kth_row {1};
    (void) kth_row;
    unsigned int kth_column {1};
    (void) kth_column;
//    (void) kth_row and kth_column;	
//    PrintVector2D(&A[0], nrows);
//    Get_Row(&MyVector1[0], &A[0], kth_row, nrows, ncols);
//    std::cout << "The vector MyVector1 is the following: " << endl;
//    PrintVector(&MyVector1[0], nrows);
//    Get_Column(&MyVector2[0], &A[0], kth_column, nrows, ncols);
//    std::cout << "The vector MyVector2 is the following: " << endl;
//    PrintVector(&MyVector2[0], nrows);
//    for (i = 0; i < nrows; i++) {
//	  MyVector1[i] = MyVector1[i] + 100;
//	  MyVector2[i] = MyVector2[i] + 100;  
//    }
//    Set_Row(&MyVector1[0], &A[0], kth_row, nrows, ncols);
//    Set_Column(&MyVector2[0], &A[0], kth_row, nrows, ncols);
//    PrintVector2D(&A[0], nrows);
//    Set_Diagonal_to_Scalar(&A[0], Scalar, nrows, ncols);
//    PrintVector2D(&A[0], nrows);
//    Scalar = 110.0;
//    Fill_Matrix_with_Scalar(&A[0], Scalar, nrows, ncols);
//    PrintVector2D(&A[0], nrows);
//	PrintVector2D(&A[0], nrows);
//    Add_Scalar_to_Diagonal(&A[0], Scalar, nrows, ncols);
//	PrintVector2D(&A[0], nrows);
//	Testing_Function_Parameters (&A[0], &X[0], &B[0], nrows, ncols);
//    Output_Results (&X[0], &B[0], nrows, ncols);
//    Test_Vec_Swap(&MyVector1[0], &MyVector2[0], Vec_Size);
//    PrintVector2D(&A[0], nrows);
//    Get_Diagonal(&X[0], &A[0], nrows, ncols);
//    PrintVector(&X[0], nrows);
//    for (i = 0; i < nrows; i++)
//      X[i] += 100.0;
//    Set_Diagonal(&X[0], &A[0], nrows, ncols);
//    PrintVector2D(&A[0], nrows);    
//    unsigned int row1 {0}, row2 {1};
//    PrintVector2D(&A[0], nrows);
//    Interchange_Rows(&A[0], row1, row2, ncols);
//    PrintVector2D(&A[0], nrows);    
//    unsigned int col1 {1}, col2 {2};
//    PrintVector2D(&A[0], nrows);
//    Interchange_Columns(&A[0], col1, col2, nrows, ncols);
//    PrintVector2D(&A[0], nrows);
    // Initializing the vector of vectors
//    unsigned int mrows {2}, mcols {2}, row {0}, col {0};
//    vector<vector<double> > A;
//    matrix<double> S;
    // Inserting elements into vector
//    for (i = 0; i < mrows; i++) {
      // Vector to store column elements
//      vector<double> v3;
//        for (j = 0; j < mcols; j++) {
//		  elem = 0.0;
//          v3.push_back(elem);
//           if (j == (mcols)) {
//				v3.push_back(elem);
//			}
//        }
        // Pushing back above 1D vector
        // to create the 2D vector
//        S.push_back(v3);
//    }
//    PrintVector2D(&A[0], nrows);
//    Get_Submatrix(&A[0], nrows, ncols, &S[0], mrows, mcols, row, col);
//    PrintVector2D(&S[0], mrows);
//    Set_Submatrix(&A[0], nrows, ncols, &S[0], mrows, mcols, row, col, Scalar);
//    PrintVector2D(&A[0], nrows); 
//  matrix<double> BMat;
  // Inserting elements into vector
//  for (i = 0; i < nrows; i++) {
    // Vector to store column elements
//    vector<double> Row3;
//    for (j = 0; j < ncols; j++) {
// Change the BMatrix elements here.
//	  elem = A[i][j];
//      Row3.push_back(elem);
//      if (j == (ncols)) {
//	    Row3.push_back(elem);
//	  }
//    }
    // Pushing back above 1D vector
    // to create the 2D vector
//    BMat.push_back(Row3);
//  }
  // Join_Rows CMatrix initialization
//  matrix<double> CMat;
//  int mcols = ncols;
    // Inserting elements into vector
//  for (i = 0; i < nrows; i++) {
    // Vector to store column elements
//    vector<double> Row4;
//    for (j = 0; j < (ncols+mcols-1); j++) {
//	  elem = 0.0;
//      Row4.push_back(elem);
//      if (j == (ncols)) {
//	    Row4.push_back(elem);
//	  }
//    }
    // Pushing back above 1D vector
    // to create the 2D vector
//    CMat.push_back(Row4);
//  }
  // Join_Columns CMatrix initialization
//  matrix<double> CMat;
//  int mcols = ncols;
    // Inserting elements into vector
//  for (i = 0; i < nrows+mcols; i++) {
    // Vector to store column elements
//    vector<double> Row4;
//    for (j = 0; j < ncols; j++) {
//	  elem = 0.0;
//      Row4.push_back(elem);
//      if (j == (ncols)) {
//	    Row4.push_back(elem);
//	  }
//    }
    // Pushing back above 1D vector
    // to create the 2D vector
//    CMat.push_back(Row4);
//  }

//  Join_Rows(&A[0], &BMatrix[0], &CMatrix[0], nrows, ncols, mcols);
//  PrintVector2D(&CMatrix[0], nrows+mcols);

//  std::cout << std::setprecision(5);
//  for (unsigned int i = 0; i < nrows; i++) {
//    for (auto it = CMatrix[i].begin();
//        it != CMatrix[i].end(); it++)
//        cout << *it << " ";
//      cout << endl;
//  }
  
//  Join_Columns(&A[0], &BMatrix[0], &CMatrix[0], nrows, ncols, mcols);
//  PrintVector2D(&CMatrix[0], nrows+mcols);
    MyVector1.clear();
    MyVector1.shrink_to_fit();
    MyVector2.clear();
    MyVector2.shrink_to_fit();
    A_Copy.clear();
    A_Copy.shrink_to_fit();
//    BMat.clear();
//    BMat.shrink_to_fit();
//    CMat.clear();
//    CMat.shrink_to_fit(); 
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
      P.push_back(0);
    }
    fclose (myFile2);
    // Initializing the vector of vectors
//    vector<vector<double> > A;
    matrix<double> A;
    // Inserting elements into vector
    for (i = 0; i < nrows; i++) {
      // Vector to store column elements
      vector<double> v1;
        for (j = 0; j < ncols; j++) {
		  elem = array_A[i*nrows+j];
          v1.push_back(elem);
            if (j == (ncols)) {
				v1.push_back(elem);
			}
        }
        // Pushing back above 1D vector
        // to create the 2D vector
        A.push_back(v1);
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
    array_A.clear();
    array_A.shrink_to_fit();
    B.clear();
    B.shrink_to_fit();
    X.clear();
    X.shrink_to_fit();
    A.clear();
    A.shrink_to_fit();
    Augmented.clear();
    Augmented.shrink_to_fit();
    P.clear();
    P.shrink_to_fit();

//    std::cin.clear(); // reset any error flags
//    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//    ignore any characters in the input buffer until we find an enter character
//    std::cin.get(); // get one more char from the user
    return 0;
  }
};

}  // close namespace

////////////////////////////////////////////////////////////////////////////////
// File: copy_matrix.c                                                        //
// Routine(s):                                                                //
//    Copy_Matrix                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Copy_Matrix(double *A, double *B, int nrows, int ncols)              //
//                                                                            //
//  Description:                                                              //
//     Copy the nrows x ncols matrix B to the nrows x ncols matrix A.         //
//     i.e.    A = B.                                                         //
//     The memory locations of the source matrix, B, and the destination      //
//     matrix, A, must not overlap, otherwise the results are installation    //
//     dependent.                                                             //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the destination matrix    //
//                  A[nrows][ncols].                                          //
//     double *B    Pointer to the first element of the source matrix         //
//                  B[nrows][ncols].                                          //
//     int    nrows The number of rows matrices A and B.                      //
//     int    ncols The number of columns of the matrices A and B.            //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[N][M],  B[N][M];                                              //
//                                                                            //
//     (your code to initialize the matrix B)                                 //
//                                                                            //
//     Copy_Matrix(&A[0][0], &B[0][0], N, M);                                 //
//     printf(" Matrix A is \n");                                             //
////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////
// File: copy_vector.c                                                        //
// Routine(s):                                                                //
//    Copy_Vector                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Copy_Vector(double *d, double *s, int n)                             //
//                                                                            //
//  Description:                                                              //
//     Copy the n dimensional vector s(source) to the n dimensional           //
//     vector d(destination).  The memory locations of the source and         //
//     destination vectors must not overlap, otherwise the results            //
//     are installation dependent.                                            //
//                                                                            //
//  Arguments:                                                                //
//      double *d  Pointer to the first element of the destination vector d.  //
//      double *s  Pointer to the first element of the source vector s.       //
//      int    n   The number of elements of the source / destination vectors.//
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double v[N],  vd[N];                                                   //
//                                                                            //
//     (your code to initialize the vector v)                                 //
//                                                                            //
//     Copy_Vector(vd, v, N);                                                 //
//     printf(" Vector vd is \n");                                            //
////////////////////////////////////////////////////////////////////////////////
 extern void Copy_Vector(double *d, double *s, unsigned int ncols)
{
   memcpy(d, s, sizeof(double) * ncols);
}
////////////////////////////////////////////////////////////////////////////////
// File: get_row.c                                                            //
// Routine(s):                                                                //
//    Get_Row                                                                 //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Get_Row(double v[], double *A, int row, int ncols)                   //
//                                                                            //
//  Description:                                                              //
//     Copy the row 'row' from the nrows x ncols matrix A to the vector v.    //
//     Note that v should be declared "double v[N]", with N >= ncols in the   //
//     calling routine.                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double v[]   Destination address of the row 'row' of the matrix A.     //
//     double *A    Pointer to the first element of the matrix A[*][ncols].   //
//     int    row   The row of A to copy to the vector v, row = 0,...,nrows-1.//
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[N];                                                 //
//                                                                            //
//     (your code to create the matrix A and the row number i)                //
//                                                                            //
//     if ( (i > 0) && ( i < M ) ) Get_Row(v, &A[0][0], i, N);                //
//     printf("The vector v is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Get_Row(double X[], vector <double> A[], unsigned int kth_row, unsigned int nrows, unsigned int ncols)
{
  unsigned int j;
  for (j = 0; j < ncols; j++)
    X[j] = A[kth_row][j];
}
////////////////////////////////////////////////////////////////////////////////
// File: get_column.c                                                         //
// Routine(s):                                                                //
//    Get_Column                                                              //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Get_Column(double v[], double *A, int col, int nrows, int ncols)     //
//                                                                            //
//  Description:                                                              //
//     Copy the column 'col' from the nrows x ncols matrix A to the vector v. //
//     Note that v should be declared "double v[N]", with N >= nrows in the   //
//     calling routine.                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double v[]   Destination address of the column col of the matrix A.    //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    col   The column of matrix A to copy to the vector v,           //
//                  0 <= col < ncols.                                         //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[N];                                                 //
//                                                                            //
//     (your code to set the matrix A and the column number i)                //
//                                                                            //
//     if ( (i >= 0) && (i < N) ) Get_Column(v, &A[0][0], i, M, N);           //
//     printf("The vector v is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Get_Column(double X[], vector <double> A[], unsigned int kth_column, unsigned int nrows, unsigned int ncols)
{
  unsigned int i;
  for (i = 0; i < nrows; i++)
    X[i] = A[i][kth_column];
}
////////////////////////////////////////////////////////////////////////////////
// File: set_row.c                                                            //
// Routine(s):                                                                //
//    Set_Row                                                                 //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Set_Row(double *A, double v[], int row, int ncols)                   //
//                                                                            //
//  Description:                                                              //
//     Copy the vector v to the row 'row' in the nrows x ncols matrix A.      //
//     Note that v should be declared double v[N], N >= ncols in the calling  //
//     routine.                                                               //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     double v[]   Source address of the row to replace the row in the       //
//                  matrix A.                                                 //
//     int    row   The row of A in which the vector v is to be copied,       //
//                  0 <= row < nrows.                                         //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[N];                                                 //
//                                                                            //
//     (your code to initialize the matrix A, vector v, and the row number i) //
//                                                                            //
//     if ( (i >= 0) && (i < M) )                                             //
//        Set_Row(&A[0][0], v, i, N);                                         //
//     printf("The matrix A is \n"); ... }                                    //
////////////////////////////////////////////////////////////////////////////////
void Set_Row(double X[], vector <double> A[], unsigned int kth_row, unsigned int nrows, unsigned int ncols)
{
  unsigned int j;
  for (j = 0; j < ncols; j++)
    A[kth_row][j] = X[j];
}
////////////////////////////////////////////////////////////////////////////////
// File: set_column.c                                                         //
// Routine(s):                                                                //
//    Set_Column                                                              //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Set_Column(double *A, double v[], int col, int nrows, int ncols)     //
//                                                                            //
//  Description:                                                              //
//     Copy the vector v to the column 'col' in the nrows x ncols matrix A.   //
//     Note that v should be declared "double v[N]", with N >= nrows in the   //
//     calling routine.                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     double v[]   Source address of the column col to replace the column in //
//                  the matrix A.                                             //
//     int    col   The column of matrix A to copy the vector v,              //
//                  0 <= col < ncols.                                         //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[M];                                                 //
//     int i;                                                                 //
//                                                                            //
//     (your code to initialize the matrix A, vector v, and the               //
//      column number i)                                                      //
//                                                                            //
//     if ( (i >= 0) && (i < N) )                                             //
//        Set_Column(&A[0][0], v, i, M, N);                                   //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Set_Column(double X[], vector <double> A[], unsigned int kth_column, unsigned int nrows, unsigned int ncols)
{
  unsigned int i;
  for (i = 0; i < nrows; i++)
    A[i][kth_column] = X[i];
}
////////////////////////////////////////////////////////////////////////////////
// File: get_diagonal.c                                                       //
// Routine(s):                                                                //
//    Get_Diagonal                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Get_Diagonal(double v[], double *A, int nrows, int ncols)            //
//                                                                            //
//  Description:                                                              //
//     Copy the diagonal of the matrix A to the vector v, i.e.                //
//     v[i] =  A[i][i], i = 0, ..., min( nrows, ncols ).                      //
//     Note that v should be declared "double v[N]", with                     //
//     N >= min(nrows, ncols) in the calling routine.                         //
//                                                                            //
//  Arguments:                                                                //
//     double v[]   Destination address of the diagonal of the matrix A.      //
//     double *A    Pointer to the first element of the source matrix A.      //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[N];                                                 //
//                                                                            //
//     (your code to set the matrix A )                                       //
//                                                                            //
//     Get_Diagonal(v, &A[0][0],  M, N);                                      //
//     printf("The diagonal is \n"); ... }                                    //
////////////////////////////////////////////////////////////////////////////////
void Get_Diagonal(double X[], vector <double> A[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i, j, n;
   n = (nrows < ncols) ? nrows: ncols;
   for (i=0; i < n; i++)
     for (j=0; j < n; j++)
       X[i] = A[i][i];
//   for (i = 0; i < n; A += (ncols + 1), i++) v[i] = *A;
}
////////////////////////////////////////////////////////////////////////////////
// File: set_diagonal.c                                                       //
// Routine(s):                                                                //
//    Set_Diagonal                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Set_Diagonal(double *A, double v[], int nrows, int ncols)            //
//                                                                            //
//  Description:                                                              //
//     Copy the vector v to the diagonal of the matrix A, i.e.                //
//     A[i][i] = v[i], where 0 <= i < min( nrows, ncols ).                    //
//     Note that v should be declared "double v[N]", N >= min( nrows, ncols ) //
//     in the calling routine.                                                //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the destination matrix A. //
//     double v[]   Source of the new diagonal for the matrix A.              //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N],  v[N];                                                 //
//                                                                            //
//     (your code to initialize the matrix A and the vector v)                //
//                                                                            //
//     Set_Diagonal(&A[0][0], v, M, N);                                       //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Set_Diagonal(double X[], vector <double> A[], unsigned int nrows, unsigned int ncols)
{
   unsigned int i, j, n;
   if (nrows < ncols) n = nrows; else n = ncols;
   for (i=0; i < n; i++)
     for (j=0; j < n; j++)
       A[i][i] = X[i];
//   for (; n > 0 ; A += (ncols + 1), n--)  *A = *v++;
}
////////////////////////////////////////////////////////////////////////////////
// File: set_diagonal_to_scalar.c                                             //
// Routine(s):                                                                //
//    Set_Diagonal_to_Scalar                                                  //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Set_Diagonal_to_Scalar(double *A, double x, int nrows, int ncols)    //
//                                                                            //
//  Description:                                                              //
//     Replace each element of the diagonal A[i][i], where 0 <= i <           //
//     min( nrows, ncols ) with the scalar x.                                 //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the destination matrix A. //
//     double x     Scalar which replaces the diagonal elements of A.         //
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
//     (your code to create the matrix A and the scalar x)                    //
//                                                                            //
//     Set_Diagonal_to_Scalar(&A[0][0], x, M, N);                             //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Set_Diagonal_to_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols)
{
   unsigned int i, n  = (nrows < ncols) ? nrows : ncols;
   for (i = 0; i < n; i++)
     A[i][i] = Scalar;
}
////////////////////////////////////////////////////////////////////////////////
// File: get_submatrix.c                                                      //
// Routine(s):                                                                //
//    Get_Submatrix                                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Get_Submatrix(double *S, int mrows, int mcols,                       //
//                                   double *A, int ncols, int row, int col)  //
//                                                                            //
//  Description:                                                              //
//     Copy the mrows and mcols of the nrows x ncols matrix A starting with   //
//     A[row][col] to the submatrix S.                                        //
//     Note that S should be declared double S[mrows][mcols] in the calling   //
//     routine.                                                               //
//                                                                            //
//  Arguments:                                                                //
//     double *S    Destination address of the submatrix.                     //
//     int    mrows The number of rows of the matrix S.                       //
//     int    mcols The number of columns of the matrix S.                    //
//     double *A    Pointer to the first element of the matrix A[nrows][ncols]//
//     int    ncols The number of columns of the matrix A.                    //
//     int    row   The row of A corresponding to the first row of S.         //
//     int    col   The column of A corresponding to the first column of S.   //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     #define NB                                                             //
//     #define MB                                                             //
//     double A[M][N],  B[MB][NB];                                            //
//     int row, col;                                                          //
//                                                                            //
//     (your code to set the matrix A, the row number row and column number   //
//      col)                                                                  //
//                                                                            //
//     if ( (row >= 0) && (col >= 0) && ((row + MB) < M) && ((col + NB) < N) )//
//        Get_Submatrix(&B[0][0], MB, NB, &A[0][0], N, row, col);             //
//     printf("The submatrix B is \n"); ... }                                 //
////////////////////////////////////////////////////////////////////////////////
void Get_Submatrix(vector<double> A[], unsigned int nrows, unsigned int ncols, vector<double> S[], unsigned int mrows, unsigned int mcols, unsigned int row, unsigned int col)
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
////////////////////////////////////////////////////////////////////////////////
// File: set_submatrix.c                                                      //
// Routine(s):                                                                //
//    Set_Submatrix                                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Set_Submatrix(double *A, int ncols, double *S, int mrows, int mcols, //
//                                                        int row, int col)   //
//                                                                            //
//  Description:                                                              //
//     Copy the mrows x mcols submatrix S into the nrows x ncols matrix A     //
//     starting at the location A[row][col].                                  //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A[n][n].       //
//     int    ncols The number of columns of the matrix A.                    //
//     double *S    Destination address of the submatrix.                     //
//     int    mrows The number of rows matrix S.                              //
//     int    mcols The number of columns of the matrix S.                    //
//     int    row   The row of A corresponding to the first row of S.         //
//     int    col   The column of A corresponding to the first column of S.   //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     #define NB                                                             //
//     #define MB                                                             //
//     double A[M][N],  B[MB][NB];                                            //
//     int row, col;                                                          //
//                                                                            //
//     (your code to initialize the matrix A, submatrix B, row number row,    //
//      and column number col )                                               //
//                                                                            //
//     if ( (row > 0) && ( row + MB < M ) && ( col > 0 ) && (col + NB < N)    //
//        Set_Submatrix(&A[0][0], N, &B[0][0], MB, NB, row, col);             //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Set_Submatrix(vector<double> A[], unsigned int nrows, unsigned int ncols, vector<double> S[], unsigned int mrows, unsigned int mcols, unsigned int row, unsigned int col, double Scalar)
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
////////////////////////////////////////////////////////////////////////////////
// File: fill_matrix_with_scalar.c                                            //
// Routine(s):                                                                //
//    Fill_Matrix_with_Scalar                                                 //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Fill_Matrix_with_Scalar(double *A, double x, int nrows, int ncols)   //
//                                                                            //
//  Description:                                                              //
//     Set each element of the matrix A to the scalar x.                      //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the destination matrix A. //
//     double x     Scalar which replaces each element of the matrix A.       //
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
//     (your code to initialize the matrix A and the scalar x)                //
//                                                                            //
//     Fill_Matrix_with_Scalar(&A[0][0], x, M, N);                            //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Fill_Matrix_with_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols)
{
   unsigned int i,j;
   for (i = 0; i < nrows; i++)
      for ( j = 0; j < ncols; j++) A[i][j] = Scalar;
}

////////////////////////////////////////////////////////////////////////////////
// File: join_by_rows.c                                                       //
// Routine(s):                                                                //
//    Join_Matrices_by_Rows                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Join_Matrices_by_Rows(double *C, double *A, int nrows, int ncols,    //
//                                                     double *B, int mcols)  //
//                                                                            //
//  Description:                                                              //
//     Copy the nrows x ncols matrix A into the nrows x (ncols + mcols)       //
//     matrix C and then copy the nrows x mcols matrix B to starting at       //
//     the location C[0][ncols], i.e. C = [A:B].                              //
//     The matrix C should be declared as double C[nrows][ncols + mcols] in   //
//     the calling routine.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *C    Pointer to the first element of the matrix C.             //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    nrows The number of rows of the matrices A and B.               //
//     int    ncols The number of columns of the matrix A.                    //
//     double *B    Pointer to the first element of the matrix B.             //
//     int    mcols The number of columns of the matrix B.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     #define NB                                                             //
//     double A[M][N],  B[M][NB], C[M][N+NB];                                 //
//                                                                            //
//     (your code to initialize the matrices A and B)                         //
//                                                                            //
//     Join_Matrices_by_Rows(&C[0][0], &A[0][0], M, N, &B[0][0], NB);         //
//     printf("The matrix C is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
int Join_Rows(vector<double> A[], vector<double> BMatrix[], vector <double> CMatrix[],  unsigned int nrows, unsigned int ncols, unsigned int mcols)
{
  unsigned int i, j;
  // Print C //
//  printf("C Before\n");
//  for (i = 0; i < nrows; i++) {
//    for (j = 0; j < (ncols+mcols); j++) printf("%6.1f",CMatrix[i][j]);
//      printf("\n");
//  }
//  printf("\n\n");	
  for (i = 0; i < nrows; i++) {
    for (j = 0; j < ncols+mcols; j++)
      if (j < ncols) { 
        CMatrix[i][j] = A[i][j];
      }
      else {
	    CMatrix[i][j] = BMatrix[i][j-ncols];	  
	  }
  }
//  printf("C After\n");
//  for (i = 0; i < nrows; i++) {
//    for (j = 0; j < (ncols+mcols); j++) printf("%6.1f",CMatrix[i][j]);
//      printf("\n");
//  }
//  printf("\n\n");
  return 0;	
}

int Join_Columns(vector<double> A[], vector<double> BMatrix[], vector <double> CMatrix[],  unsigned int nrows, unsigned int ncols, unsigned int mcols)
{
  unsigned int i, j;
  // Print C //
//  printf("C Before\n");
//  for (i = 0; i < (nrows+mcols); i++) {
//    for (j = 0; j < ncols; j++) printf("%6.1f",CMatrix[i][j]);
//      printf("\n");
//  }
//  printf("\n\n");	
  for (i = 0; i < (nrows+mcols); i++) {
    for (j = 0; j < ncols; j++)
      if (i < nrows) { 
        CMatrix[i][j] = A[i][j];
      }
      else {
	    CMatrix[i][j] = BMatrix[i-nrows][j];	  
	  }
  }
//  printf("C After\n");
//  for (i = 0; i < (nrows+mcols); i++) {
//    for (j = 0; j < ncols; j++) printf("%6.1f",CMatrix[i][j]);
//      printf("\n");
//  }
//  printf("\n\n");
  return 0;	
}
////////////////////////////////////////////////////////////////////////////////
// File: interchange_rows.c                                                   //
// Routine(s):                                                                //
//    Interchange_Rows                                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Interchange_Rows(double *A, int row1, int row2, int ncols)           //
//                                                                            //
//  Description:                                                              //
//     Interchange the rows 'row1' and 'row2' of the  nrows x ncols matrix A. //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    row1  The row of A which is to be interchanged with row row2.   //
//     int    row2  The row of A which is to be interchanged with row row1.   //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N];                                                        //
//     int i, j;                                                              //
//                                                                            //
//  (your code to initialize the matrix A, the row number i and row number j) //
//                                                                            //
//     if ( (i >= 0) && ( i < M ) && (j > 0) && ( j < M ) )                   //
//        Interchange_Rows(&A[0][0], i, j, N);                                //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Interchange_Rows(vector<double> A[], unsigned int row1, unsigned int row2, unsigned int ncols)
{
   unsigned int i, j;
   double temp;
   for (j = 0; j < ncols; j++)
   {
     temp = A[row1][j];
     A[row1][j] = A[row2][j];
     A[row2][j] = temp;
   }  
}
////////////////////////////////////////////////////////////////////////////////
// File: interchange_cols.c                                                   //
// Routine(s):                                                                //
//    Interchange_Columns                                                     //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Interchange_Columns(double *A, int col1, int col2, int nrows,        //
//                                                                 int ncols) //
//                                                                            //
//  Description:                                                              //
//     Interchange the columns 'col1' and 'col2' of the  nrows x ncols        //
//     matrix A.                                                              //
//                                                                            //
//  Arguments:                                                                //
//     double *A    Pointer to the first element of the matrix A.             //
//     int    col1  The column of A which is to be interchanged with col2.    //
//     int    col2  The column of A which is to be interchanged with col1.    //
//     int    nrows The number of rows matrix A.                              //
//     int    ncols The number of columns of the matrix A.                    //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     #define M                                                              //
//     double A[M][N];                                                        //
//     int i,j;                                                               //
//                                                                            //
//     (your code to initialize the matrix A, the column number i and column  //
//       number j)                                                            //
//                                                                            //
//     if ( (i >= 0) && ( i < N ) && ( j >= 0 ) && (j < N) )                  //
//        Interchange_Columns(&A[0][0], i, j, M, N);                          //
//     printf("The matrix A is \n"); ...                                      //
////////////////////////////////////////////////////////////////////////////////
void Interchange_Columns(vector<double> A[], unsigned int col1, unsigned int col2, unsigned int nrows, unsigned int ncols)
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
#endif
