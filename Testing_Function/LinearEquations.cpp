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

#include "LinearEquations.h"

#define a(i,j) a[(i)*nrows+(j)]

template<typename T>
using matrix = std::vector< std::vector<T> >;


/*
Iterate over vector of vectors and for each of the
nested vector print its contents
*/
void PrintVector2D(const vector<double> A[], unsigned int nrows)
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
      auto it = A[i].begin();
        it != A[i].end(); it++)
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
//    std::cout << "The X vector should be the following: -0.436 0.430 5.12" << endl;

//	Testing_Function (&A[0], &X[0], &B[0], nrows, ncols);
//    Output_Results (&X[0], &B[0], nrows, ncols, tolerance, max_iter);

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

} // close namespace

int Output_Results(double X[], double B[], unsigned int nrows, unsigned int ncols)
{
  unsigned int i;
  std::cout << std::setprecision (5) << endl;
  std::cout << "Displaying Output_Results" << endl;
  std::cout << "The vector X is the following: " << endl;
  PrintVector(&X[0], nrows);
  std::cout << "The vector B is the following: " << endl;
  PrintVector(&B[0], ncols);
  // Create a new file named "C.dat"
  std::ofstream outputFile("C.dat");
  if (outputFile.is_open()) {
    // Write some text into the file
    outputFile << "%%MatrixMarket_Output_vector_C.dat matrix coordinate pattern general" << endl;
    outputFile << 1 << " " << nrows << " " << nrows;
    outputFile << endl;
    for (i = 0; i < ncols; i++)
      outputFile <<  1 << " " << i+1 << " " << X[i] << endl;
    // Close the file
    outputFile.close();
  } else {
    std::cout << "Error writing to file." << std::endl;
  }
  return 0;
}

int Output_Data(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter)
{
  unsigned int i;
  (void) max_iter;
  (void) tolerance;
  (void) eigenvalue;
  std::cout << std::setprecision (7) << endl;
  std::cout << "******************** Solve Ax = B ********************" << endl;
  // Displaying the 2D vector
  std::cout << "The vector A is the following: " << endl;
  PrintVector2D(&A[0], ncols);
  std::cout << "The vector X is the following: " << endl;
  PrintVector(&X[0], ncols);
  std::cout << "The vector B is the following: " << endl;
  PrintVector(&B[0], ncols);
  // Create a new file named "C.dat"
  std::ofstream outputFile("C.dat");
  if (outputFile.is_open()) {
    // Write some text into the file
    outputFile << "%%MatrixMarket_Output_vector_C.dat matrix coordinate pattern general" << endl;
    outputFile << 1 << " " << ncols << " " << nrows;
    outputFile << endl;
    for (i = 0; i < ncols; i++)
      outputFile <<  1 << " " << i+1 << " " << X[i] << endl;
    // Close the file
    outputFile.close();
  } else {
    std::cout << "Error writing to file." << std::endl;
  }
  return 0;
}

#endif
