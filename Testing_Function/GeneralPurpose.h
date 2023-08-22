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
#define a(i,j) a[(i)*nrows+(j)]
void PrintVector2D(const vector<double> A[], unsigned int nrows);
void PrintVector(const double TempVector[], unsigned int nrows);
int Testing_General_Purpose_Parameters(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols);
void Copy_Matrix(vector<double> A[], vector<double> B[], unsigned int nrows, unsigned int ncols);
void Copy_Vector(double *d, double *s, unsigned int ncols);
void Get_Row(double X[], vector <double> A[], unsigned int kth_row, unsigned int nrows, unsigned int ncols);
void Get_Column(double X[], vector <double> A[], unsigned int kth_column, unsigned int nrows, unsigned int ncols);
void Set_Row(double X[], vector <double> A[], unsigned int kth_row, unsigned int nrows, unsigned int ncols);
void Set_Column(double X[], vector <double> A[], unsigned int kth_column, unsigned int nrows, unsigned int ncols);
void Get_Diagonal(double X[], vector <double> A[], unsigned int nrows, unsigned int ncols);
void Set_Diagonal(double X[], vector <double> A[], unsigned int nrows, unsigned int ncols);
void Set_Diagonal_to_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols);
void Get_Submatrix(vector<double> A[], unsigned int nrows, unsigned int ncols, vector<double> S[], unsigned int mrows, unsigned int mcols, unsigned int row, unsigned int col);
void Set_Submatrix(vector<double> A[], unsigned int nrows, unsigned int ncols, vector<double> S[], unsigned int mrows, unsigned int mcols, unsigned int row, unsigned int col, double Scalar);
void Fill_Matrix_with_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols);
int Join_Rows(vector<double> A[], vector<double> BMatrix[], vector <double> CMatrix[],  unsigned int nrows, unsigned int ncols, unsigned int mcols);
int Join_Columns(vector<double> A[], vector<double> BMatrix[], vector <double> CMatrix[],  unsigned int nrows, unsigned int ncols, unsigned int mcols);
void Interchange_Rows(vector<double> A[], unsigned int row1, unsigned int row2, unsigned int ncols);
void Interchange_Columns(vector<double> A[], unsigned int col1, unsigned int col2, unsigned int nrows, unsigned int ncols);
int Vec_Swap(double MyVector1[], double MyVector2[], unsigned int Vec_Size);
int Test_Vec_Swap(double MyVector1[], double MyVector2[], unsigned int Vec_Size);
#endif
