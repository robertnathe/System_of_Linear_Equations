#ifndef ArithmeticRoutines_H
#define ArithmeticRoutines_H

#include <string.h>                                 // required for memcpy()
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cfloat>
#include <math.h>		// required for fabs() and for sqrt()

using namespace std;
#define a(i,j) a[(i)*nrows+(j)]

void PrintVector2D(const vector<double> TempMatrix[], unsigned int nrows);
void PrintVector(const double TempVector[], unsigned int nrows);
int Testing_Arithmetic_Matrix_Parameters(vector<double> A[], double X[], double B[], vector <double> AMat[], vector<double> BMat[], vector<double> CMat[], double MyVector1[], double MyVector2[], double MyVector3[], double Scalar, unsigned int nrows, unsigned int ncols);
int Add_Matrices(vector<double> CMat[], vector<double> AMat[], vector<double> BMat[], unsigned int nrows, unsigned int ncols);
int Subtract_Matrices(vector<double> CMat[], vector<double> AMat[], vector<double> BMat[], unsigned int nrows, unsigned int ncols);
extern int Add_Scalar_to_Diagonal(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols);
extern int Subtract_Scalar_from_Diagonal(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols);
void Multiply_Matrix_by_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols);
void Divide_Matrix_by_Scalar(vector<double> A[], double Scalar, unsigned int nrows, unsigned int ncols);
void Multiply_Matrices(vector<double> CMat[], vector <double> AMat[], unsigned int nrows, unsigned int ncols, vector <double> BMat[], unsigned int mcols);
int Multiply_Matrix_by_Vector(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols);
int array_Multiply_Matrix_by_Vector(double array_A[], double X[], double B[], unsigned int nrows, unsigned int ncols);
void Row_Transformation(vector<double> A[], double Scalar, unsigned int row1, unsigned int row2, unsigned int ncols);
void Column_Transformation(vector<double> A[], double Scalar, unsigned int col1, unsigned int col2, unsigned int nrows);
void Multiply_Row_by_Scalar(vector<double> A[], double Scalar, unsigned int row, unsigned int ncols);
void Multiply_Column_by_Scalar(vector<double> A[], double Scalar, unsigned int col, unsigned int nrows);
double Sum_over_Row(vector<double> A[], unsigned int ncols, unsigned int row);
double Sum_over_Rows(vector<double> A[], double X[], unsigned int nrows, unsigned int ncols);
double Sum_over_Column(vector<double> A[], unsigned int nrows, unsigned int col);
double Sum_over_Columns(vector<double> A[], double X[], unsigned int nrows, unsigned int ncols);
double Trace_of_Matrix(vector<double> A[], unsigned int nrows);
int Zero_Matrix(vector<double> A[], unsigned int nrows, unsigned int ncols);
int Identity_Matrix(vector<double> A[], unsigned int nrows);
void Identity_Matrix_ut(vector<double> A[], unsigned int n);
void Identity_Matrix_lt(vector<double> A[], unsigned int n);
double Bilinear_Function(double u[], vector<double> A[], double v[], unsigned int nrows, unsigned int ncols);
int Add_Vectors(double MyVector3[], double MyVector2[], double MyVector1[], unsigned int nrows); 
int Subtract_Vectors(double MyVector3[], double MyVector2[], double MyVector1[], unsigned int nrows);
int Multiply_Vector_by_Scalar(double MyVector1[], double Scalar, unsigned int nrows);
void Divide_Vector_by_Scalar(double X[], double TempScalar, unsigned int ncols);
int Multiply_Vector_by_Vector(vector<double> A[], unsigned int nrows, unsigned int ncols, double MyVector1[], double MyVector2[]);
int Zero_Vector(double MyVector1[], unsigned int nrows);
int Canonical_Basis_Vector(double MyVector1[], unsigned int i, unsigned int nrows);

int TransposeMatrix(vector<double> A[], unsigned int nrows, unsigned int ncols);
int TransposeSquareMatrix(vector<double> A[], unsigned int nrows, unsigned int ncols);

double Inner_Product(double u[], double v[], unsigned int nrows);
double Vector_Max_Norm(double TempVector[], unsigned int ncols);

double Vector_L2_Norm(double v[], unsigned int ncols);

#endif
