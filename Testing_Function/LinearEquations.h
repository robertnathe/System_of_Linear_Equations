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
#define A(i,j) A[(i)*nrows+(j)]

void PrintVector2D(const vector<double> A[], unsigned int nrows);
void PrintVector(const double TempVector[], unsigned int nrows);
int Output_Results(double X[], double B[], unsigned int nrows, unsigned int ncols);
int Output_Data(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols, double eigenvalue, double tolerance, unsigned int max_iter);

template<typename T>
using matrix = std::vector< std::vector<T> >;
#endif
