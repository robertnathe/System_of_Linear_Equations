#ifndef MAKE_GPP_MyIncludeFile
#define MAKE_GPP_MyIncludeFile
#include <string.h>                                 // required for memcpy()
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cfloat>
#include <math.h>		// required for fabs() and for sqrt()
using namespace std;
#include "LinearEquations.h"
//#include "LinearEquations.cpp"
#include "GeneralPurpose.h"
//#include "GeneralPurpose.cpp"	
#include "ArithmeticRoutines.h"
//#include "ArithmeticRoutines.cpp"
#define A(i,j) A[(i)*nrows+(j)]
template<typename T>
using matrix = std::vector< std::vector<T> >;
#endif
