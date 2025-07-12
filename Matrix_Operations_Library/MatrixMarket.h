#ifndef MATRIX_MARKET_H
#define MATRIX_MARKET_H

#include <string>
#include "Matrix.h"
#include "Vector.h"

class MatrixMarket {
public:
    static Matrix<double> readMatrix(const std::string& filename);
    static Vector<double> readVector(const std::string& filename);
    static void writeVector(const std::string& filename, const Vector<double>& vec);
};

#endif
