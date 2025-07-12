#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include "Matrix.h"
#include "Vector.h"
#include "MatrixMarket.h"

int main() {
    try {
        Matrix<double> A = MatrixMarket::readMatrix("A.dat");
        Vector<double> B = MatrixMarket::readVector("B.dat");

        if (A.getRows() != B.getSize()) {
            throw std::runtime_error("Matrix and vector dimensions are not compatible.");
        }

        Vector<double> X(B.getSize());

        // Perform matrix operations here

        std::cout << "Matrix A:\n" << A << std::endl;
        std::cout << "Vector B:\n" << B << std::endl;
        std::cout << "Vector X (result):\n" << X << std::endl;

        MatrixMarket::writeVector("X.dat", X);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
