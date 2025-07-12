#include "MatrixMarket.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <tuple>

Matrix<double> MatrixMarket::readMatrix(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line) && line[0] == '%') {
        // Skip comments
    }

    int rows, cols, entries;
    std::stringstream ss(line);
    ss >> rows >> cols >> entries;

    Matrix<double> mat(rows, cols);

    for (int i = 0; i < entries; ++i) {
        int row, col;
        double value;
        file >> row >> col >> value;
        mat[row - 1][col - 1] = value;
    }

    return mat;
}

Vector<double> MatrixMarket::readVector(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line) && line[0] == '%') {
        // Skip comments
    }

    int rows, cols, entries;
    std::stringstream ss(line);
    ss >> rows >> cols >> entries;

    if (rows != 1 && cols != 1) {
        throw std::runtime_error("Matrix market file is not a vector: " + filename);
    }

    int size = (rows == 1) ? cols : rows;
    Vector<double> vec(size);

    for (int i = 0; i < entries; ++i) {
        int row, col;
        double value;
        file >> row >> col >> value;
        if (rows == 1) {
            vec[col - 1] = value;
        } else {
            vec[row - 1] = value;
        }
    }

    return vec;
}

void MatrixMarket::writeVector(const std::string& filename, const Vector<double>& vec) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    file << "%%MatrixMarket matrix array real general\n";
    file << vec.getSize() << " 1\n";
    for (int i = 0; i < vec.getSize(); ++i) {
        file << vec[i] << '\n';
    }
}
