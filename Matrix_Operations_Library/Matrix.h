#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>

template <typename T>
class Matrix {
public:
    Matrix(int rows = 0, int cols = 0);

    int getRows() const;
    int getCols() const;

    std::vector<T>& operator[](int index);
    const std::vector<T>& operator[](int index) const;

    void resize(int rows, int cols);

private:
    int rows;
    int cols;
    std::vector<std::vector<T>> data;
};

#include "Matrix.tpp"

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat);

#endif // MATRIX_H
