#include "Matrix.h"

template <typename T>
Matrix<T>::Matrix(int rows, int cols) : rows(rows), cols(cols), data(rows, std::vector<T>(cols)) {}

template <typename T>
int Matrix<T>::getRows() const {
    return rows;
}

template <typename T>
int Matrix<T>::getCols() const {
    return cols;
}

template <typename T>
std::vector<T>& Matrix<T>::operator[](int index) {
    return data[index];
}

template <typename T>
const std::vector<T>& Matrix<T>::operator[](int index) const {
    return data[index];
}

template <typename T>
void Matrix<T>::resize(int newRows, int newCols) {
    rows = newRows;
    cols = newCols;
    data.resize(rows);
    for (int i = 0; i < rows; ++i) {
        data[i].resize(cols);
    }
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat) {
    for (int i = 0; i < mat.getRows(); ++i) {
        for (int j = 0; j < mat.getCols(); ++j) {
            os << std::setw(10) << mat[i][j];
        }
        os << '\n';
    }
    return os;
}
