// main.cpp - CSR Optimized Sparse Matrix Operations
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

// ==================== Vector class ====================
template <typename T>
class Vector {
public:
    Vector(int size = 0) : data(size, 0) {}
    int getSize() const { return (int)data.size(); }
    T& operator[](int index) { return data[index]; }
    const T& operator[](int index) const { return data[index]; }
    void resize(int newSize) { data.assign(newSize, 0); }
private:
    std::vector<T> data;
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& vec) {
    for (int i = 0; i < vec.getSize(); ++i) os << vec[i] << '\n';
    return os;
}

// ==================== CSR Sparse Matrix class ====================
template <typename T>
class Matrix {
public:
    Matrix(int rows = 0, int cols = 0) : rows(rows), cols(cols) {
        row_ptr.assign(rows + 1, 0);
    }

    // CSR Construction: Converts temporary triplets to compressed format
    void buildFromTriplets(int r, int c, const std::vector<int>& rows_vec, 
                           const std::vector<int>& cols_vec, const std::vector<T>& vals_vec) {
        rows = r; cols = c;
        row_ptr.assign(rows + 1, 0);
        values = vals_vec;
        col_indices = cols_vec;

        // Calculate row offsets
        for (int i = 0; i < (int)rows_vec.size(); ++i) {
            row_ptr[rows_vec[i] + 1]++;
        }
        for (int i = 0; i < rows; ++i) {
            row_ptr[i + 1] += row_ptr[i];
        }
    }

    void scale(T factor) {
        for (auto& val : values) val *= factor;
    }

    // High-performance SpMV using CSR
    Vector<T> multiply(const Vector<T>& x) const {
        if (cols != x.getSize()) throw std::runtime_error("Dimension mismatch.");
        Vector<T> result(rows);
        for (int i = 0; i < rows; ++i) {
            T sum = 0;
            for (int j = row_ptr[i]; j < row_ptr[i+1]; ++j) {
                sum += values[j] * x[col_indices[j]];
            }
            result[i] = sum;
        }
        return result;
    }

    int getRows() const { return rows; }
    int getCols() const { return cols; }
    const std::vector<T>& getValues() const { return values; }
    const std::vector<int>& getColIndices() const { return col_indices; }
    const std::vector<int>& getRowPtr() const { return row_ptr; }

private:
    int rows, cols;
    std::vector<T> values;         // Non-zero values
    std::vector<int> col_indices;  // Column index of each value
    std::vector<int> row_ptr;      // Offset into values/col_indices for each row
};

// ==================== MatrixMarket utilities ====================
class MatrixMarket {
public:
    static Matrix<double> readMatrix(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) throw std::runtime_error("Could not open: " + filename);
        std::string line;
        while (std::getline(file, line) && line[0] == '%');
        int rows, cols, entries;
        std::stringstream(line) >> rows >> cols >> entries;
        
        std::vector<int> r_idx(entries), c_idx(entries);
        std::vector<double> vals(entries);
        for (int i = 0; i < entries; ++i) {
            file >> r_idx[i] >> c_idx[i] >> vals[i];
            r_idx[i]--; c_idx[i]--; // 0-based
        }
        Matrix<double> mat;
        mat.buildFromTriplets(rows, cols, r_idx, c_idx, vals);
        return mat;
    }

    static void writeMatrix(const std::string& filename, const Matrix<double>& mat) {
        std::ofstream file(filename);
        file << "%%MatrixMarket matrix coordinate real general\n";
        file << mat.getRows() << " " << mat.getCols() << " " << mat.getValues().size() << "\n";
        const auto& vals = mat.getValues();
        const auto& cols = mat.getColIndices();
        const auto& ptr = mat.getRowPtr();
        for (int i = 0; i < mat.getRows(); ++i) {
            for (int j = ptr[i]; j < ptr[i+1]; ++j) {
                file << (i + 1) << " " << (cols[j] + 1) << " " 
                     << std::fixed << std::setprecision(10) << vals[j] << "\n";
            }
        }
    }

    static Vector<double> readVector(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;
        while (std::getline(file, line) && line[0] == '%');
        int r, c, e;
        std::stringstream(line) >> r >> c >> e;
        Vector<double> vec(r == 1 ? c : r);
        for (int i = 0; i < e; ++i) {
            int vr, vc; double vv;
            file >> vr >> vc >> vv;
            vec[(r == 1) ? (vc - 1) : (vr - 1)] = vv;
        }
        return vec;
    }

    static void writeVector(const std::string& filename, const Vector<double>& vec) {
        std::ofstream file(filename);
        file << "%%MatrixMarket matrix coordinate real general\n";
        file << "1 " << vec.getSize() << " " << vec.getSize() << "\n";
        for (int i = 0; i < vec.getSize(); ++i) {
            file << "1 " << (i + 1) << " " << std::fixed << std::setprecision(10) << vec[i] << "\n";
        }
    }
};

int main() {
    try {
        Matrix<double> A = MatrixMarket::readMatrix("A.mtx");
        Vector<double> B = MatrixMarket::readVector("B.mtx");
        A.scale(2.0);
        MatrixMarket::writeMatrix("output.mtx", A);
        Vector<double> X = A.multiply(B);
        std::cout << "Successfully processed CSR Matrix. Result Vector X:\n" << X << std::endl;
        MatrixMarket::writeVector("X.mtx", X);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
