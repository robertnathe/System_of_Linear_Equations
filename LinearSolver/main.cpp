#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstddef>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdio>

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;

// Function prototypes
bool convergedUsingL1Norm(const vector<double>& X, const vector<double>& XPrevious, double tolerance, size_t size);
int Tridiagonal(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols);
int Jacobi_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter);
int Jacobi_L2_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter);
int Gauss_Seidel_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter);
int SOR_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter, double omega);
int Write1DArrayToMatrixMarketFile(const double B[], size_t num_rows);
int Write2DArrayToMatrixMarketFile(const double array_A[], size_t num_rows, size_t num_cols);
int WriteMatrixMarketMatrix(const matrix<double>& A, size_t num_rows, size_t num_cols);
int WriteMatrixMarketVector(const vector<double>& X, size_t num_cols);
int WriteTridiagonalToVectorMarketFile(double B[], vector<double>& X, size_t num_rows);
int ReadMatrixMarketMatrix(matrix<double>& A, size_t& num_rows, size_t& num_cols);
int ReadMatrixMarketVector(vector<double>& B, size_t& num_rows);
void PrintExecutionTime(const chrono::time_point<chrono::system_clock>& start, const chrono::time_point<chrono::system_clock>& end);
void PrintVector(const string& label, const vector<double>& vec);
void PrintMatrix(const string& label, const matrix<double>& A);
bool isTridiagonal(const matrix<double>& A, size_t num_rows, size_t num_cols);

// New functions for residual computation
vector<double> computeResidual(const matrix<double>& A, const vector<double>& X, const vector<double>& B, size_t num_rows, size_t num_cols) {
    vector<double> residual(num_rows, 0.0);
    for (size_t i = 0; i < num_rows; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < num_cols; ++j) {
            sum += A[i][j] * X[j];
        }
        residual[i] = sum - B[i];
    }
    return residual;
}

double computeL2Norm(const vector<double>& vec) {
    double sum = 0.0;
    for (double val : vec) {
        sum += val * val;
    }
    return sqrt(sum);
}

bool isTridiagonal(const matrix<double>& A, size_t num_rows, size_t num_cols) {
    if (num_rows != num_cols) return false;
    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_cols; ++j) {
            if ((i > 0 && j < i - 1) || j > i + 1) {
                if (A[i][j] != 0.0) return false;
            }
        }
    }
    return true;
}

bool convergedUsingL1Norm(const vector<double>& X, const vector<double>& XPrevious, double tolerance, size_t size) {
    double distance_max = 0.0;
    for (size_t i = 0; i < size; ++i) {
        distance_max = max(distance_max, abs(X[i] - XPrevious[i]));
    }
    return distance_max < tolerance;
}

int Tridiagonal(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols) {
    vector<double> lowerDiagonal(num_cols);
    vector<double> mainDiagonal(num_cols);
    vector<double> upperDiagonal(num_cols);
    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_cols; ++j) {
            if (i == j) {
                mainDiagonal[i] = A[i][j];
            } else if (i == j - 1) {
                upperDiagonal[i] = A[i][j];
            } else if (i == j + 1) {
                lowerDiagonal[i] = A[i][j];
            }
        }
    }
    for (size_t i = 1; i < num_cols; ++i) {
        double factor = lowerDiagonal[i] / mainDiagonal[i - 1];
        mainDiagonal[i] -= factor * upperDiagonal[i - 1];
        B[i] -= factor * B[i - 1];
    }
    X[num_cols - 1] = B[num_cols - 1] / mainDiagonal[num_cols - 1];
    for (size_t i = num_cols - 1; i-- > 0;) {
        X[i] = (B[i] - upperDiagonal[i] * X[i + 1]) / mainDiagonal[i];
    }
    WriteTridiagonalToVectorMarketFile(B.data(), X, num_rows);
    return 0;
}

// Updated Jacobi_L1_Norm using residual norm (renamed for clarity)
int Jacobi_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter) {
    if (num_rows != num_cols) {
        cerr << "Error: Jacobi method requires a square matrix." << endl;
        return -1;
    }
    vector<double> X_prev(num_rows, 0.0);
    size_t iter = 0;
    double B_norm = computeL2Norm(B);
    double rel_tolerance = tolerance * B_norm;
    double residual_norm;
    do {
        X_prev = X;
        for (size_t i = 0; i < num_rows; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < num_cols; ++j) {
                if (i != j) sum += A[i][j] * X_prev[j];
            }
            X[i] = (B[i] - sum) / A[i][i];
        }
        vector<double> residual = computeResidual(A, X, B, num_rows, num_cols);
        residual_norm = computeL2Norm(residual);
        iter++;
    } while (residual_norm > rel_tolerance && iter < Max_Iter);
    if (residual_norm <= rel_tolerance) {
        cout << "Converged in " << iter << " iterations. (Relative L2 residual: " << residual_norm / B_norm << " < " << tolerance << ")" << endl;
        return static_cast<int>(iter);
    } else {
        cerr << "Max iterations reached. (Relative L2 residual: " << residual_norm / B_norm << " > " << tolerance << ")" << endl;
        return -1;
    }
}

// Updated Jacobi_L2_Norm using residual norm
int Jacobi_L2_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter) {
    if (num_rows != num_cols || num_rows == 0) {
        cerr << "Error: Matrix must be square." << endl;
        return -1;
    }
    vector<double> X_prev(num_rows, 0.0);
    size_t iter = 0;
    double B_norm = computeL2Norm(B);
    double rel_tolerance = tolerance * B_norm;
    double residual_norm;
    do {
        X_prev = X;
        for (size_t i = 0; i < num_rows; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < num_cols; ++j) {
                if (i != j) sum += A[i][j] * X_prev[j];
            }
            X[i] = (B[i] - sum) / A[i][i];
        }
        vector<double> residual = computeResidual(A, X, B, num_rows, num_cols);
        residual_norm = computeL2Norm(residual);
        iter++;
    } while (residual_norm > rel_tolerance && iter < Max_Iter);
    if (residual_norm <= rel_tolerance) {
        cout << "Converged in " << iter << " iterations. (Relative L2 residual: " << residual_norm / B_norm << " < " << tolerance << ")" << endl;
        return static_cast<int>(iter);
    } else {
        cerr << "Max iterations reached. (Relative L2 residual: " << residual_norm / B_norm << " > " << tolerance << ")" << endl;
        return -1;
    }
}

// Updated Gauss_Seidel_L1_Norm using residual norm
int Gauss_Seidel_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter) {
    if (num_rows != num_cols) {
        cerr << "Error: Matrix must be square." << endl;
        return -1;
    }
    vector<double> X_prev(num_rows, 0.0);
    size_t iter = 0;
    double B_norm = computeL2Norm(B);
    double rel_tolerance = tolerance * B_norm;
    double residual_norm;
    do {
        X_prev = X;
        for (size_t i = 0; i < num_rows; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < i; ++j) {
                sum += A[i][j] * X[j];
            }
            for (size_t j = i + 1; j < num_cols; ++j) {
                sum += A[i][j] * X_prev[j];
            }
            X[i] = (B[i] - sum) / A[i][i];
        }
        vector<double> residual = computeResidual(A, X, B, num_rows, num_cols);
        residual_norm = computeL2Norm(residual);
        iter++;
    } while (residual_norm > rel_tolerance && iter < Max_Iter);
    if (residual_norm <= rel_tolerance) {
        cout << "Converged in " << iter << " iterations. (Relative L2 residual: " << residual_norm / B_norm << " < " << tolerance << ")" << endl;
        return static_cast<int>(iter);
    } else {
        cerr << "Max iterations reached. (Relative L2 residual: " << residual_norm / B_norm << " > " << tolerance << ")" << endl;
        return -1;
    }
}

// Updated SOR_L1_Norm using residual norm
int SOR_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter, double omega) {
    if (num_rows != num_cols) {
        cerr << "Error: Matrix must be square." << endl;
        return -1;
    }
    vector<double> X_prev(num_rows, 0.0);
    size_t iter = 0;
    double B_norm = computeL2Norm(B);
    double rel_tolerance = tolerance * B_norm;
    double residual_norm;
    do {
        X_prev = X;
        for (size_t i = 0; i < num_rows; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < i; ++j) {
                sum += A[i][j] * X[j];
            }
            for (size_t j = i + 1; j < num_cols; ++j) {
                sum += A[i][j] * X_prev[j];
            }
            double gs_update = (B[i] - sum) / A[i][i];
            X[i] = (1.0 - omega) * X_prev[i] + omega * gs_update;
        }
        vector<double> residual = computeResidual(A, X, B, num_rows, num_cols);
        residual_norm = computeL2Norm(residual);
        iter++;
    } while (residual_norm > rel_tolerance && iter < Max_Iter);
    if (residual_norm <= rel_tolerance) {
        cout << "Converged in " << iter << " iterations. (Relative L2 residual: " << residual_norm / B_norm << " < " << tolerance << ")" << endl;
        return static_cast<int>(iter);
    } else {
        cerr << "Max iterations reached. (Relative L2 residual: " << residual_norm / B_norm << " > " << tolerance << ")" << endl;
        return -1;
    }
}

int Write1DArrayToMatrixMarketFile(const double B[], size_t num_rows) {
    ofstream outfile("B_out.dat");
    if (!outfile.is_open()) {
        cerr << "Error opening file: B_new_out.dat" << endl;
        return 1;
    }
    outfile << "%%MatrixMarket matrix coordinate real general\n";
    outfile << 1 << " " << num_rows << " " << num_rows << "\n";
    outfile << fixed << setprecision(6);
    for (size_t i = 0; i < num_rows; ++i) {
        outfile << 1 << " " << (i + 1) << " " << B[i] << "\n";
    }
    outfile.close();
    return 0;
}

int Write2DArrayToMatrixMarketFile(const double array_A[], size_t num_rows, size_t num_cols) {
    ofstream outfile("A_out.dat");
    if (!outfile.is_open()) {
        cerr << "Error opening file for writing: A_out.dat" << endl;
        return 1;
    }
    outfile << "%%%%MatrixMarket matrix coordinate real general\n"
            << num_rows << " " << num_cols << " " << num_rows * num_cols << "\n";
    ostringstream file_buffer;
    file_buffer << scientific << setprecision(6);
    cout << "A:\n";
    for (size_t i = 0; i < num_rows; ++i) {
        ostringstream console_line;
        console_line << fixed << setprecision(5);
        for (size_t j = 0; j < num_cols; ++j) {
            const double val = array_A[i * num_cols + j];
            file_buffer << i+1 << " " << j+1 << " " << val << "\n";
            console_line << setw(6) << val;
            if (j < num_cols-1) console_line << "    ";
        }
        console_line << "\n";
        cout << console_line.str();
    }
    outfile << file_buffer.str();
    outfile.close();
    return 0;
}

int WriteMatrixMarketMatrix(const matrix<double>& A, size_t num_rows, size_t num_cols) {
    FILE* fptr = fopen("A_out.dat", "w");
    if (!fptr) return -1;
    vector<tuple<size_t, size_t, double>> entries;
    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_cols; ++j) {
            if (A[i][j] != 0.0) {
                entries.emplace_back(i+1, j+1, A[i][j]);
            }
        }
    }
    fprintf(fptr, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fptr, "%zu %zu %zu\n", num_rows, num_cols, entries.size());
    for (const auto& entry : entries) {
        fprintf(fptr, "%zu %zu %.15g\n", get<0>(entry), get<1>(entry), get<2>(entry));
    }
    fclose(fptr);
    return 0;
}

int WriteMatrixMarketVector(const vector<double>& X, size_t num_cols) {
    FILE* fptr = fopen("X_out.dat", "w");
    if (!fptr) return -1;
    fprintf(fptr, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fptr, "1 %zu %zu\n", num_cols, X.size());
    for (size_t i = 0; i < X.size(); ++i) {
        fprintf(fptr, "1 %zu %.15g\n", i + 1, X[i]);
    }
    fclose(fptr);
    return 0;
}

int WriteTridiagonalToVectorMarketFile(double B[], vector<double>& X, size_t num_rows) {
    Write1DArrayToMatrixMarketFile(B, num_rows);
    WriteMatrixMarketVector(X, X.size());
    return 0;
}

int ReadMatrixMarketMatrix(matrix<double>& A, size_t& num_rows, size_t& num_cols) {
    ifstream file("A.dat");
    if (!file.is_open()) return -1;
    string line;
    while (getline(file, line) && line[0] == '%');
    size_t num_entries;
    istringstream header_stream(line);
    if (!(header_stream >> num_rows >> num_cols >> num_entries)) return -1;
    A.resize(num_rows, vector<double>(num_cols, 0.0));
    for (size_t i = 0; i < num_entries; ++i) {
        size_t row, col;
        double value;
        if (!(file >> row >> col >> value)) return -1;
        row--; col--;
        A[row][col] = value;
    }
    return 0;
}

int ReadMatrixMarketVector(vector<double>& B, size_t& num_rows) {
    ifstream file("B.dat");
    if (!file.is_open()) {
        cerr << "Error opening vector file B.dat" << endl;
        return -1;
    }
    string line;
    bool is_coordinate = false;
    while (getline(file, line) && line[0] == '%') {
        if (line.find("matrix coordinate") != string::npos) {
            is_coordinate = true;
        }
    }
    size_t rows = 0, cols = 1, entries = 0;
    istringstream dim_stream(line);
    if (is_coordinate) {
        size_t declared_rows, declared_cols;
        if (!(dim_stream >> declared_rows >> declared_cols >> entries)) {
            cerr << "Error parsing coordinate dimensions" << endl;
            return -1;
        }
        num_rows = (declared_rows == 1) ? declared_cols : declared_rows;
    } else {
        if (!(dim_stream >> rows >> cols)) {
            cerr << "Error parsing array dimensions" << endl;
            return -1;
        }
        num_rows = rows;
        entries = rows * cols;
    }
    if (cols != 1) {
        cerr << "Error: Vector must have 1 column (found " << cols << ")" << endl;
        return -1;
    }
    B.resize(num_rows, 0.0);
    for (size_t i = 0; i < entries; ++i) {
        if (is_coordinate) {
            size_t row, col;
            double value;
            if (!(file >> row >> col >> value)) {
                cerr << "Error reading entry " << i << endl;
                return -1;
            }
            if (row == 1) {
                if (col > num_rows) {
                    cerr << "Column index out of bounds: " << col << endl;
                    return -1;
                }
                B[col-1] = value;
            } else {
                if (row > num_rows) {
                    cerr << "Row index out of bounds: " << row << endl;
                    return -1;
                }
                B[row-1] = value;
            }
        } else {
            double value;
            if (!(file >> value)) {
                cerr << "Error reading entry " << i << endl;
                return -1;
            }
            if (i >= num_rows) {
                cerr << "Too many entries for vector size " << num_rows << endl;
                return -1;
            }
            B[i] = value;
        }
    }
    return 0;
}

void PrintExecutionTime(const chrono::time_point<chrono::system_clock>& start, const chrono::time_point<chrono::system_clock>& end) {
    chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time: " << elapsed.count() << "s\n";
}

void PrintVector(const string& label, const vector<double>& vec) {
    cout << label << ":\n[";
    for (size_t i = 0; i < vec.size(); ++i) {
        cout << fixed << setprecision(5) << vec[i];
        if (i != vec.size()-1) cout << ", ";
    }
    cout << "]\n" << endl;
}

void PrintMatrix(const string& label, const matrix<double>& A) {
    cout << label << ":\n";
    for (const auto& row : A) {
        for (double value : row) {
            cout << fixed << setprecision(5) << setw(10) << value;
        }
        cout << endl;
    }
    cout << endl;
}

int main() {
    size_t num_rows = 0, num_cols = 0;
    double omega {1.25};
    matrix<double> A;
    vector<double> B;
    if (ReadMatrixMarketMatrix(A, num_rows, num_cols) != 0) {
        cerr << "Failed to read matrix A." << endl;
        return -1;
    }
    size_t vector_rows = 0;
    if (ReadMatrixMarketVector(B, vector_rows) != 0) {
        cerr << "Failed to read vector B." << endl;
        return -1;
    }
    if (num_rows != vector_rows) {
        cerr << "Error: Matrix and vector dimensions do not match." << endl;
        return -1;
    }

    vector<double> flatA;
    flatA.reserve(num_rows * num_cols);
    for (const auto& row : A) {
        for (double val : row) {
            flatA.push_back(val);
        }
    }
    Write2DArrayToMatrixMarketFile(flatA.data(), num_rows, num_cols);
    cout << endl;

    PrintVector("B", B);
    chrono::time_point<chrono::system_clock> start, end;
    const int Max_Iter = 1000;
    const double tolerance = 0.00001;

    // Tridiagonal Solver with residual check
    if (isTridiagonal(A, num_rows, num_cols)) {
        vector<double> B_trid = B;
        vector<double> X_trid(num_cols, 0.0);
        start = chrono::system_clock::now();
        Tridiagonal(A, B_trid, X_trid, num_rows, num_cols);
        end = chrono::system_clock::now();
        cout << "Tridiagonal Method:" << endl;
        PrintExecutionTime(start, end);
        PrintVector("Solution", X_trid);
        vector<double> residual = computeResidual(A, X_trid, B, num_rows, num_cols);
        double residual_norm = computeL2Norm(residual);
        cout << "Residual L2 Norm: " << residual_norm << endl << endl;
    } else {
        cout << "Matrix is not tridiagonal. Skipping Tridiagonal solver.\n" << endl;
    }

    // Jacobi L1 Norm Solver (now using residual)
    vector<double> B_jacobi_l1 = B;
    vector<double> X_jacobi_l1(num_cols, 0.0);
    cout << "Jacobi (Residual Norm) Method:" << endl;
    start = chrono::system_clock::now();
    int status = Jacobi_L1_Norm(A, B_jacobi_l1, X_jacobi_l1, num_rows, num_cols, tolerance, Max_Iter);
    end = chrono::system_clock::now();
    PrintExecutionTime(start, end);
    if (status >= 0) {
        PrintVector("Solution", X_jacobi_l1);
    } else {
        cerr << "Jacobi (Residual Norm) did not converge." << endl;
    }

    // Gauss-Seidel Solver (now using residual)
    vector<double> B_gs_l1 = B;
    vector<double> X_gs_l1(num_cols, 0.0);
    cout << "\nGauss-Seidel (Residual Norm) Method:" << endl;
    start = chrono::system_clock::now();
    status = Gauss_Seidel_L1_Norm(A, B_gs_l1, X_gs_l1, num_rows, num_cols, tolerance, Max_Iter);
    end = chrono::system_clock::now();
    PrintExecutionTime(start, end);
    if (status >= 0) {
        PrintVector("Solution", X_gs_l1);
    } else {
        cerr << "Gauss-Seidel (Residual Norm) did not converge." << endl;
    }

    // SOR Solver (now using residual)
    vector<double> B_sor_l1 = B;
    vector<double> X_sor_l1(num_cols, 0.0);
    cout << "\nSOR (Residual Norm) Method:" << endl;
    start = chrono::system_clock::now();
    status = SOR_L1_Norm(A, B_sor_l1, X_sor_l1, num_rows, num_cols, tolerance, Max_Iter, omega);
    end = chrono::system_clock::now();
    PrintExecutionTime(start, end);
    if (status >= 0) {
        PrintVector("Solution", X_sor_l1);
    } else {
        cerr << "SOR (Residual Norm) did not converge." << endl;
    }

    return 0;
}
