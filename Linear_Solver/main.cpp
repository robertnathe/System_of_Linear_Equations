#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <cstddef> // for size_t
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdio>

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;

bool convergedUsingL1Norm(const vector<double>& X, const vector<double>& XPrevious, double tolerance, size_t size);
int Tridiagonal(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols);
int Tridiagonal_LU_Decomposition(double lowerDiagonal[], double mainDiagonal[], double upperDiagonal[], size_t num_cols);
int Tridiagonal_LU_Solve(double lowerDiagonal[], double mainDiagonal[], double upperDiagonal[], double B[], double x[], size_t num_cols);
int Jacobi_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, 
                  size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter);
int Jacobi_L2_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter);
int Gauss_Seidel_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter);
int SOR_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter, double omega);
int Write1DArrayToMatrixMarketFile(const double B[], size_t num_rows); 
int Write2DArrayToMatrixMarketFile(const double array_A[], size_t num_rows, size_t num_cols);
int WriteMatrixMarketMatrix(const matrix<double>& A, size_t num_rows, size_t num_cols);
int WriteMatrixMarketVector(const vector<double>& X, size_t num_cols);
int WriteTridiagonalToVectorMarketFile(double B[], vector<double>& X, size_t num_rows);
vector<double> CreateVectorFilledWithValue(size_t num_rows);
int ReadMatrixMarketMatrix(matrix<double>& A, size_t& num_rows, size_t& num_cols);
int ReadMatrixMarketVector(vector<double>& B, size_t& num_rows);
void PrintExecutionTime(const chrono::time_point<chrono::system_clock>& start,
                        const chrono::time_point<chrono::system_clock>& end);
void PrintVector(const string& label, const vector<double>& vec);
void PrintMatrix(const string& label, const matrix<double>& A);

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

int Tridiagonal_LU_Decomposition(double lowerDiagonal[], double mainDiagonal[], double upperDiagonal[], size_t num_cols) {
   for (size_t i = 0; i < num_cols-1; ++i) {
      if (mainDiagonal[i] == 0.0) return -1;
      lowerDiagonal[i] /= mainDiagonal[i];
      mainDiagonal[i+1] -= lowerDiagonal[i] * upperDiagonal[i];
   }
   if (mainDiagonal[num_cols-1] == 0.0) return -1;
   return 0;
}

int Tridiagonal_LU_Solve(double lowerDiagonal[], double mainDiagonal[],
                       double upperDiagonal[], double B[], double x[], size_t num_cols) {
   for (size_t i = 0; i < num_cols; ++i)
       if (mainDiagonal[i] == 0.0) return -1;
   x[0] = B[0];
   for (size_t i = 1; i < num_cols; ++i)
       x[i] = B[i] - lowerDiagonal[i-1] * x[i-1];
   x[num_cols-1] /= mainDiagonal[num_cols-1];
   for (size_t i = num_cols-1; i-- > 0;) {
      x[i] -= upperDiagonal[i] * x[i+1];
      x[i] /= mainDiagonal[i];
   }
   return 0;
}

int Jacobi_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, 
                  size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter) {
    if (num_rows != num_cols) {
        cerr << "Error: Jacobi method requires a square matrix." << endl;
        return -1;
    }
    if (A.size() != num_rows || (num_rows > 0 && A[0].size() != num_cols)) {
        cerr << "Error: Matrix dimensions mismatch." << endl;
        return -1;
    }
    if (B.size() != num_rows) {
        cerr << "Error: Vector B size does not match matrix dimensions." << endl;
        return -1;
    }
    vector<double> X_prev = X;
    double last_max_diff = 0.0;
    bool converged = false;
    for (size_t iter = 0; iter < Max_Iter; ++iter) {
        vector<double> X_new(num_rows, 0.0);
        double current_max_diff = 0.0;
        for (size_t i = 0; i < num_rows; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < num_cols; ++j) {
                if (i != j) sum += A[i][j] * X_prev[j];
            }
            X_new[i] = (B[i] - sum) / A[i][i];
            current_max_diff = max(current_max_diff, abs(X_new[i] - X_prev[i]));
        }
        last_max_diff = current_max_diff;
        X_prev = X_new;
        if (last_max_diff < tolerance) {
            converged = true;
            X = X_new;
            cout << "Converged in " << iter + 1 << " iterations (L∞ norm: " 
                 << last_max_diff << " < " << tolerance << ")" << endl;
            break;
        }
    }
    if (!converged) {
        X = X_prev;
        cerr << "Max iterations reached (L∞ norm: " 
             << last_max_diff << " > " << tolerance << ")" << endl;
        return 1;
    }
    return 0;
}

int Jacobi_L2_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, 
                 size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter) {
    if (num_rows != num_cols || num_rows == 0) {
        cerr << "Error: Matrix must be square." << endl;
        return -1;
    }
    vector<double> X_prev(num_rows, 0.0);
    size_t iter = 0;
    double residual;
    do {
        residual = 0.0;
        X_prev = X;
        for (size_t i = 0; i < num_rows; ++i) {
            double sum = 0.0;
            for (size_t j = 0; j < num_cols; ++j) {
                if (i != j) sum += A[i][j] * X_prev[j];
            }
            X[i] = (B[i] - sum) / A[i][i];
            residual += pow(X[i] - X_prev[i], 2);
        }
        residual = sqrt(residual);
        iter++;
    } while (residual > tolerance && iter < Max_Iter);
    if (residual <= tolerance) {
		cout << "Converged in " << iter << " iterations. (L2 norm: " 
                 << residual << " < " << tolerance << ")" << endl;
        return static_cast<int>(iter);
    } else {
        cerr << "Max iterations reached." << endl;
        return -1;
    }
}

int Gauss_Seidel_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter) {
    if (num_rows != num_cols) {
        cerr << "Error: Matrix must be square." << endl;
        return -1;
    }
    vector<double> inv_diag(num_rows);
    for (size_t i = 0; i < num_rows; ++i) {
        inv_diag[i] = 1.0 / A[i][i];
    }
    vector<double> X_prev(num_cols);
    size_t iter = 0;
    double max_diff;
    do {
        X_prev = X; // Save the current X to compute differences later
        max_diff = 0.0;
        for (size_t i = 0; i < num_rows; ++i) {
            double sum = 0.0;
            // Use updated X for j < i
            for (size_t j = 0; j < i; ++j) {
                sum += A[i][j] * X[j];
            }
            // Use X_prev for j > i
            for (size_t j = i + 1; j < num_cols; ++j) {
                sum += A[i][j] * X_prev[j];
            }
            X[i] = (B[i] - sum) * inv_diag[i];
            double diff = abs(X[i] - X_prev[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
        iter++;
        if (max_diff < tolerance) {
            cout << "Converged in " << iter << " iterations (L∞ norm: " 
                 << max_diff << " < " << tolerance << ")" << endl;
            return 0;
        }
    } while (iter < Max_Iter);
    cerr << "Max iterations reached (L∞ norm: " 
         << max_diff << " > " << tolerance << ")" << endl;
    return 1;
}

int SOR_L1_Norm(const matrix<double>& A, vector<double>& B, vector<double>& X, size_t num_rows, size_t num_cols, double tolerance, size_t Max_Iter, double omega) {
    if (num_rows != num_cols) {
        cerr << "Error: Matrix must be square." << endl;
        return -1;
    }
    vector<double> inv_diag(num_rows);
    for (size_t i = 0; i < num_rows; ++i) {
        inv_diag[i] = 1.0 / A[i][i];
    }
    vector<double> X_prev(num_cols);
    size_t iter = 0;
    double max_diff;
    do {
        X_prev = X; // Save the current X to compute differences later
        max_diff = 0.0;
        for (size_t i = 0; i < num_rows; ++i) {
            double sum = 0.0;
            // Use updated X for j < i
            for (size_t j = 0; j < i; ++j) {
                sum += A[i][j] * X[j];
            }
            // Use X_prev for j > i
            for (size_t j = i + 1; j < num_cols; ++j) {
                sum += A[i][j] * X_prev[j];
            }
            double gs_update = (B[i] - sum) * inv_diag[i];
            X[i] = (1.0 - omega) * X_prev[i] + omega * gs_update;
            double diff = abs(X[i] - X_prev[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }
        iter++;
        if (max_diff < tolerance) {
            cout << "Converged in " << iter << " iterations (L∞ norm: " 
                 << max_diff << " < " << tolerance << ")" << endl;
            return 0;
        }
    } while (iter < Max_Iter);
    cerr << "Max iterations reached (L∞ norm: " 
         << max_diff << " > " << tolerance << ")" << endl;
    return 1;
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
    //cout << "B:\n";
    for (size_t i = 0; i < num_rows; ++i) {
        outfile << 1 << " " << (i + 1) << " " << B[i] << "\n";
        //printf("%7.6lf    ", B[i]);
    }
    //cout << "\n";
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
        fprintf(fptr, "%zu %zu %.15g\n", 
                get<0>(entry), 
                get<1>(entry), 
                get<2>(entry));
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
    // Read header
    bool is_coordinate = false;
    while (getline(file, line) && line[0] == '%') {
        if (line.find("matrix coordinate") != string::npos) {
            is_coordinate = true;
        }
    }
    // Parse dimensions
    size_t rows = 0, cols = 1, entries = 0;
    istringstream dim_stream(line);
    if (is_coordinate) {
        // Coordinate format: rows cols entries
        size_t declared_rows, declared_cols;
        if (!(dim_stream >> declared_rows >> declared_cols >> entries)) {
            cerr << "Error parsing coordinate dimensions" << endl;
            return -1;
        }
        num_rows = (declared_rows == 1) ? declared_cols : declared_rows;
    } else {
        // Array format: rows cols
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
    // Read entries
    for (size_t i = 0; i < entries; ++i) {
        if (is_coordinate) {
            size_t row, col;
            double value;
            if (!(file >> row >> col >> value)) {
                cerr << "Error reading entry " << i << endl;
                return -1;
            }
            if (row == 1) {  // Row vector
                if (col > num_rows) {
                    cerr << "Column index out of bounds: " << col << endl;
                    return -1;
                }
                B[col-1] = value;
            } else {         // Column vector
                if (row > num_rows) {
                    cerr << "Row index out of bounds: " << row << endl;
                    return -1;
                }
                B[row-1] = value;
            }
        } else {
            // Array format - read values directly
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

void PrintExecutionTime(const chrono::time_point<chrono::system_clock>& start,
                        const chrono::time_point<chrono::system_clock>& end) {
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
        cerr << "F-ailed to read vector B." << endl;
        return -1;
    }
    if (num_rows != vector_rows) {
        cerr << "Error: Matrix and vector dimensions do not match." << endl;
        return -1;
    }
    
    // Convert matrix to flat array for Write2DArrayToMatrixMarketFile
    vector<double> flatA;
    flatA.reserve(num_rows * num_cols);
    for (const auto& row : A) {
        for (double val : row) {
            flatA.push_back(val);
        }
    }
    // Write matrix using flat array version
    Write2DArrayToMatrixMarketFile(flatA.data(), num_rows, num_cols);
    cout << endl;
    
    // PrintMatrix("A", A);
    PrintVector("B", B);
    chrono::time_point<chrono::system_clock> start, end;
    const int Max_Iter = 1000;
    const double tolerance = 0.00001;
    // --- Tridiagonal Solver ---
    vector<double> B_trid = B; // Copy for Tridiagonal
    vector<double> X_trid(num_cols, 0.0);
    start = chrono::system_clock::now();
    Tridiagonal(A, B_trid, X_trid, num_rows, num_cols);
    end = chrono::system_clock::now();
    cout << "Tridiagonal Method:" << endl;
    PrintExecutionTime(start, end);
    PrintVector("Solution", X_trid);
    // --- Jacobi L1 Norm Solver ---
    vector<double> B_jacobi_l1 = B; // Copy for Jacobi L1
    vector<double> X_jacobi_l1(num_cols, 0.0);
    cout << "\nJacobi (L1 Norm) Method:" << endl;
    start = chrono::system_clock::now();
    int status = Jacobi_L1_Norm(A, B_jacobi_l1, X_jacobi_l1, num_rows, num_cols, tolerance, Max_Iter);
    end = chrono::system_clock::now();
    PrintExecutionTime(start, end);
    if (status == 0) {
        PrintVector("Solution", X_jacobi_l1);
    } else {
        cerr << "Jacobi (L1 Norm) did not converge." << endl;
    }
    // --- Jacobi L2 Norm Solver ---
    vector<double> B_jacobi_l2 = B; // Copy for Jacobi L2
    vector<double> X_jacobi_l2(num_cols, 0.0);
    cout << "\nJacobi (L2 Norm) Method:" << endl;
    start = chrono::system_clock::now();
    status = Jacobi_L2_Norm(A, B_jacobi_l2, X_jacobi_l2, num_rows, num_cols, tolerance, Max_Iter);
    end = chrono::system_clock::now();
    PrintExecutionTime(start, end);
    if (status >= 0) {
        PrintVector("Solution", X_jacobi_l2);
    } else {
        cerr << "Jacobi (L2 Norm) did not converge." << endl;
    }
    
    //WriteMatrixMarketVector(X_jacobi_l2,num_cols);
    //WriteMatrixMarketMatrix(A,num_rows,num_cols);
    
    vector<double> B_gs_l1 = B; // Copy for Gauss_Seidel L1
    vector<double> X_gs_l1(num_cols, 0.0);
    cout << "\nGauss-Seidel (L1 Norm) Method:" << endl;
    start = std::chrono::system_clock::now();  
    status = Gauss_Seidel_L1_Norm(A, B_gs_l1, X_gs_l1, num_rows, num_cols, tolerance, Max_Iter);
    end = std::chrono::system_clock::now();
    PrintExecutionTime(start, end); 
    if (status >= 0) {
        PrintVector("Solution", X_gs_l1);
    } else {
        cerr << "Gauss_Seidel (L1 Norm) did not converge." << endl;
    }
    vector<double> B_sor_l1 = B; // Copy for SOR L1
    vector<double> X_sor_l1(num_cols, 0.0);
    cout << "\nSOR (L1 Norm) Method:" << endl;
    start = std::chrono::system_clock::now();  
    status = SOR_L1_Norm(A, B_sor_l1, X_sor_l1, num_rows, num_cols, tolerance, Max_Iter, omega);
    end = std::chrono::system_clock::now();
    PrintExecutionTime(start, end); 
    if (status >= 0) {
        PrintVector("Solution", X_sor_l1);
    } else {
        cerr << "SOR (L1 Norm) did not converge." << endl;
    }
    return 0;
}
