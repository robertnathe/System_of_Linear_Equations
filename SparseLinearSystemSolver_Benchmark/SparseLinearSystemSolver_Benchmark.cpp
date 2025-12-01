#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <string>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <cstring>

// A small tolerance for floating-point comparisons.
static constexpr double PIVOT_TOLERANCE = 1e-12;
static constexpr std::size_t MAX_PRINT_SIZE = 10;

// Compressed Sparse Row (CSR) matrix format for efficient storage and computation.
struct SparseMatrix {
    std::vector<double> values;
    std::vector<std::size_t> col_indices;
    std::vector<std::size_t> row_ptr;
    std::vector<std::size_t> diag_indices;
    std::size_t num_rows = 0;
    std::size_t num_cols = 0;

    SparseMatrix() = default;
    SparseMatrix(std::size_t rows, std::size_t cols)
        : num_rows(rows), num_cols(cols) {}
};

// Represents a matrix entry for sorting and building the CSR matrix.
struct Triplet {
    std::size_t row, col;
    double value;

    bool operator<(const Triplet& other) const {
        if (row != other.row) return row < other.row;
        return col < other.col;
    }
};

// Reads a matrix from a MatrixMarket file into a sparse CSR format.
[[nodiscard]] int ReadMatrixMarketMatrix(SparseMatrix& A, bool& is_tridiagonal) {
    FILE* file = fopen("A.mtx", "r");
    if (!file) {
        perror("Error opening matrix file A.mtx");
        return -1;
    }

    char line[1024];
    do {
        if (fgets(line, sizeof(line), file) == nullptr) {
            fprintf(stderr, "Error reading matrix file header or empty file.\n");
            fclose(file);
            return -1;
        }
    } while (line[0] == '%');

    std::size_t num_rows, num_cols, num_entries;
    if (sscanf(line, "%zu %zu %zu", &num_rows, &num_cols, &num_entries) != 3) {
        fprintf(stderr, "Error parsing matrix header.\n");
        fclose(file);
        return -1;
    }

    std::vector<Triplet> triplets;
    triplets.reserve(num_entries);
    is_tridiagonal = (num_rows == num_cols && num_rows > 0);

    for (std::size_t i = 0; i < num_entries; ++i) {
        std::size_t r, c;
        double v;
        if (fscanf(file, "%zu %zu %lf", &r, &c, &v) != 3) {
            fprintf(stderr, "Error reading matrix entry #%zu.\n", i + 1);
            fclose(file);
            return -1;
        }
        if (r == 0 || r > num_rows || c == 0 || c > num_cols) {
            fprintf(stderr, "Matrix index out of bounds.\n");
            fclose(file);
            return -1;
        }
        triplets.push_back({r - 1, c - 1, v});
        if (is_tridiagonal && std::abs(static_cast<ptrdiff_t>(r) - static_cast<ptrdiff_t>(c)) > 1) {
            if (std::abs(v) > PIVOT_TOLERANCE) {
                is_tridiagonal = false;
            }
        }
    }
    fclose(file);
    
    std::sort(triplets.begin(), triplets.end());

    if (!triplets.empty()) {
        std::size_t write_idx = 0;
        for (std::size_t read_idx = 1; read_idx < triplets.size(); ++read_idx) {
            if (triplets[read_idx].row == triplets[write_idx].row && triplets[read_idx].col == triplets[write_idx].col) {
                triplets[write_idx].value += triplets[read_idx].value;
            } else {
                ++write_idx;
                if (write_idx != read_idx) {
                    triplets[write_idx] = triplets[read_idx];
                }
            }
        }
        triplets.resize(write_idx + 1);
    }

    A = SparseMatrix(num_rows, num_cols);
    A.values.reserve(triplets.size());
    A.col_indices.reserve(triplets.size());
    A.row_ptr.reserve(num_rows + 1);
    A.diag_indices.assign(num_rows, std::numeric_limits<std::size_t>::max());

    A.row_ptr.push_back(0);
    std::size_t current_row = 0;
    for (const auto& t : triplets) {
        if (std::abs(t.value) > PIVOT_TOLERANCE) {
            while (current_row < t.row) {
                A.row_ptr.push_back(A.values.size());
                current_row++;
            }
            if (t.row == t.col) {
                A.diag_indices[t.row] = A.values.size();
            }
            A.values.push_back(t.value);
            A.col_indices.push_back(t.col);
        }
    }

    while (current_row < A.num_rows) {
        A.row_ptr.push_back(A.values.size());
        current_row++;
    }
    A.values.shrink_to_fit();
    A.col_indices.shrink_to_fit();
    return 0;
}

[[nodiscard]] int ReadMatrixMarketVector(std::vector<double>& B) {
    FILE* file = fopen("B.mtx", "r");
    if (!file) {
        perror("Error opening vector file B.mtx");
        return -1;
    }

    char line[1024];
    do {
        if (fgets(line, sizeof(line), file) == nullptr) {
            fprintf(stderr, "Error reading vector file header or empty file.\n");
            fclose(file);
            return -1;
        }
    } while (line[0] == '%');

    std::size_t rows, cols, entries;
    if (sscanf(line, "%zu %zu %zu", &rows, &cols, &entries) != 3) {
        fprintf(stderr, "Error parsing vector dimensions.\n");
        fclose(file);
        return -1;
    }
    if (cols != 1 && rows != 1) {
        fprintf(stderr, "Error: B.mtx must represent a vector (1 column or 1 row).\n");
        fclose(file);
        return -1;
    }

    const bool is_column_vector = (cols == 1);
    const std::size_t num_elements = is_column_vector ? rows : cols;
    B.assign(num_elements, 0.0);

    for (std::size_t i = 0; i < entries; ++i) {
        std::size_t r, c;
        double value;
        if (fscanf(file, "%zu %zu %lf", &r, &c, &value) != 3) {
            fprintf(stderr, "Error reading vector entry #%zu.\n", i + 1);
            fclose(file);
            return -1;
        }
        const std::size_t index = is_column_vector ? (r - 1) : (c - 1);
        if (index < num_elements) {
            B[index] += value;
        } else {
            fprintf(stderr, "Vector index out of bounds.\n");
            fclose(file);
            return -1;
        }
    }
    fclose(file);
    return 0;
}

[[nodiscard]] int WriteMatrixMarketVector(const char* filename, const std::vector<double>& vec) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error opening file for writing: " << filename << '\n';
        return -1;
    }
    // MatrixMarket header for a dense vector (array format)
    file << "%%MatrixMarket matrix array real general\n";
    file << "1 " << vec.size() << " "<< vec.size() << '\n';
    file << std::fixed << std::setprecision(15);
    int i = 1;
    for (const auto& val : vec) {
        file << 1 << " " << i << " " << val << '\n';
        i++;
    }
    return 0;
}

[[nodiscard]] int TridiagonalSolve(
    const std::size_t n,
    std::vector<double> lower_diag,
    std::vector<double> main_diag,
    std::vector<double> upper_diag,
    std::vector<double> rhs,
    std::vector<double>& X)
{
    if (n == 0) return 0;
    X.assign(n, 0.0);
    if (n == 1) {
        if (std::abs(main_diag[0]) < PIVOT_TOLERANCE) {
            std::cerr << "Error: Tridiagonal algorithm failed due to zero pivot.\n";
            return -1;
        }
        X[0] = rhs[0] / main_diag[0];
        return 0;
    }

    for (std::size_t i = 1; i < n; ++i) {
        if (std::abs(main_diag[i - 1]) < PIVOT_TOLERANCE) {
             std::cerr << "Error: Tridiagonal algorithm failed due to zero pivot.\n";
             return -1;
        }
        const double factor = lower_diag[i - 1] / main_diag[i - 1];
        main_diag[i] -= factor * upper_diag[i - 1];
        rhs[i] -= factor * rhs[i - 1];
    }

    if (std::abs(main_diag[n - 1]) < PIVOT_TOLERANCE) {
          std::cerr << "Error: Tridiagonal algorithm failed due to zero pivot.\n";
          return -1;
    }
    X[n - 1] = rhs[n - 1] / main_diag[n - 1];
    for (ptrdiff_t i = n - 2; i >= 0; --i) {
        X[i] = (rhs[i] - upper_diag[i] * X[i + 1]) / main_diag[i];
    }
    return 0;
}

[[nodiscard]] int Jacobi(const SparseMatrix& A, const std::vector<double>& B, std::vector<double>& X, double tolerance, std::size_t Max_Iter) {
    const std::size_t n = A.num_rows;
    if (n != A.num_cols) {
        std::cerr << "Error: Jacobi method requires a square matrix.\n";
        return -1;
    }
    std::vector<double> D_inv(n);
    for(std::size_t i = 0; i < n; ++i) {
        const std::size_t diag_idx = A.diag_indices[i];
        if (diag_idx == std::numeric_limits<std::size_t>::max()) {
            std::cerr << "Error: Jacobi method requires a diagonal element for each row (" << i << ").\n";
            return -1;
        }
        const double diag_val = A.values[diag_idx];
        if (std::abs(diag_val) < PIVOT_TOLERANCE) {
            std::cerr << "Error: Jacobi method requires non-zero diagonal elements (row " << i << ").\n";
            return -1;
        }
        D_inv[i] = 1.0 / diag_val;
    }

    std::vector<double> X_next(n);
    for (std::size_t iter = 0; iter < Max_Iter; ++iter) {
        double max_diff = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double off_diag_sum = 0.0;
            const std::size_t row_start = A.row_ptr[i];
            const std::size_t row_end = A.row_ptr[i+1];
            const std::size_t diag_idx = A.diag_indices[i];
            for (std::size_t j_idx = row_start; j_idx < row_end; ++j_idx) {
                if (j_idx != diag_idx) {
                    off_diag_sum += A.values[j_idx] * X[A.col_indices[j_idx]];
                }
            }
            X_next[i] = (B[i] - off_diag_sum) * D_inv[i];
            max_diff = std::max(max_diff, std::abs(X_next[i] - X[i]));
        }
        X.swap(X_next);
        if (max_diff < tolerance) {
            std::cout << "Converged in " << iter + 1 << " iterations (L-infinity norm of change < " << tolerance << ")\n";
            return static_cast<int>(iter + 1);
        }
    }
    std::cerr << "Max iterations reached without convergence.\n";
    return -1;
}

[[nodiscard]] int SOR(const SparseMatrix& A, const std::vector<double>& B, std::vector<double>& X, double tolerance, std::size_t Max_Iter, double omega) {
    const std::size_t n = A.num_rows;
    if (n != A.num_cols) {
        std::cerr << "Error: SOR requires a square matrix.\n";
        return -1;
    }
    if (omega <= 0.0 || omega >= 2.0) {
        std::cerr << "Warning: SOR omega is typically in (0, 2) for convergence. Value used: " << omega << "\n";
    }
    
    std::vector<double> D_inv(n);
    for(std::size_t i = 0; i < n; ++i) {
        const std::size_t diag_idx = A.diag_indices[i];
        if (diag_idx == std::numeric_limits<std::size_t>::max()) {
            std::cerr << "Error: SOR requires a diagonal element for each row (" << i << ").\n";
            return -1;
        }
        const double diag_val = A.values[diag_idx];
        if (std::abs(diag_val) < PIVOT_TOLERANCE) {
            std::cerr << "Error: SOR requires non-zero diagonal elements (row " << i << ").\n";
            return -1;
        }
        D_inv[i] = 1.0 / diag_val;
    }
    
    for (std::size_t iter = 0; iter < Max_Iter; ++iter) {
        double max_diff = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double old_xi = X[i];
            
            double off_diag_sum = 0.0;
            const std::size_t row_start = A.row_ptr[i];
            const std::size_t row_end = A.row_ptr[i+1];
            const std::size_t diag_idx = A.diag_indices[i];
            
            for (std::size_t j_idx = row_start; j_idx < diag_idx; ++j_idx) {
                off_diag_sum += A.values[j_idx] * X[A.col_indices[j_idx]];
            }
            for (std::size_t j_idx = diag_idx + 1; j_idx < row_end; ++j_idx) {
                off_diag_sum += A.values[j_idx] * X[A.col_indices[j_idx]];
            }
            
            const double new_xi_gs = (B[i] - off_diag_sum) * D_inv[i];
            const double new_xi = (1.0 - omega) * old_xi + omega * new_xi_gs;
            
            max_diff = std::max(max_diff, std::abs(new_xi - old_xi));
            X[i] = new_xi;
        }
        
        if (max_diff < tolerance) {
            std::cout << "Converged in " << iter + 1 << " iterations (L-infinity norm of change < " << tolerance << ")\n";
            return static_cast<int>(iter + 1);
        }
    }
    std::cerr << "Max iterations reached without convergence.\n";
    return -1;
}

std::vector<double> computeResidual(const SparseMatrix& A, const std::vector<double>& X, const std::vector<double>& B) {
    const std::size_t n = A.num_rows;
    std::vector<double> residual(n);
    for (std::size_t i = 0; i < n; ++i) {
        double ax_i = 0.0;
        const std::size_t row_start = A.row_ptr[i];
        const std::size_t row_end = A.row_ptr[i+1];
        for (std::size_t j_idx = row_start; j_idx < row_end; ++j_idx) {
            ax_i += A.values[j_idx] * X[A.col_indices[j_idx]];
        }
        residual[i] = B[i] - ax_i;
    }
    return residual;
}

double computeL2Norm(const std::vector<double>& vec) noexcept {
    double sum_sq = 0.0;
    for (const double val : vec) {
        sum_sq += val * val;
    }
    return std::sqrt(sum_sq);
}

void PrintExecutionTime(const std::chrono::steady_clock::time_point& start, const std::chrono::steady_clock::time_point& end) {
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << "ms\n";
}

void PrintVector(const std::string& label, const std::vector<double>& vec) {
    std::cout << label << ":\n";
    if (vec.size() > MAX_PRINT_SIZE) {
        std::cout << "[ ... vector of size " << vec.size() << " ... ]\n\n";
        return;
    }
    std::cout << "[";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        std::cout << std::fixed << std::setprecision(5) << vec[i] << (i < vec.size() - 1 ? ", " : "");
    }
    std::cout << "]\n\n";
}

void PrintMatrixInfo(const std::string& label, const SparseMatrix& A) {
    std::cout << label << ":\n"
              << "[ Sparse matrix " << A.num_rows << "x" << A.num_cols 
              << ", " << A.values.size() << " non-zero entries ]\n\n";
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cout.tie(nullptr);

    SparseMatrix A;
    bool matrix_is_tridiagonal;
    if (ReadMatrixMarketMatrix(A, matrix_is_tridiagonal) != 0) return 1;

    std::vector<double> B;
    if (ReadMatrixMarketVector(B) != 0) return 1;

    if (A.num_rows != B.size() || A.num_cols != A.num_rows) {
        std::cerr << "Error: Matrix must be square and dimensions compatible with vector.\n";
        return 1;
    }
    const std::size_t n = A.num_cols;

    PrintMatrixInfo("A", A);
    PrintVector("B", B);

    constexpr std::size_t Max_Iter = 1000;
    constexpr double tolerance = 1e-5;
    constexpr double omega = 1.25;

    if (matrix_is_tridiagonal) {
        std::cout << "Tridiagonal Method:\n";
        std::vector<double> X_trid;
        
        std::vector<double> lower(n > 0 ? n - 1 : 0, 0.0);
        std::vector<double> main_diag(n, 0.0);
        std::vector<double> upper(n > 0 ? n - 1 : 0, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j_idx = A.row_ptr[i]; j_idx < A.row_ptr[i+1]; ++j_idx) {
                const std::size_t j = A.col_indices[j_idx];
                const double val = A.values[j_idx];
                if (j == i) main_diag[i] = val;
                else if (j + 1 == i) lower[i-1] = val;
                else if (j == i + 1) upper[i] = val;
            }
        }
        
        const auto start = std::chrono::steady_clock::now();
        int status = TridiagonalSolve(n, std::move(lower), std::move(main_diag), std::move(upper), B, X_trid);
        const auto end = std::chrono::steady_clock::now();
        
        PrintExecutionTime(start, end);
        if (status == 0) {
            PrintVector("Solution", X_trid);
            const auto residual = computeResidual(A, X_trid, B);
            const double residual_norm = computeL2Norm(residual);
            std::cout << "Residual L2 Norm: " << residual_norm << "\n\n";
            if (WriteMatrixMarketVector("X_trid_out.mtx", X_trid) != 0) return 1;
        }
    } else {
        std::cout << "Matrix is not tridiagonal. Skipping Tridiagonal solver.\n\n";
    }

    auto run_iterative_solver = [&](const std::string& name, const char* out_filename, auto solve_fn) -> int {
        std::vector<double> X(n, 0.0);
        std::cout << name << ":\n";
        
        const auto start = std::chrono::steady_clock::now();
        int status = solve_fn(X);
        const auto end = std::chrono::steady_clock::now();
        
        PrintExecutionTime(start, end);
        if (status >= 0) {
            PrintVector("Solution", X);
            const auto residual = computeResidual(A, X, B);
            const double residual_norm = computeL2Norm(residual);
            std::cout << "Residual L2 Norm: " << residual_norm << "\n\n";
            if (WriteMatrixMarketVector(out_filename, X) != 0) return 1;
        }
        return 0;
    };

    if (run_iterative_solver("Jacobi Method", "X_jacobi_out.mtx", 
        [&](std::vector<double>& X){ return Jacobi(A, B, X, tolerance, Max_Iter); }) != 0) return 1;

    if (run_iterative_solver("\nGauss-Seidel Method", "X_gs_out.mtx", 
        [&](std::vector<double>& X){ return SOR(A, B, X, tolerance, Max_Iter, 1.0); }) != 0) return 1;
    
    char sor_label[128];
    snprintf(sor_label, sizeof(sor_label), "\nSOR Method (omega=%.3f)", omega);
    if (run_iterative_solver(sor_label, "X_sor_out.mtx", 
        [&](std::vector<double>& X){ return SOR(A, B, X, tolerance, Max_Iter, omega); }) != 0) return 1;

    return 0;
}
