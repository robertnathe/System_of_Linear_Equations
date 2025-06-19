import numpy as np
import time
import sys

Matrix = np.ndarray

def print_execution_time(start, end):
    elapsed = end - start
    print(f"Elapsed time: {elapsed:.6f}s\n")

def print_vector_2d(temp_matrix: Matrix) -> None:
    print("Displaying the 2D vector:")
    print(np.array2string(temp_matrix, precision=5, suppress_small=True))

def print_vector(label: str, vec: np.ndarray) -> None:
    print(f"{label}:\n{vec}")

def read_matrix_market_matrix(file_path: str) -> tuple[Matrix, int, int]:
    with open(file_path, 'r') as file:
        # Skip comments and read the header
        header = next(line for line in file if not line.startswith('%'))
        num_rows, num_cols, num_entries = map(int, header.split())
        matrix = np.zeros((num_rows, num_cols))
        
        # Read matrix entries efficiently
        for _ in range(num_entries):
            row, col, value = map(float, file.readline().split())
            matrix[int(row) - 1, int(col) - 1] = value
            
    return matrix, num_rows, num_cols

def read_matrix_market_vector(file_path: str) -> tuple[np.ndarray, int]:
    with open(file_path, 'r') as file:
        header = next(line for line in file if not line.startswith('%'))
        declared_rows, declared_cols, num_entries = map(int, header.split())

        print(f"Declared rows: {declared_rows}, Declared cols: {declared_cols}, Num entries: {num_entries}")  # Debug output
        
        if declared_rows == 1 and declared_cols > 1:  # Check for row vector
            vector = np.zeros(declared_cols)
            for _ in range(num_entries):
                _, col, value = map(float, file.readline().split())
                vector[int(col) - 1] = value
            return vector, declared_cols
        else:
            raise ValueError("Expected a row vector.")

def write_matrix_market_file(file_path: str, data: np.ndarray, is_vector: bool = False) -> None:
    with open(file_path, 'w') as outfile:
        if is_vector:
            # Write as a row vector
            outfile.write("%%MatrixMarket matrix coordinate real general\n")
            outfile.write(f"1 {len(data)} {len(data)}\n")  # Adjusting to write as a row vector
            for i, value in enumerate(data, start=1):
                outfile.write(f"1 {i} {value:.6f}\n")
        else:
            num_rows, num_cols = data.shape
            outfile.write("%%MatrixMarket Matrix coordinate real general\n")
            outfile.write(f"{num_rows} {num_cols} {np.count_nonzero(data)}\n")
            for i in range(num_rows):
                for j in range(num_cols):
                    if data[i, j] != 0.0:  # Only write non-zero entries
                        outfile.write(f"{i+1} {j+1} {data[i, j]:.6f}\n")

def write_tridiagonal_to_vector_market_file(B: np.ndarray, X: np.ndarray) -> None:
    write_matrix_market_file("B_out.dat", B, is_vector=True)
    write_matrix_market_file("X_out.dat", X, is_vector=True)

def converged_using_l1_norm(X, X_previous, tolerance):
    distance_max = np.max(np.abs(X - X_previous))
    return distance_max < tolerance

def tridiagonal(A, B, X):
    num_rows, num_cols = A.shape
    if num_rows != num_cols:
        raise ValueError("Matrix must be square")

    lower_diagonal = np.copy(np.diag(A, -1))  # Create a copy of the lower diagonal
    main_diagonal = np.copy(np.diag(A))        # Create a copy of the main diagonal
    upper_diagonal = np.copy(np.diag(A, 1))    # Create a copy of the upper diagonal

    for i in range(1, num_rows):  # Iterate based on num_rows
        factor = lower_diagonal[i - 1] / main_diagonal[i - 1]  # Adjusted to i-1
        main_diagonal[i] -= factor * upper_diagonal[i - 1]
        B[i] -= factor * B[i - 1]

    X[-1] = B[-1] / main_diagonal[-1]
    for i in range(num_rows - 2, -1, -1):
        X[i] = (B[i] - upper_diagonal[i] * X[i + 1]) / main_diagonal[i]

    return X

def tridiagonal_lu_decomposition(lower_diagonal, main_diagonal, upper_diagonal, num_cols):
    for i in range(num_cols - 1):
        if main_diagonal[i] == 0.0:
            return -1
        lower_diagonal[i] /= main_diagonal[i]
        main_diagonal[i + 1] -= lower_diagonal[i] * upper_diagonal[i]
    if main_diagonal[num_cols - 1] == 0.0:
        return -1
    return 0

def tridiagonal_lu_solve(lower_diagonal, main_diagonal, upper_diagonal, B, num_cols):
    X = B[:]
    for i in range(num_cols):
        if main_diagonal[i] == 0.0:
            raise ValueError("Division by zero")
        if i > 0:
            X[i] -= lower_diagonal[i - 1] * X[i - 1]
    X[-1] /= main_diagonal[-1]
    for i in range(num_cols - 1, -1, -1):
        if main_diagonal[i] == 0.0:
            raise ValueError("Division by zero")
        X[i] /= main_diagonal[i]
        if i < num_cols - 1:
            X[i] -= upper_diagonal[i] * X[i + 1]
    return X

def jacobi_l1_norm(A, B, X, num_rows, num_cols, tolerance, max_iter):
    if num_rows != num_cols:
        print("Error: Jacobi method requires a square matrix.")
        return -1
    if A.shape[0] != num_rows or (num_rows > 0 and A.shape[1] != num_cols):
        print("Error: Matrix dimensions mismatch.")
        return -1
    if B.size != num_rows:
        print("Error: Vector B size does not match matrix dimensions.")
        return -1

    X_prev = np.copy(X)
    last_max_diff = 0.0
    converged = False

    for iter in range(max_iter):
        X_new = np.zeros(num_rows)
        current_max_diff = 0.0
        
        for i in range(num_rows):
            sum_val = np.sum(A[i, :] * X_prev) - A[i, i] * X_prev[i]
            X_new[i] = (B[i] - sum_val) / A[i, i]
            current_max_diff = max(current_max_diff, abs(X_new[i] - X_prev[i]))
        
        last_max_diff = current_max_diff
        X_prev = X_new
        
        if last_max_diff < tolerance:
            converged = True
            X[:] = X_new  # Update X with the new values
            print(f"Converged in {iter + 1} iterations (L∞ norm: {last_max_diff} < {tolerance})")
            break

    if not converged:
        X[:] = X_prev  # Restore previous values
        print(f"Max iterations reached (L∞ norm: {last_max_diff} > {tolerance})")
        return 1

    return 0

def jacobi_l2_norm(A, B, X, num_rows, num_cols, tolerance, max_iter):
    if num_rows != num_cols or num_rows == 0:
        print("Error: Matrix must be square.")
        return -1

    X_prev = np.zeros(num_rows)
    iter_count = 0
    residual = 0.0

    while True:
        residual = 0.0
        X_prev[:] = X  # Copy previous values
        
        for i in range(num_rows):
            sum_val = np.sum(A[i, :] * X_prev) - A[i, i] * X_prev[i]
            X[i] = (B[i] - sum_val) / A[i, i]
            residual += (X[i] - X_prev[i]) ** 2
        
        residual = np.sqrt(residual)
        iter_count += 1
        
        if residual <= tolerance or iter_count >= max_iter:
            break

    if residual <= tolerance:
        print(f"Converged in {iter_count} iterations. (L2 norm: {residual} < {tolerance})")
        return iter_count
    else:
        print("Max iterations reached.")
        return -1
        
def gauss_seidel_l1_norm(A, B, X, num_rows, num_cols, tolerance, max_iter):
    num_rows, num_cols = A.shape
    if num_rows != num_cols:
        raise ValueError("Matrix must be square")

    inv_diag = 1 / np.diag(A)
    for _ in range(max_iter):
        X_prev = np.copy(X)
        max_diff = 0.0
        for i in range(num_rows):
            s = sum(A[i, j] * X[j] for j in range(num_cols) if j != i)
            X[i] = (B[i] - s) * inv_diag[i]
            max_diff = max(max_diff, np.abs(X[i] - X_prev[i]))
        if max_diff < tolerance:
            print(f"Converged in {_+1} iterations (L∞ norm: {max_diff} < {tolerance})")
            return X, 0  # Return the solution and status code
    print(f"Max iterations reached (L∞ norm: {max_diff} > {tolerance})")
    return X, 1  # Return the solution and status code

def sor_l1_norm(A, B, X, num_rows, num_cols, tolerance, max_iter, omega):
    num_rows, num_cols = A.shape
    if num_rows != num_cols:
        raise ValueError("Matrix must be square")

    inv_diag = 1 / np.diag(A)
    for _ in range(max_iter):
        X_prev = np.copy(X)
        max_diff = 0.0
        for i in range(num_rows):
            s = sum(A[i, j] * X[j] for j in range(num_cols) if j != i)
            gs_update = (B[i] - s) * inv_diag[i]
            X[i] = (1.0 - omega) * X_prev[i] + omega * gs_update
            max_diff = max(max_diff, np.abs(X[i] - X_prev[i]))
        if max_diff < tolerance:
            print(f"Converged in {_+1} iterations (L∞ norm: {max_diff} < {tolerance})")
            return X, 0  # Return the solution and status code
    print(f"Max iterations reached (L∞ norm: {max_diff} > {tolerance})")
    return X, 1  # Return the solution and status code

def read_matrix_market_matrix(file_path):
    with open(file_path, 'r') as file:
        # Skip comments and read the header
        header = next(line for line in file if not line.startswith('%'))
        num_rows, num_cols, num_entries = map(int, header.split())
        A = np.zeros((num_rows, num_cols))
        
        for _ in range(num_entries):
            row, col, value = map(float, file.readline().split())
            A[int(row) - 1, int(col) - 1] = value
            
    return A, num_rows, num_cols

def read_matrix_market_vector(file_path):
    with open(file_path, 'r') as file:
        header = next(line for line in file if not line.startswith('%'))
        declared_rows, declared_cols, num_entries = map(int, header.split())
        
        if declared_rows == 1 and declared_cols > 1:  # Row vector
            B = np.zeros(declared_cols)
            for _ in range(num_entries):
                _, col, value = map(float, file.readline().split())
                B[int(col) - 1] = value
            return B, declared_cols
        else:
            raise ValueError("Expected a row vector.")

def print_vector(label, vec):
    print(f"{label}:\n{vec}")

def print_execution_time(start, end):
    elapsed = end - start
    print(f"Elapsed time: {elapsed:.6f}s")

def main():
    omega = 1.25
    tolerance = 0.00001
    max_iter = 1000
    A, num_rows, num_cols = read_matrix_market_matrix("A.dat")
    B, vector_rows = read_matrix_market_vector("B.dat")

    if num_rows != vector_rows:
        print("Error: Matrix and vector dimensions do not match.")
        return

    # --- Tridiagonal Solver ---
    B_trid = B.copy()
    X_trid = np.zeros(num_cols)
    start = time.time()
    X_trid = tridiagonal(A, B_trid, X_trid)  # Corrected function call
    end = time.time()
    print("Tridiagonal Method:")
    print_execution_time(start, end)
    print_vector("Solution", X_trid)

    # --- Jacobi L1 Norm Solver ---
    B_jacobi_l1 = B.copy()
    X_jacobi_l1 = np.zeros(num_cols)
    print("\nJacobi (L1 Norm) Method:")
    start = time.time()
    status = jacobi_l1_norm(A, B_jacobi_l1, X_jacobi_l1, num_rows, num_cols, tolerance, max_iter)
    end = time.time()
    print_execution_time(start, end)
    if status == 0:
        print_vector("Solution", X_jacobi_l1)
    else:
        print("Jacobi (L1 Norm) did not converge.")

    # --- Jacobi L2 Norm Solver ---
    B_jacobi_l2 = B.copy()
    X_jacobi_l2 = np.zeros(num_cols)
    print("\nJacobi (L2 Norm) Method:")
    start = time.time()
    status = jacobi_l2_norm(A, B_jacobi_l2, X_jacobi_l2, num_rows, num_cols, tolerance, max_iter)  # Implement this function
    end = time.time()
    print_execution_time(start, end)
    if status >= 0:
        print_vector("Solution", X_jacobi_l2)
    else:
        print("Jacobi (L2 Norm) did not converge.")

    # --- Gauss-Seidel L1 Norm Solver ---
    B_gs_l1 = B.copy()
    X_gs_l1 = np.zeros(num_cols)
    print("\nGauss-Seidel (L1 Norm) Method:")
    start = time.time()
    X_gs_l1, status = gauss_seidel_l1_norm(A, B_gs_l1, X_gs_l1, num_rows, num_cols, tolerance, max_iter)  # Unpack the values
    end = time.time()
    print_execution_time(start, end)
    if status == 0:
        print_vector("Solution", X_gs_l1)
    else:
        print("Gauss-Seidel (L1 Norm) did not converge.")

    # --- SOR L1 Norm Solver ---
    B_sor_l1 = B.copy()
    X_sor_l1 = np.zeros(num_cols)
    print("\nSOR (L1 Norm) Method:")
    start = time.time()
    X_sor_l1, status = sor_l1_norm(A, B_sor_l1, X_sor_l1, num_rows, num_cols, tolerance, max_iter, omega)  # Unpack the values
    end = time.time()
    print_execution_time(start, end)
    if status == 0:
        print_vector("Solution", X_sor_l1)
    else:
        print("SOR (L1 Norm) did not converge.")

if __name__ == "__main__":
    main()