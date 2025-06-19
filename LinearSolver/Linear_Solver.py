import numpy as np
import time
import math
import sys

def compute_residual(A, X, B):
    return A @ X - B

def compute_l2_norm(vec):
    return np.linalg.norm(vec)

def is_tridiagonal(A):
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        return False
    for i in range(n):
        for j in range(n):
            if (i > 0 and j < i - 1) or (j > i + 1):
                if A[i, j] != 0.0:
                    return False
    return True

def converged_using_l1_norm(X, X_prev, tolerance):
    distance_max = np.max(np.abs(X - X_prev))
    return distance_max < tolerance

def tridiagonal(A, B):
    n = A.shape[0]
    lower = np.zeros(n)
    main = np.zeros(n)
    upper = np.zeros(n)
    for i in range(n):
        main[i] = A[i, i]
        if i < n - 1:
            upper[i] = A[i, i + 1]
        if i > 0:
            lower[i] = A[i, i - 1]

    # Forward elimination
    for i in range(1, n):
        factor = lower[i] / main[i - 1]
        main[i] -= factor * upper[i - 1]
        B[i] -= factor * B[i - 1]

    X = np.zeros(n)
    X[-1] = B[-1] / main[-1]

    # Back substitution
    for i in range(n - 2, -1, -1):
        X[i] = (B[i] - upper[i] * X[i + 1]) / main[i]

    write_tridiagonal_to_matrix_market_file(B, X)
    return X

def jacobi_l1_norm(A, B, X, tolerance, max_iter):
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        print("Error: Jacobi method requires a square matrix.", file=sys.stderr)
        return -1, X
    X_prev = np.zeros_like(X)
    iter_count = 0
    B_norm = compute_l2_norm(B)
    rel_tolerance = tolerance * B_norm

    while iter_count < max_iter:
        X_prev[:] = X
        for i in range(n):
            s = np.dot(A[i, :], X_prev) - A[i, i] * X_prev[i]
            X[i] = (B[i] - s) / A[i, i]
        residual = compute_residual(A, X, B)
        residual_norm = compute_l2_norm(residual)
        iter_count += 1
        if residual_norm <= rel_tolerance:
            print(f"Converged in {iter_count} iterations. (Relative L2 residual: {residual_norm / B_norm} < {tolerance})")
            return iter_count, X
    print(f"Max iterations reached. (Relative L2 residual: {residual_norm / B_norm} > {tolerance})", file=sys.stderr)
    return -1, X

def gauss_seidel_l1_norm(A, B, X, tolerance, max_iter):
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        print("Error: Matrix must be square.", file=sys.stderr)
        return -1, X
    X_prev = np.zeros_like(X)
    iter_count = 0
    B_norm = compute_l2_norm(B)
    rel_tolerance = tolerance * B_norm

    while iter_count < max_iter:
        X_prev[:] = X
        for i in range(n):
            s1 = np.dot(A[i, :i], X[:i])
            s2 = np.dot(A[i, i+1:], X_prev[i+1:])
            X[i] = (B[i] - s1 - s2) / A[i, i]
        residual = compute_residual(A, X, B)
        residual_norm = compute_l2_norm(residual)
        iter_count += 1
        if residual_norm <= rel_tolerance:
            print(f"Converged in {iter_count} iterations. (Relative L2 residual: {residual_norm / B_norm} < {tolerance})")
            return iter_count, X
    print(f"Max iterations reached. (Relative L2 residual: {residual_norm / B_norm} > {tolerance})", file=sys.stderr)
    return -1, X

def sor_l1_norm(A, B, X, tolerance, max_iter, omega):
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        print("Error: Matrix must be square.", file=sys.stderr)
        return -1, X
    X_prev = np.zeros_like(X)
    iter_count = 0
    B_norm = compute_l2_norm(B)
    rel_tolerance = tolerance * B_norm

    while iter_count < max_iter:
        X_prev[:] = X
        for i in range(n):
            s1 = np.dot(A[i, :i], X[:i])
            s2 = np.dot(A[i, i+1:], X_prev[i+1:])
            gs_update = (B[i] - s1 - s2) / A[i, i]
            X[i] = (1 - omega) * X_prev[i] + omega * gs_update
        residual = compute_residual(A, X, B)
        residual_norm = compute_l2_norm(residual)
        iter_count += 1
        if residual_norm <= rel_tolerance:
            print(f"Converged in {iter_count} iterations. (Relative L2 residual: {residual_norm / B_norm} < {tolerance})")
            return iter_count, X
    print(f"Max iterations reached. (Relative L2 residual: {residual_norm / B_norm} > {tolerance})", file=sys.stderr)
    return -1, X

def write_1d_array_to_matrix_market_file(B):
    with open("B_out.dat", "w") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write(f"1 {len(B)} {len(B)}\n")
        for i, val in enumerate(B):
            f.write(f"1 {i+1} {val:.6f}\n")

def write_2d_array_to_matrix_market_file(array_A, num_rows, num_cols):
    with open("A_out.dat", "w") as f:
        f.write("%%%%MatrixMarket matrix coordinate real general\n")
        f.write(f"{num_rows} {num_cols} {num_rows * num_cols}\n")
        for i in range(num_rows):
            for j in range(num_cols):
                val = array_A[i * num_cols + j]
                f.write(f"{i+1} {j+1} {val:.6e}\n")

def write_matrix_market_matrix(A):
    entries = []
    num_rows, num_cols = A.shape
    for i in range(num_rows):
        for j in range(num_cols):
            if A[i, j] != 0.0:
                entries.append((i+1, j+1, A[i, j]))
    with open("A_out.dat", "w") as f:
        f.write("%%%%MatrixMarket matrix coordinate real general\n")
        f.write(f"{num_rows} {num_cols} {len(entries)}\n")
        for (r, c, v) in entries:
            f.write(f"{r} {c} {v:.15g}\n")

def write_matrix_market_vector(X):
    with open("X_out.dat", "w") as f:
        f.write("%%%%MatrixMarket matrix coordinate real general\n")
        f.write(f"1 {len(X)} {len(X)}\n")
        for i, val in enumerate(X):
            f.write(f"1 {i+1} {val:.15g}\n")

def write_tridiagonal_to_matrix_market_file(B, X):
    write_1d_array_to_matrix_market_file(B)
    write_matrix_market_vector(X)

def read_matrix_market_matrix(filename="A.dat"):
    with open(filename, "r") as f:
        line = f.readline()
        while line.startswith('%'):
            line = f.readline()
        parts = line.strip().split()
        if len(parts) < 3:
            raise ValueError("Invalid header in matrix file")
        num_rows, num_cols, num_entries = map(int, parts)
        A = np.zeros((num_rows, num_cols))
        for _ in range(num_entries):
            line = f.readline()
            r, c, val = line.strip().split()
            r, c = int(r) - 1, int(c) - 1
            A[r, c] = float(val)
    return A

def read_matrix_market_vector(filename="B.dat"):
    with open(filename, "r") as f:
        line = f.readline()
        while line.startswith('%'):
            line = f.readline()
        parts = line.strip().split()
        if len(parts) == 3:
            # coordinate format
            rows, cols, entries = map(int, parts)
            if rows == 1:
                num_rows = cols
            else:
                num_rows = rows
            B = np.zeros(num_rows)
            for _ in range(entries):
                line = f.readline()
                r, c, val = line.strip().split()
                r, c = int(r), int(c)
                val = float(val)
                if r == 1:
                    B[c - 1] = val
                else:
                    B[r - 1] = val
        else:
            # array format
            rows, cols = map(int, parts)
            if cols != 1:
                raise ValueError("Vector must have 1 column")
            B = np.zeros(rows)
            for i in range(rows):
                val = float(f.readline().strip())
                B[i] = val
    return B

def print_execution_time(start, end):
    elapsed = end - start
    print(f"Elapsed time: {elapsed:.6f}s")

def print_vector(label, vec):
    print(f"{label}:\n[{', '.join(f'{x:.5f}' for x in vec)}]\n")

def print_matrix(label, A):
    print(f"{label}:")
    for row in A:
        print(' '.join(f"{v:10.5f}" for v in row))
    print()

def main():
    omega = 1.25
    try:
        A = read_matrix_market_matrix("A.dat")
    except Exception as e:
        print(f"Failed to read matrix A: {e}", file=sys.stderr)
        return -1

    try:
        B = read_matrix_market_vector("B.dat")
    except Exception as e:
        print(f"Failed to read vector B: {e}", file=sys.stderr)
        return -1

    if A.shape[0] != B.shape[0]:
        print("Error: Matrix and vector dimensions do not match.", file=sys.stderr)
        return -1

    flatA = A.flatten()
    write_2d_array_to_matrix_market_file(flatA, A.shape[0], A.shape[1])
    print()

    print_vector("B", B)

    Max_Iter = 1000
    tolerance = 1e-5

    # Tridiagonal Solver
    if is_tridiagonal(A):
        B_trid = B.copy()
        start = time.time()
        X_trid = tridiagonal(A, B_trid)
        end = time.time()
        print("Tridiagonal Method:")
        print_execution_time(start, end)
        print_vector("Solution", X_trid)
        residual = compute_residual(A, X_trid, B)
        residual_norm = compute_l2_norm(residual)
        print(f"Residual L2 Norm: {residual_norm}\n")
    else:
        print("Matrix is not tridiagonal. Skipping Tridiagonal solver.\n")

    # Jacobi L1 Norm Solver
    B_jacobi_l1 = B.copy()
    X_jacobi_l1 = np.zeros_like(B)
    print("Jacobi (Residual Norm) Method:")
    start = time.time()
    status, X_jacobi_l1 = jacobi_l1_norm(A, B_jacobi_l1, X_jacobi_l1, tolerance, Max_Iter)
    end = time.time()
    print_execution_time(start, end)
    if status >= 0:
        print_vector("Solution", X_jacobi_l1)
    else:
        print("Jacobi (Residual Norm) did not converge.", file=sys.stderr)

    # Gauss-Seidel Solver
    B_gs_l1 = B.copy()
    X_gs_l1 = np.zeros_like(B)
    print("\nGauss-Seidel (Residual Norm) Method:")
    start = time.time()
    status, X_gs_l1 = gauss_seidel_l1_norm(A, B_gs_l1, X_gs_l1, tolerance, Max_Iter)
    end = time.time()
    print_execution_time(start, end)
    if status >= 0:
        print_vector("Solution", X_gs_l1)
    else:
        print("Gauss-Seidel (Residual Norm) did not converge.", file=sys.stderr)

    # SOR Solver
    B_sor_l1 = B.copy()
    X_sor_l1 = np.zeros_like(B)
    print("\nSOR (Residual Norm) Method:")
    start = time.time()
    status, X_sor_l1 = sor_l1_norm(A, B_sor_l1, X_sor_l1, tolerance, Max_Iter, omega)
    end = time.time()
    print_execution_time(start, end)
    if status >= 0:
        print_vector("Solution", X_sor_l1)
    else:
        print("SOR (Residual Norm) did not converge.", file=sys.stderr)

    return 0

if __name__ == "__main__":
    main()
