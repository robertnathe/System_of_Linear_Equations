"""
An optimized Python program for solving systems of linear equations (AX=B)
using direct (Tridiagonal) and iterative (Jacobi, Gauss-Seidel, SOR) methods
for sparse matrices in Matrix Market format, leveraging a CSR data structure.
"""
import sys
import time
from typing import Tuple, Optional, Dict, Any, Callable

import numpy as np
import numpy.typing as npt
import numba

# --- Type Aliases for Clarity ---
# CSR sparse matrix format: a tuple holding indptr, indices, data, and shape
SparseCsr = Tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64], Tuple[int, int]]
# Solver function type for iterative methods
IterativeSolverFunc = Callable[..., Tuple[int, npt.NDArray[np.float64]]]

# --- Constants ---
ZERO_TOL = 1e-12

# --- Sparse Matrix Operations using CSR ---

def _coo_to_csr(
    rows: npt.NDArray[np.int32],
    cols: npt.NDArray[np.int32],
    data: npt.NDArray[np.float64],
    shape: Tuple[int, int]
) -> SparseCsr:
    """
    Converts a COO representation to CSR, summing duplicate entries efficiently.
    This is a high-performance implementation using vectorized NumPy operations.
    """
    num_rows = shape[0]

    # Sort by row, then by column. This groups all entries for a given row together.
    sorter = np.lexsort((cols, rows))
    rows, cols, data = rows[sorter], cols[sorter], data[sorter]

    # Identify unique (row, col) pairs in the sorted data to sum duplicates.
    # A change in row or col indicates a new unique element.
    unique_mask = np.concatenate(([True], (rows[1:] != rows[:-1]) | (cols[1:] != cols[:-1])))

    # The 'indptr' for the summed data can be found using the unique mask's locations.
    sum_indptr = np.concatenate((np.flatnonzero(unique_mask), [len(data)]))

    # Sum the data for each unique (row, col) using the boundaries found.
    summed_data = np.add.reduceat(data, sum_indptr[:-1])

    # Get the unique rows and cols corresponding to the summed data.
    unique_rows = rows[unique_mask]
    unique_cols = cols[unique_mask]

    # Now, build the final CSR indptr from the unique, sorted rows.
    # This counts the number of non-zero elements in each row.
    indptr = np.zeros(num_rows + 1, dtype=np.int32)
    np.add.at(indptr, unique_rows + 1, 1)
    indptr = np.cumsum(indptr)

    return indptr, unique_cols.astype(np.int32), summed_data, shape


@numba.jit(nopython=True, fastmath=True, cache=True)
def _csr_matvec_kernel(
    indptr: npt.NDArray[np.int32],
    indices: npt.NDArray[np.int32],
    data: npt.NDArray[np.float64],
    X: npt.NDArray[np.float64],
    Y: npt.NDArray[np.float64],
) -> None:
    """Computes Y = AX using CSR format. Y is an output parameter."""
    for i in range(len(indptr) - 1):
        y_i = 0.0
        for j_ptr in range(indptr[i], indptr[i+1]):
            y_i += data[j_ptr] * X[indices[j_ptr]]
        Y[i] = y_i

def _compute_residual_csr(
    A_csr: SparseCsr,
    X: npt.NDArray[np.float64],
    B: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """Computes the residual vector (B - AX) for a sparse matrix in CSR format."""
    Y = np.zeros_like(B)
    _csr_matvec_kernel(A_csr[0], A_csr[1], A_csr[2], X, Y)
    return B - Y

def is_tridiagonal_csr(A_csr: SparseCsr) -> bool:
    """Checks if a sparse matrix in CSR format is tridiagonal using vectorization."""
    indptr, indices, _, shape = A_csr
    n_rows = shape[0]
    # Efficiently get the row index for each element in 'indices'/'data'
    rows = np.repeat(np.arange(n_rows, dtype=np.int32), np.diff(indptr))
    return np.all(np.abs(rows - indices) <= 1)

def _extract_tridiagonals_csr(
    A_csr: SparseCsr
) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Extracts main, upper, and lower diagonals from a CSR matrix."""
    indptr, indices, data, shape = A_csr
    n = shape[0]
    main_diag = np.zeros(n, dtype=np.float64)
    upper_diag = np.zeros(n - 1, dtype=np.float64)
    lower_diag = np.zeros(n - 1, dtype=np.float64)

    for i in range(n):
        for j_ptr in range(indptr[i], indptr[i+1]):
            col = indices[j_ptr]
            val = data[j_ptr]
            if col == i:
                main_diag[i] += val
            elif col == i + 1:
                if i < n - 1:
                    upper_diag[i] += val
            elif col == i - 1:
                if i > 0:
                    lower_diag[i - 1] += val
    return main_diag, upper_diag, lower_diag

def _extract_diagonal_csr(A_csr: SparseCsr) -> npt.NDArray[np.float64]:
    """Extracts the main diagonal from a CSR matrix."""
    indptr, indices, data, shape = A_csr
    n = shape[0]
    diag = np.zeros(n, dtype=np.float64)
    for i in range(n):
        for j_ptr in range(indptr[i], indptr[i+1]):
            if indices[j_ptr] == i:
                diag[i] += data[j_ptr]
                break # Assume one diagonal entry per row after summing duplicates
    return diag

# --- Linear Equation Solvers ---

def solve_tridiagonal(
    A_csr: SparseCsr,
    B: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """Solves a tridiagonal system AX = B using the Thomas algorithm."""
    n = A_csr[3][0]
    main_diag, upper_diag, lower_diag = _extract_tridiagonals_csr(A_csr)
    a, b, c = lower_diag.copy(), main_diag.copy(), upper_diag.copy()
    d = B.copy()

    for i in range(1, n):
        if abs(b[i - 1]) < ZERO_TOL:
            raise ValueError(f"Zero pivot at index {i-1} during forward elimination.")
        factor = a[i - 1] / b[i - 1]
        b[i] -= factor * c[i - 1]
        d[i] -= factor * d[i - 1]

    if abs(b[-1]) < ZERO_TOL:
        raise ValueError("Zero pivot at last element during back substitution.")
    X = np.zeros(n, dtype=np.float64)
    X[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        if abs(b[i]) < ZERO_TOL:
             raise ValueError(f"Zero pivot at index {i} during back substitution.")
        X[i] = (d[i] - c[i] * X[i + 1]) / b[i]

    return X

@numba.jit(nopython=True, fastmath=True, cache=True)
def _jacobi_kernel(
    indptr: npt.NDArray[np.int32],
    indices: npt.NDArray[np.int32],
    data: npt.NDArray[np.float64],
    D_inv: npt.NDArray[np.float64],
    B: npt.NDArray[np.float64],
    X_old: npt.NDArray[np.float64],
    X_new: npt.NDArray[np.float64]
) -> None:
    """Numba-jitted kernel for a single Jacobi sweep using CSR format."""
    n = len(B)
    for i in range(n):
        off_diag_sum = 0.0
        for j_ptr in range(indptr[i], indptr[i+1]):
            col = indices[j_ptr]
            if col != i:
                off_diag_sum += data[j_ptr] * X_old[col]
        X_new[i] = (B[i] - off_diag_sum) * D_inv[i]

def solve_jacobi(
    A_csr: SparseCsr, D: npt.NDArray[np.float64], B: npt.NDArray[np.float64], X0: npt.NDArray[np.float64],
    tolerance: float, max_iter: int, b_norm: float
) -> Tuple[int, npt.NDArray[np.float64]]:
    """Solves AX = B using the Jacobi method for a sparse CSR matrix."""
    if np.any(np.abs(D) < ZERO_TOL):
        raise ValueError("Zero or near-zero diagonal element found.")
    D_inv = 1.0 / D

    X_old = X0.copy()
    X_new = np.empty_like(X_old)
    rel_tolerance_norm = tolerance * b_norm
    indptr, indices, data, _ = A_csr

    for iter_count in range(1, max_iter + 1):
        _jacobi_kernel(indptr, indices, data, D_inv, B, X_old, X_new)
        residual_norm = np.linalg.norm(_compute_residual_csr(A_csr, X_new, B))
        if residual_norm <= rel_tolerance_norm:
            return iter_count, X_new
        X_old, X_new = X_new, X_old # Swap buffers for next iteration

    return -1, X_old

@numba.jit(nopython=True, fastmath=True, cache=True)
def _sor_gs_kernel(
    indptr: npt.NDArray[np.int32],
    indices: npt.NDArray[np.int32],
    data: npt.NDArray[np.float64],
    D_inv: npt.NDArray[np.float64],
    B: npt.NDArray[np.float64],
    X: npt.NDArray[np.float64],
    omega: float,
) -> None:
    """Numba-jitted kernel for a single SOR/Gauss-Seidel sweep using CSR."""
    n = len(B)
    for i in range(n):
        off_diag_sum = 0.0
        for j_ptr in range(indptr[i], indptr[i+1]):
            col = indices[j_ptr]
            if i != col:
                off_diag_sum += data[j_ptr] * X[col]
        update = (B[i] - off_diag_sum) * D_inv[i]
        X[i] = (1.0 - omega) * X[i] + omega * update

def solve_sor_gs(
    A_csr: SparseCsr, D: npt.NDArray[np.float64], B: npt.NDArray[np.float64], X0: npt.NDArray[np.float64],
    tolerance: float, max_iter: int, omega: float, b_norm: float
) -> Tuple[int, npt.NDArray[np.float64]]:
    """Solves AX=B using SOR (omega > 1), or Gauss-Seidel (omega=1) with CSR."""
    X = X0.copy()
    if np.any(np.abs(D) < ZERO_TOL):
        raise ValueError("Zero or near-zero diagonal element found.")
    D_inv = 1.0 / D
    indptr, indices, data, _ = A_csr
    rel_tolerance_norm = tolerance * b_norm

    for iter_count in range(1, max_iter + 1):
        _sor_gs_kernel(indptr, indices, data, D_inv, B, X, omega)
        residual_norm = np.linalg.norm(_compute_residual_csr(A_csr, X, B))
        if residual_norm <= rel_tolerance_norm:
            return iter_count, X
    return -1, X

# --- File I/O ---
def write_mm_vector(filename: str, vec: npt.NDArray[np.float64]) -> None:
    """Writes a NumPy vector to a file in Matrix Market coordinate format."""
    try:
        n = len(vec)
        with open(filename, "w", encoding="utf-8") as f:
            # Header for a coordinate format vector (N rows, 1 col, N non-zero entries)
            f.write("%%MatrixMarket matrix coordinate real general\n")
            f.write(f"{n} 1 {n}\n")
            
            # Write row (i+1), column (1), and value
            for i in range(n):
                f.write(f"{i + 1} 1 {vec[i]:.15g}\n")
    except IOError as e:
        print(f"Error writing vector to file {filename}: {e}", file=sys.stderr)

def read_mm_matrix_to_coo(filename: str) -> Optional[Tuple[npt.NDArray, npt.NDArray, npt.NDArray, Tuple[int, int]]]:
    """Reads a matrix from a file in Matrix Market coordinate format into COO arrays."""
    try:
        with open(filename, "r", encoding="utf-8") as f:
            banner = f.readline()
            if not (banner.lower().startswith("%%matrixmarket") and "coordinate" in banner.lower()):
                raise ValueError("Unsupported matrix format: expected MatrixMarket coordinate.")
            line = f.readline()
            while line.startswith('%'):
                line = f.readline()
            num_rows, num_cols, _ = map(int, line.split())
            data = np.loadtxt(f, dtype=[('r', 'i4'), ('c', 'i4'), ('v', 'f8')], ndmin=1)
            return (data['r'] - 1, data['c'] - 1, data['v'], (num_rows, num_cols))
    except (IOError, ValueError, IndexError) as e:
        print(f"Error reading matrix file {filename}: {e}", file=sys.stderr)
        return None

def read_mm_vector(filename: str) -> Optional[npt.NDArray[np.float64]]:
    """Reads a vector from a file in Matrix Market array or coordinate format."""
    try:
        with open(filename, "r", encoding="utf-8") as f:
            banner = f.readline()
            if not banner.lower().startswith("%%matrixmarket"):
                raise ValueError("Not a valid Matrix Market file.")
            line = f.readline()
            while line.startswith('%'):
                line = f.readline()
            if "array" in banner.lower():
                return np.loadtxt(f, dtype=np.float64, ndmin=1)
            else: # Handle coordinate format vectors as dense
                num_rows, num_cols, _ = map(int, line.split())
                size = max(num_rows, num_cols)
                B = np.zeros(size, dtype=np.float64)
                data = np.loadtxt(f, dtype=[('r', 'i4'), ('c', 'i4'), ('v', 'f8')], ndmin=1)
                if data.size > 0:
                    indices = data['r'] - 1 if num_cols == 1 else data['c'] - 1
                    B[indices] = data['v']
                return B
    except (IOError, ValueError, IndexError) as e:
        print(f"Error reading vector file {filename}: {e}", file=sys.stderr)
        return None

# --- Main Execution ---

def print_vector(label: str, vec: npt.NDArray[np.float64]) -> None:
    """Prints a labeled vector with formatted floating-point numbers."""
    print(f"{label}:")
    if vec.size == 0:
        print("[]\n")
        return
    if vec.size > 10:
        head = ", ".join(f'{x:.5f}' for x in vec[:5])
        tail = ", ".join(f'{x:.5f}' for x in vec[-5:])
        print(f"[{head}, ..., {tail}]\n")
    else:
        print(f"[{', '.join(f'{x:.5f}' for x in vec)}]\n")

def main() -> int:
    """Main function to read matrices, run solvers, and report results."""
    try:
        coo_matrix_data = read_mm_matrix_to_coo("A.mtx")
        if coo_matrix_data is None: return 1
        B = read_mm_vector("B.mtx")
        if B is None: return 1
    except FileNotFoundError as e:
        print(f"Error: Required file not found - {e}", file=sys.stderr)
        return 1

    coo_rows, coo_cols, coo_data, shape = coo_matrix_data
    num_rows, num_cols = shape
    if num_rows != num_cols:
        print("Error: Input matrix 'A' must be square.", file=sys.stderr)
        return 1
    if num_rows != B.shape[0]:
        print(f"Error: Dimension mismatch. A is {shape}, B is {B.shape}.", file=sys.stderr)
        return 1
    print_vector("B", B)

    # --- Pre-processing: Convert to CSR ---
    print("Pre-processing matrix from COO to CSR format...")
    start_time = time.perf_counter()
    A_csr = _coo_to_csr(coo_rows, coo_cols, coo_data, shape)
    D = _extract_diagonal_csr(A_csr)
    end_time = time.perf_counter()
    print(f"Pre-processing time: {end_time - start_time:.6f}s\n")

    # --- Configuration ---
    MAX_ITER = 1000
    TOLERANCE = 1e-5
    OMEGA = 1.25 # Optimal omega for SOR is problem-dependent
    X0 = np.zeros_like(B, dtype=np.float64)
    b_norm_val = np.linalg.norm(B)
    b_norm = b_norm_val if b_norm_val > ZERO_TOL else 1.0

    # --- Tridiagonal Solver (Direct Method) ---
    print("--- Tridiagonal Method ---")
    start_time = time.perf_counter()
    if is_tridiagonal_csr(A_csr):
        try:
            X_trid = solve_tridiagonal(A_csr, B)
            end_time = time.perf_counter()
            print(f"Elapsed time: {end_time - start_time:.6f}s")
            print_vector("Solution", X_trid)
            residual_norm = np.linalg.norm(_compute_residual_csr(A_csr, X_trid, B))
            print(f"Final L2 Residual Norm: {residual_norm:.2e}\n")
            write_mm_vector("X_out.mtx", X_trid)
        except ValueError as e:
            end_time = time.perf_counter()
            print(f"Elapsed time: {end_time - start_time:.6f}s")
            print(f"Solver failed: {e}\n", file=sys.stderr)
    else:
        end_time = time.perf_counter()
        print(f"Elapsed time: {end_time - start_time:.6f}s")
        print("Matrix is not tridiagonal. Skipping solver.\n")

    # --- Iterative Solvers ---
    solvers: Dict[str, Tuple[IterativeSolverFunc, Dict[str, Any]]] = {
        "Jacobi": (solve_jacobi, {"A_csr": A_csr, "D": D}),
        "Gauss-Seidel": (solve_sor_gs, {"A_csr": A_csr, "D": D, "omega": 1.0}),
        "SOR": (solve_sor_gs, {"A_csr": A_csr, "D": D, "omega": OMEGA})
    }
    for name, (solver_func, kwargs) in solvers.items():
        print(f"--- {name} Method ---")
        start_time = time.perf_counter()
        try:
            iter_count, X_sol = solver_func(B=B, X0=X0, tolerance=TOLERANCE, max_iter=MAX_ITER, b_norm=b_norm, **kwargs)
            end_time = time.perf_counter()
            print(f"Elapsed time: {end_time - start_time:.6f}s")
            if iter_count != -1:
                print(f"Converged in {iter_count} iterations.")
                print_vector("Solution", X_sol)
            else:
                print(f"Failed to converge within {MAX_ITER} iterations.", file=sys.stderr)
            res_norm = np.linalg.norm(_compute_residual_csr(A_csr, X_sol, B))
            rel_residual = res_norm / b_norm if b_norm > ZERO_TOL else res_norm
            print(f"Final Relative L2 Residual: {rel_residual:.2e}\n")
        except ValueError as e:
            end_time = time.perf_counter()
            print(f"Elapsed time: {end_time - start_time:.6f}s")
            print(f"Solver failed: {e}\n", file=sys.stderr)

    return 0

if __name__ == "__main__":
    sys.exit(main())
