Matrix Solvers

Overview

Several numerical methods and utilities are implemented for matrix operations.

Direct Solvers:

Tridiagonal Thomas Algorithm: Efficiently solves tridiagonal systems.

Iterative Solvers:

Jacobi Method (L1 and L2 Norms): A simple iterative method for solving linear systems.

Gauss-Seidel Method (L1 Norm): An improved iterative method over Jacobi.

Successive Overrelaxation (SOR) Method (L1 Norm): An optimized version of Gauss-Seidel with adjustable relaxation parameter.

Features

Input/Output Formats: Supports the Matrix Market file format for reading and writing matrices and vectors.

Performance Measurement: Includes timing utilities to measure the execution time of each algorithm.

Flexibility: Allows for easy modification and extension of algorithms.

Files

main.cpp: Contains the main function demonstrating the usage of all implemented algorithms.

Matrix_Operations_Library: A collection of utility functions for matrix and vector operations.

Input files containing the coefficient matrix and the right-hand side vector, respectively.

Output file containing the computed solution vector.
