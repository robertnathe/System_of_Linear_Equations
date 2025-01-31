Matrix Solvers and Neural Network Utilities

Overview

This repository contains implementations of several numerical methods and utilities for matrix operations and neural network training. It includes:

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

Neural Network Utilities: Includes basic functions for neural network training, including stochastic gradient descent and forward propagation.

Files

main.cpp: Contains the main function demonstrating the usage of all implemented algorithms.

Matrix_Operations_Library: A collection of utility functions for matrix and vector operations.

A.dat and B.dat: Input files containing the coefficient matrix and the right-hand side vector, respectively.

X_out.dat: Output file containing the computed solution vector.

