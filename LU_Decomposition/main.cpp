#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <math.h>		// required for fabs() and for sqrt()
#include <cmath>
using namespace std;
int Input_Data (void);
int Jacobi_Solve(double array_A[], double B[], double X_k[], int row,
		      int col, double Tolerance, int Maximum_Number_of_Iterations);
int Gauss_Seidel_Solve(double array_A[], double B[], double X_k[], int row,
		      int col, double Tolerance, int Maximum_Number_of_Iterations);
int SOR_Solve(double array_A[], double B[], double X_k[], int row,
		      int col, double Tolerance, int Maximum_Number_of_Iterations);
int print_results (double array[], double B[], double X_k[], int row, int col, double tolerance, int Maximum_Number_of_Iterations);

int main ()
{
  Input_Data();
  return 0;
}

int Input_Data(void)
{
  double tolerance = 0.001;
  int i, j, kk;
  int Maximum_Number_of_Iterations = 5;
  FILE *myFile;
  myFile = fopen ("A.dat", "r");
  //read file into array
  int dimArray[3];
  int row, col, number_of_entries_A;
  int i_index, j_index;
  double value;
  if (myFile == NULL)
    {
      printf ("Error Reading File\n");
      exit (0);
    }
  fscanf (myFile, "%*s %*s %*s %*s %*s");
  for (i = 0; i < 3; i++)
    {
      fscanf (myFile, "%d,", &dimArray[i]);
    }
  row = dimArray[0];
  col = dimArray[1];
  number_of_entries_A = dimArray[2];
  kk = row * col;
  double array_A[kk];
  for (i = 0; i < number_of_entries_A; i++)
    {
      fscanf (myFile, "%d,", &i_index);
      i_index--;
      fscanf (myFile, "%d,", &j_index);
      j_index--;
      fscanf (myFile, "%lf,", &value);
//    Change program to use the single index array_A
      array_A[i_index * row + j_index] = value;
    }
  fclose (myFile);
  FILE *myFile2;
  myFile2 = fopen ("B.dat", "r");
  //read file into array
  int dim_B_Array[3];
  int row_B, col_B, number_of_entries_B;
  if (myFile2 == NULL)
    {
      printf ("Error Reading File\n");
      exit (0);
    }
  fscanf (myFile2, "%*s %*s %*s %*s %*s");
  for (i = 0; i < 3; i++)
    {
      fscanf (myFile2, "%d,", &dim_B_Array[i]);
    }
  row_B = dim_B_Array[0];
  col_B = dim_B_Array[1];
  number_of_entries_B = dim_B_Array[2];
  double B[col_B], X_k[col_B];
  for (i = 0; i < number_of_entries_B; i++)
    {
      fscanf (myFile2, "%d,", &i_index);
      i_index--;
      fscanf (myFile2, "%d,", &j_index);
      j_index--;
      fscanf (myFile2, "%lf,", &value);
      B[j_index] = value;
    }
  fclose (myFile2);
  Jacobi_Solve (&array_A[0], &B[0], &X_k[0], row, col, tolerance, Maximum_Number_of_Iterations);
//  Gauss_Seidel_Solve (&array_A[0], &B[0], &X_k[0], row, col, tolerance, Maximum_Number_of_Iterations);
//  SOR_Solve (&array_A[0], &B[0], &X_k[0], row, col, tolerance, Maximum_Number_of_Iterations);
  return 0;
}

int Jacobi_Solve (double array_A[], double B[], double X_k[], int row, int col, double tolerance, int Maximum_Number_of_Iterations)
{
  bool jacobi_convergence_check = true;
  int i, j, jacobi_iteration = 1;
  double sum = 0;
  double X_k1[col];
  printf ("Prog: Jacobi_MM_Solution.cpp\n\n\n");
  for (i = 0; i < row; i++)
  {
    X_k1[i] = 0;
  }
//  Begin Gauss-Seidel iteration modification
//  for (i = 0; i<row; i++)
//  {
//    sum = 0;
//    for (j=0; j<col; j++)
//    {
//      if (j < i)
//      {
//        sum = sum + array_A[i*row+j] * X_k1[j];
//      }
//      if (j > i)
//      {
//        sum = sum + array_A[i*row+j] * X_k[j];
//      }
//    }
//    X_k1[i] = -sum / array_A[i*row+i] + B[i] / array_A[i*row+i];
//  }
//  End Gauss-Seidel iteration modification

//  Begin Jacobi iteration modification
  printf ("jacobi_iteration %d   \n", jacobi_iteration);
  for (i = 0; i<row; i++)
  {
    sum = 0;
    for (j=0; j<col; j++)
    {
      if (i != j)
      {
        sum = sum + array_A[i*row+j] * X_k[j];
      }
    }
    X_k1[i] = -sum / array_A[i*row+i] + B[i] / array_A[i*row+i];

    printf ("%6.3f   ", X_k[i]);
    printf ("%6.3f   \n", X_k1[i]);
  }
//  End Jacobi iteration modification
  while (jacobi_convergence_check == true)
  {
    for (i = 0; i<row; i++)
    {
      X_k[i] = X_k1[i];
    }
    for (i = 0; i<row; i++)
    {
      sum = 0;
      for (j=0; j<col; j++)
      {
        if (i != j)
        {
          sum = sum + array_A[i*row+j] * X_k[j];
        }
      }
      X_k1[i] = -sum / array_A[i*row+i] + B[i] / array_A[i*row+i];
    }
    jacobi_iteration++;
    if (jacobi_iteration >= Maximum_Number_of_Iterations)
    {
      jacobi_convergence_check = false;
    }
    printf ("jacobi_iteration %d   \n", jacobi_iteration);
    for (i = 0; i<row; i++)
    {
      printf ("%6.3f   ", X_k[i]);
      printf ("%6.3f   \n", X_k1[i]);
    }
  }
//  print_results (&array_A[0], &B[0], &X_k[0], row, col, tolerance, Maximum_Number_of_Iterations);
  return 0;
}

int Gauss_Seidel_Solve (double array_A[], double B[], double X_k[], int row, int col, double tolerance, int Maximum_Number_of_Iterations)
{
  bool Gauss_Seidel_convergence_check = true;
  int i, j, Gauss_Seidel_iteration = 1;
  double sum = 0;
  double X_k1[col];
  printf ("Prog: Gauss_Seidel_MM_Solution.cpp\n\n\n");
  for (i = 0; i < row; i++)
  {
    X_k1[i] = 0;
  }
//  Begin Gauss-Seidel iteration modification
  for (i = 0; i<row; i++)
  {
    sum = 0;
    for (j=0; j<col; j++)
    {
      if (j < i)
      {
        sum = sum + array_A[i*row+j] * X_k1[j];
      }
      if (j > i)
      {
        sum = sum + array_A[i*row+j] * X_k[j];
      }
    }
    X_k1[i] = -sum / array_A[i*row+i] + B[i] / array_A[i*row+i];
  }
//  End Gauss-Seidel iteration modification

//  Begin Jacobi iteration modification
  printf ("Gauss_Seidel_iteration %d   \n", Gauss_Seidel_iteration);
  for (i = 0; i<row; i++)
  {
    sum = 0;
    for (j=0; j<col; j++)
    {
      if (j < i)
      {
        sum = sum + array_A[i*row+j] * X_k1[j];
      }
      if (j > i)
      {
        sum = sum + array_A[i*row+j] * X_k[j];
      }
    }
    X_k1[i] = -sum / array_A[i*row+i] + B[i] / array_A[i*row+i];
  }
  for (i = 0; i<row; i++)
  {
    printf ("%6.3f   ", X_k[i]);
    printf ("%6.3f   \n", X_k1[i]);
  }
//  End Jacobi iteration modification
  while (Gauss_Seidel_convergence_check == true)
  {
    for (i = 0; i<row; i++)
    {
      X_k[i] = X_k1[i];
    }
    for (i = 0; i<row; i++)
    {
      sum = 0;
      for (j=0; j<col; j++)
      {
        if (j < i)
        {
          sum = sum + array_A[i*row+j] * X_k1[j];
        }
        if (j > i)
        {
          sum = sum + array_A[i*row+j] * X_k[j];
        }
      }
      X_k1[i] = -sum / array_A[i*row+i] + B[i] / array_A[i*row+i];
    }
    Gauss_Seidel_iteration++;
    if (Gauss_Seidel_iteration >= Maximum_Number_of_Iterations)
    {
      Gauss_Seidel_convergence_check = false;
    }
    printf ("Gauss_Seidel_iteration %d   \n", Gauss_Seidel_iteration);
    for (i = 0; i<row; i++)
    {
      printf ("%6.3f   ", X_k[i]);
      printf ("%6.3f   \n", X_k1[i]);
    }
  }
//  print_results (&array_A[0], &B[0], &X_k[0], row, col, tolerance, Maximum_Number_of_Iterations);
  return 0;
}

int SOR_Solve (double array_A[], double B[], double X_k[], int row, int col, double tolerance, int Maximum_Number_of_Iterations)
{
  double SOR_omega = 1.25;
  bool SOR_convergence_check = true;
  int i, j, SOR_iteration = 1;
  double sum = 0;
  double X_k1[col];
  printf ("Prog: SOR_MM_Solution.cpp\n\n\n");
  for (i = 0; i < row; i++)
  {
    X_k[i] = 1.0;
    X_k1[i] = 0.0;
  }
//  Begin SOR iteration modification
  for (i = 0; i<row; i++)
  {
    sum = 0;
    for (j=0; j<col; j++)
    {
      if (j < i)
      {
        sum = sum + array_A[i*row+j] * X_k1[j];
      }
      if (j > i)
      {
        sum = sum + array_A[i*row+j] * X_k[j];
      }
    }
    // Begin SOR formula
    X_k1[i] = X_k[i] + SOR_omega * ( (B[i] - sum) / (array_A[i*row+i]) - X_k[i] );
    // End SOR formula
//    X_k1[i] = -sum / array_A[i*row+i] + B[i] / array_A[i*row+i];
  }
//  End Gauss-Seidel iteration modification

//  Begin Jacobi iteration modification
  printf ("SOR_iteration %d   \n", SOR_iteration);
  for (i = 0; i<row; i++)
  {
    sum = 0;
    for (j=0; j<col; j++)
    {
      if (j < i)
      {
        sum = sum + array_A[i*row+j] * X_k1[j];
      }
      if (j > i)
      {
        sum = sum + array_A[i*row+j] * X_k[j];
      }
    }
    // Begin SOR formula
    X_k1[i] = X_k[i] + SOR_omega * ( (B[i] - sum) / (array_A[i*row+i]) - X_k[i] );
    // End SOR formula
//    X_k1[i] = -sum / array_A[i*row+i] + B[i] / array_A[i*row+i];

  }
  for (i = 0; i<row; i++)
  {
    printf ("%6.3f   ", X_k[i]);
    printf ("%6.3f   \n", X_k1[i]);
  }
//  End Jacobi iteration modification
  while (SOR_convergence_check == true)
  {
    for (i = 0; i<row; i++)
    {
      X_k[i] = X_k1[i];
    }
    for (i = 0; i<row; i++)
    {
      sum = 0;
      for (j=0; j<col; j++)
      {
        if (j < i)
        {
          sum = sum + array_A[i*row+j] * X_k1[j];
        }
        if (j > i)
        {
          sum = sum + array_A[i*row+j] * X_k[j];
        }
      }
    // Begin SOR formula
    X_k1[i] = X_k[i] + SOR_omega * ( (B[i] - sum) / (array_A[i*row+i]) - X_k[i] );
    // End SOR formula
//      X_k1[i] = -sum / array_A[i*row+i] + B[i] / array_A[i*row+i];
    }
    SOR_iteration++;
    if (SOR_iteration >= Maximum_Number_of_Iterations)
    {
      SOR_convergence_check = false;
    }
    printf ("SOR_iteration %d   \n", SOR_iteration);
    for (i = 0; i<row; i++)
    {
      printf ("%6.3f   ", X_k[i]);
      printf ("%6.3f   \n", X_k1[i]);
    }
  }
//  print_results (&array_A[0], &B[0], &X_k[0], row, col, tolerance, Maximum_Number_of_Iterations);
  return 0;
}

int print_results (double array[], double x0[], double x_k[], int row, int col, double tolerance, int Maximum_Number_of_Iterations)
{
  int i, j;
  printf ("******************** Solve Ax = B ********************\n\n");
  printf ("where A = \n");
  for (i = 0; i < row; i++)
    {
      for (j = 0; j < col; j++)
	{
	  printf ("%6.3f   ", array[i * row + j]);
	}
      printf ("\n");
    }
  printf ("\n");
  printf ("and B = \n");
  for (i = 0; i < col; i++)
    printf ("%6.3f   ", x0[i]);
  printf ("\n\n");
  printf ("The solution is X = \n");
  for (i = 0; i < col; i++)
    {
      printf ("%6.3f   ", x_k[i]);
      // Assigning values to the arrays X and pivot to fix warning message
    }
  printf ("\n\n");
  printf ("\n\n\n");
  printf ("Tolerance: %6.3f\n", tolerance);
  printf ("Maximum_Number_of_Iterations: %d\n", Maximum_Number_of_Iterations);
  return 0;
}				// static void
