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
int print_results (double array[], double B[], double X_k[], int row, int col, double tolerance, int iteration, int Maximum_Number_of_Iterations);

int main ()
{
  Input_Data();
  return 0;
}
int Input_Data(void)
{
  double tolerance = 0.0001;
  int i;
  int Maximum_Number_of_Iterations = 100;
  FILE *myFile;
  myFile = fopen ("A.dat", "r");
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
  double *array_A;
  array_A = new double[row*col];
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
  int dim_B_Array[3];
  int col_B, number_of_entries_B;
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
  col_B = dim_B_Array[1];
  number_of_entries_B = dim_B_Array[2];
  double *B;
  B = new double[col_B];
  double *X_k;
  X_k = new double[col_B];
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
//  Jacobi_Solve (&array_A[0], &B[0], &X_k[0], row, col, tolerance, Maximum_Number_of_Iterations);
//  Gauss_Seidel_Solve (&array_A[0], &B[0], &X_k[0], row, col, tolerance, Maximum_Number_of_Iterations);
  SOR_Solve (&array_A[0], &B[0], &X_k[0], row, col, tolerance, Maximum_Number_of_Iterations);
  delete [] array_A;
  delete [] B;
  delete [] X_k;
  return 0;
}

int Jacobi_Solve (double array_A[], double B[], double X_k[], int row, int col, double tolerance, int Maximum_Number_of_Iterations)
{
  bool jacobi_convergence_check = true;
  int i, j, iteration = 1;
  double sum = 0;
  double *X_k1;
  X_k1 = new double[col];
  double tolerance_variable;
  double distance = 0.0, distance_squared = 0.0;
  double temporary_variable_one = 0.0;
  double temporary_variable_two = 0.0;
  double norm_variable_squared = 0.0;
  double norm_variable = 0.0;
  printf ("Prog: Jacobi_MM_Solution.cpp\n\n\n");
// X_k1[i is initialized to 0.0
  for (i = 0; i < row; i++)
  {
    X_k1[i] = 0;
  }
//  Begin Jacobi iteration modification
//  printf ("iteration %d   \n", iteration);
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
    iteration++;
// Begin tolerance evaluation function
    if (iteration >= Maximum_Number_of_Iterations)
    {
      jacobi_convergence_check = false;
    }
    for (j=0;j<col;j++)
    {
      temporary_variable_one = X_k1[j] - X_k[j];
      distance_squared += temporary_variable_one * temporary_variable_one;
      temporary_variable_two = X_k[j] * X_k[j];
      norm_variable_squared += temporary_variable_two;
    }
    distance = sqrt(distance_squared);
    norm_variable = sqrt(norm_variable_squared);
    tolerance_variable = distance / norm_variable;
    distance_squared = 0.0;
    norm_variable_squared = 0.0;
    if ( (tolerance_variable) < tolerance)
    {
      break;
    }
    if (norm_variable > 100.0)
    {
      break;
    }
// End tolerance evaluation function
  }
  print_results (&array_A[0], &B[0], &X_k1[0], row, col, tolerance, iteration, Maximum_Number_of_Iterations);
  delete [] X_k1;
  return 0;
}

int Gauss_Seidel_Solve (double array_A[], double B[], double X_k[], int row, int col, double tolerance, int Maximum_Number_of_Iterations)
{
  bool Gauss_Seidel_convergence_check = true;
  int i, j, iteration = 1;
  double sum = 0;
  double *X_k1;
  X_k1 = new double[col];
  double tolerance_variable;
  double distance = 0.0, distance_squared = 0.0;
  double temporary_variable_one = 0.0;
  double temporary_variable_two = 0.0;
  double norm_variable_squared = 0.0;
  double norm_variable = 0.0;
  printf ("Prog: Gauss_Seidel_MM_Solution.cpp\n\n\n");
// X_k[i] is initialized to 0.0
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
//  printf ("Gauss_Seidel_iteration %d   \n", iteration);
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
    iteration++;
    // Begin tolerance evaluation function
    if (iteration >= Maximum_Number_of_Iterations)
    {
      Gauss_Seidel_convergence_check = false;
    }
    for (j=0;j<col;j++)
    {
      temporary_variable_one = X_k1[j] - X_k[j];
      distance_squared += temporary_variable_one * temporary_variable_one;
      temporary_variable_two = X_k[j] * X_k[j];
      norm_variable_squared += temporary_variable_two;
    }
    distance = sqrt(distance_squared);
    norm_variable = sqrt(norm_variable_squared);
    tolerance_variable = distance / norm_variable;
    distance_squared = 0.0;
    norm_variable_squared = 0.0;
    if ( (tolerance_variable) < tolerance)
    {
      break;
    }
    if (norm_variable > 100.0)
    {
      break;
    }
// End tolerance evaluation function
    if (iteration >= Maximum_Number_of_Iterations)
    {
      Gauss_Seidel_convergence_check = false;
    }
  }
  print_results (&array_A[0], &B[0], &X_k1[0], row, col, tolerance, iteration, Maximum_Number_of_Iterations);
  delete [] X_k1;
  return 0;
}

int SOR_Solve (double array_A[], double B[], double X_k[], int row, int col, double tolerance, int Maximum_Number_of_Iterations)
{
  double SOR_omega = 1.25;
  bool SOR_convergence_check = true;
  int i, j, iteration = 1;
  double sum = 0;
  double *X_k1;
  X_k1 = new double[col];
  double tolerance_variable;
  double distance = 0.0, distance_squared = 0.0;
  double temporary_variable_one = 0.0;
  double temporary_variable_two = 0.0;
  double norm_variable_squared = 0.0;
  double norm_variable = 0.0;
  printf ("Prog: SOR_MM_Solution.cpp\n\n\n");
// Initialize X_k[i] to 1.0
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
//  printf ("iteration %d   \n", iteration);
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
    iteration++;
    // Begin tolerance evaluation function
    if (iteration >= Maximum_Number_of_Iterations)
    {
      SOR_convergence_check = false;
    }
    for (j=0;j<col;j++)
    {
      temporary_variable_one = X_k1[j] - X_k[j];
      distance_squared += temporary_variable_one * temporary_variable_one;
      temporary_variable_two = X_k[j] * X_k[j];
      norm_variable_squared += temporary_variable_two;
    }
    distance = sqrt(distance_squared);
    norm_variable = sqrt(norm_variable_squared);
    tolerance_variable = distance / norm_variable;
    distance_squared = 0.0;
    norm_variable_squared = 0.0;
//    printf("tolerance_variable %6.4f \n\n", tolerance_variable);
    if ( (tolerance_variable) < tolerance)
    {
      break;
    }
    if (norm_variable > 100.0)
    {
      break;
    }
// End tolerance evaluation function
    if (iteration >= Maximum_Number_of_Iterations)
    {
      SOR_convergence_check = false;
    }
  }
  print_results (&array_A[0], &B[0], &X_k1[0], row, col, tolerance, iteration, Maximum_Number_of_Iterations);
  delete [] X_k1;
  return 0;
}

int print_results (double array[], double x0[], double x[], int row, int col, double Tolerance, int iteration, int Maximum_Number_of_Iterations)
{
  int i, j;
  printf ("******************** Solve Ax = B ********************\n\n");
  printf ("where A = \n");
  for (i = 0; i < row; i++)
    {
      for (j = 0; j < col; j++)
	{
	  printf ("%6.4f   ", array[i * row + j]);
	}
      printf ("\n");
    }
  printf ("\n");
  printf ("and B = \n");
  for (i = 0; i < col; i++)
    printf ("%6.4f   ", x0[i]);
  printf ("\n\n");
  FILE *myFile3;
  myFile3 = fopen("C.dat","w+");
  if (myFile3 == NULL)
  {
    printf("Error writing to file.\n");
    exit(0);
  }
  fprintf(myFile3,"%%MatrixMarket_Output_vector_C.dat matrix coordinate pattern general\n");
  fprintf (myFile3,"%d %d %d\n", 1, row, col);
  printf ("The solution is X = \n");
  for (i = 0; i < col; i++)
  {
    fprintf(myFile3, "%d %d %lf\n", 1, i+1, x[i]);
    printf ("%6.4f    ", x[i]);
  }
  fclose(myFile3);
  printf ("\n\n");
  printf ("The ");
  printf ("%d", iteration);
  printf ("th iteration resulted in a vector with the allowed tolerance.\n");
  printf ("\n");
  printf ("The tolerance is ");
  printf ("%6.4f \n", Tolerance);
  printf ("The Maximum_Number_of_Iterations is %d\n", Maximum_Number_of_Iterations);
  return 0;
}
