#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <math.h>		// required for fabs() and for sqrt()
#include <cmath>
using namespace std;
int Input_Data (void);
int Tridiagonal_Input(double subdiagonal[], double diagonal[], double superdiagonal[], double array_A[], double B[], double X[], int nrows, int ncols);
int print_results (double subdiagonal[], double diagonal[], double superdiagonal[], double array_A[], double B[], double X[], int nrows, int ncols);

int main ()
{
  Input_Data();
  return 0;
}
int Input_Data(void)
{
  int i;
  FILE *myFile;
  myFile = fopen ("A.dat", "r");
  int dimArray[3];
  int nrows, ncols, number_of_entries_A;
  int i_index, j_index;
  double value;
  while (myFile == NULL)
  {
    printf ("Error Reading File\n");
    exit (0);
  }
  fscanf (myFile, "%*s %*s %*s %*s %*s");
  for (i = 0; i < 3; i++)
  {
  fscanf (myFile, "%d,", &dimArray[i]);
  }
  nrows = dimArray[0];
  ncols = dimArray[1];
  number_of_entries_A = dimArray[2];
  double *array_A;
  array_A = new double[nrows*ncols];
  for (i = 0; i < number_of_entries_A; i++)
  {
    fscanf (myFile, "%d,", &i_index);
    i_index--;
    fscanf (myFile, "%d,", &j_index);
    j_index--;
    fscanf (myFile, "%lf,", &value);
//    Change program to use the single index array_A
    array_A[i_index * ncols + j_index] = value;
  }
  fclose (myFile);
  FILE *myFile2;
  myFile2 = fopen ("B.dat", "r");
  int dim_B_Array[3];
  int col_B, number_of_entries_B;
  while (myFile2 == NULL)
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
  double *X;
  X = new double[col_B];
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
  int j;
  double *superdiagonal;
  superdiagonal = new double[col_B-1];
  double *diagonal;
  diagonal = new double[col_B];
  double *subdiagonal;
  subdiagonal = new double[col_B-1];
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
	{
	  if (i == j)
	  {
	    diagonal[i] = array_A[i * ncols + j];
      }
      else if  (i == j-1)
      {
	    superdiagonal[i] = array_A[i*ncols+j];
	  }
	  else if (i == j+1)
	  {
	    subdiagonal[i] = array_A[i*ncols+j];
	  }
	}
	printf("\n");
  }
  Tridiagonal_Input (&subdiagonal[0], &diagonal[0], &superdiagonal[0], &array_A[0], &B[0], &X[0], nrows, ncols);
  delete [] array_A;
  delete [] B;
  delete [] X;
  delete [] subdiagonal;
  delete [] diagonal;
  delete [] superdiagonal;
  return 0;
}

int Tridiagonal_Input(double subdiagonal[], double diagonal[], double superdiagonal[], double array_A[], double B[], double X[], int nrows, int ncols)
{
  int i;
  int err = 0;
  double m;
//         Check that all diagonal elements are nonzero.
//         If a diagonal element is zero then U is singular, so return
//         signalling an error.
  for (i = 0; i < ncols; i++) if (diagonal[i] == 0.0) err = -1;
  //forward elimination
  for (i=1;i<=ncols-1;i++)
  {
    if(i==0){
      m = subdiagonal[i]/diagonal[i-1];
    }
    if (i != 0){
      m = (subdiagonal[i])/diagonal[i-1];
    }
    if (i == 0)
    {
      diagonal[i] = diagonal[i] - (m * superdiagonal[i-1]);
    }
    if (i != 0)
    {
      diagonal[i] = diagonal[i] - (m * superdiagonal[i-1]);
    }
    B[i] = B[i] - (m * B[i-1]);
  }
  //backward substitution
  X[ncols-1]=B[ncols]/diagonal[ncols];
  for(i=ncols-1;i>=0;i--)
  {
    X[i]=(B[i]-superdiagonal[i]*X[i+1])/diagonal[i];
  }
  if (err == -1)
    printf("Matrix A fails the LU decomposition\n");
  print_results (&subdiagonal[0], &diagonal[0], &superdiagonal[0], &array_A[0], &B[0], &X[0], nrows, ncols);
  return 0;
}

int print_results (double subdiagonal[], double diagonal[], double superdiagonal[], double array_A[], double B[], double X[], int nrows, int ncols)
{
  int i, j;
  printf("\n");
  printf ("******************** Solve Ax = B ********************\n\n");
  printf ("where A = \n");
  for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < ncols; j++)
	{
	  printf ("%6.4f   ", array_A[i * ncols + j]);
	}
      printf ("\n");
    }
  printf ("\n");
  printf ("and B = \n");
  for (i = 0; i < ncols; i++)
    printf ("%6.4f   ", B[i]);
  printf ("\n\n");
  FILE *myFile3;
  myFile3 = fopen("X.dat","w+");
  if (myFile3 == NULL)
  {
    printf("Error writing to file.\n");
    exit(0);
  }
  fprintf(myFile3,"%%MatrixMarket_Output_vector_X.dat matrix coordinate pattern general\n");
  fprintf (myFile3,"%d %d %d\n", 1, nrows, ncols);
  printf ("The solution is X = \n");
  for (i = 0; i < ncols; i++)
  {
    fprintf(myFile3, "%d %d %lf\n", 1, i+1, X[i]);
    printf ("%6.4f    ", X[i]);
  }
  fclose(myFile3);
  printf ("\n\n");
  return 0;
}
