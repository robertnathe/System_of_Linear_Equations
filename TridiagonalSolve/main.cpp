#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <math.h>		// required for fabs() and for sqrt()
#include <cmath>
#include <iomanip>
using namespace std;

int ST_Input(double array_A[], double B[], double X[], int nrows, int ncols);
int Tridiagonal_Input(double subdiagonal[], double diagonal[], double superdiagonal[], double array_A[], double B[], double X[], int nrows, int ncols);
int print_results (double array_A[], double B[], double X[], int nrows, int ncols);

int ST_Input(double array_A[], double B[], double X[], int nrows, int ncols)
{
  int i,j;
// Begin code here.
  double *superdiagonal;
  // Change from ncols-1 to ncols
  superdiagonal = new double[ncols];
  double *diagonal;
  diagonal = new double[ncols];
  double *subdiagonal;
  // Change from ncols-1 to ncols
  subdiagonal = new double[ncols];
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
//	std::cout << endl;
  }                                                                
// End code here.
  Tridiagonal_Input (&subdiagonal[0], &diagonal[0], &superdiagonal[0], &array_A[0], &B[0], &X[0], nrows, ncols);
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
  print_results (&array_A[0], &B[0], &X[0], nrows, ncols);
  return 0;
}

int print_results (double array_A[], double B[], double X[], int nrows, int ncols)
{
  int i, j;
  std::cout << "******************** Solve Ax = B ********************" << endl << endl;
  std::cout << "where A = " << endl;
  std::cout << std::setprecision(5);
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
	{
	  std::cout << array_A[i * ncols + j] << "   ";
	}
	std::cout << endl;
  }
  std::cout << endl;
  std::cout << "and B = " << endl;
  for (i = 0; i < ncols; i++)
    std::cout << B[i] << "   "; 
  std::cout << endl << endl;
  FILE *myFile3;
  myFile3 = fopen("C.dat","w+");
  if (myFile3 == NULL)
  {
	std::cout << "Error writing to file." << endl;
    exit(0);
  }
  fprintf(myFile3,"%%MatrixMarket_Output_vector_C.dat matrix coordinate pattern general\n");
  fprintf (myFile3,"%d %d %d\n", 1, nrows, ncols);
  std::cout << "The solution is X = " << endl;
  for (i = 0; i < ncols; i++)
  {
    fprintf(myFile3, "%d %d %lf\n", 1, i+1, X[i]);
    std::cout << X[i] << "   ";  
  }
  fclose(myFile3);
  std::cout << endl << endl;
  return 0;
}

class Input {
  public:
  
  unsigned int Input_Begin(void)  
  {
  int i,j;
  FILE *myFile;
  myFile = fopen ("sparse_array_A.dat", "r");
  int dimArray[3];
  int nrows, ncols, number_of_entries_A;
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
  nrows = dimArray[0];
  ncols = dimArray[1];
  number_of_entries_A = dimArray[2];
  double *array_A;
  array_A = new double[nrows*ncols];
  for (i = 0; i < nrows; i++)
  {
    for (j = 0; j < ncols; j++)
	{
	  array_A[i * ncols + j] = 0.0;
	}
  }
  for (i = 0; i < number_of_entries_A; i++)
  {
    fscanf (myFile, "%d,", &i_index);
    i_index--;
    fscanf (myFile, "%d,", &j_index);
    j_index--;
    fscanf (myFile, "%lf,", &value);
    array_A[i_index * ncols + j_index] = value;
  }
  fclose (myFile);
  FILE *myFile2;
  myFile2 = fopen ("sparse_array_B.dat", "r");
  int dim_B_Array[3];
  int number_of_entries_B;
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
  ncols = dimArray[1];
  number_of_entries_B = dim_B_Array[2];
  double *B;
  B = new double[ncols];
  for (j = 0; j < ncols; j++)
	{
	  B[j] = 0.0;
  }
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
  double *X;
  X = new double[ncols];
  ST_Input(&array_A[0], &B[0], &X[0], nrows, ncols);
  delete [] array_A;
  delete [] B;
  delete [] X;
  return 0;
  }
};

int main ()
{
  Input Input_One;
  Input_One.Input_Begin();
  return 0;
}
