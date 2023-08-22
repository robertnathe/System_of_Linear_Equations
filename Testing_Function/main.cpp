#include "MyIncludeFile.h"
#include <chrono>
template<typename T>
using matrix = std::vector< std::vector<T> >;
// Need to implement these functions.

int Full_Pivoting_Solution(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols);
int Rook_Pivoting_Solution(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols);
int Row_Pivoting_Solution(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols);
int Column_Pivoting_Solution(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols);

void PrintVector2D(const vector<double> A[], unsigned int nrows);
void PrintVector(const double TempVector[], unsigned int nrows);

int Full_Pivoting_Solution(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols)
{
    unsigned int i {0}, j {0}, k {0};
    unsigned int pi,pj;
    double max;
    double TempVar {0.0};
    for (k=0;k<nrows;k++){
        pi=0;pj=0;max=0.0;
        //find pivot in submatrix a(k:n,k:n)
        for (i=k;i<nrows;i++) {
            for (j=k;j<nrows;j++) {
                if (fabs(A[i][j])>max){
                    max = fabs(A[i][j]);
                    pi=i;
                    pj=j;
                }
            }
        }
        if ( (A[pi][pj]) >= (A[k][k]) ) {
        //Swap Row
        for (j=0;j<nrows;j++){
          std::swap(A[pi][j], A[k][j]);
          TempVar = B[pi];
          B[pi] = B[k];
          B[k] = TempVar;
        }
        //Swap Col
        for (i=0;i<ncols;i++){
          std::swap(A[i][pj], A[i][k]);
          TempVar = B[pj];
          B[pj] = B[k];
          B[k] = TempVar;          
        }
    }
    //END PIVOT
  }
  PrintVector(&B[0], ncols);
  PrintVector2D(&A[0], ncols);
  return 0;
}

int Rook_Pivoting_Solution(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols)
{
  // Find the best pivot
  unsigned int i {0}, p {0}, k {0};
  double TempVar {0.0};
  double maxPivotRow {0.0}, maxPivotColumn {0.0};
  for (k = 0; k < ncols; k++) 
  {
	maxPivotRow = 0.0;  
    maxPivotColumn = 0.0;
    for (i = k; i < nrows; i++)
    {
      if ((abs(A[i][k]) > maxPivotRow)) {
	    maxPivotRow = abs(A[i][k]);
	    p = i;
      }
    }
    for (i = k; i < ncols; i++)
    {
      if ((abs(A[k][i]) > maxPivotColumn)) {
	    maxPivotColumn = abs(A[k][i]);
	    p = i;
      }
    }
    if (maxPivotRow >= maxPivotColumn) {
	  // Swap rows k and p
      for (i = k; i < nrows; i++) {
        std::swap(A[p][i], A[k][i]);
        TempVar = B[p];
        B[p] = B[k];
        B[k] = TempVar;
      }
	}
	else {
	  // Swap rows k and p
      for (i = k; i < ncols; i++) {
        std::swap(A[i][p], A[i][k]);
        TempVar = B[p];
        B[p] = B[k];
        B[k] = TempVar;
      }
    }
  }
  PrintVector(&B[0], ncols);
  PrintVector2D(&A[0], ncols);
  return 0;
}

int Row_Pivoting_Solution(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols)
{
  // Find the best pivot
  unsigned int i {0}, p {0}, k {0};
  double TempVar {0.0}, maxPivotRow {0.0};
  for (k = 0; k < nrows; k++) 
  {
    maxPivotRow = 0.0;
    for (i = k; i < nrows; i++)
    {
      if ((abs(A[i][k]) > maxPivotRow)) {
	    maxPivotRow = abs(A[i][k]);
	    p = i;
      }
    }
    // Swap rows k and p
    for (i = k; i < nrows; i++) {
      std::swap(A[p][i], A[k][i]);
      TempVar = B[p];
      B[p] = B[k];
      B[k] = TempVar;
    }
  PrintVector(&B[0], nrows);
  PrintVector2D(&A[0], nrows);
  }
  return 0;
}

int Column_Pivoting_Solution(vector<double> A[], double X[], double B[], unsigned int nrows, unsigned int ncols)
{
  // Find the best pivot
  unsigned int i {0}, p {0}, k {0};
  double TempVar {0.0}, maxPivotColumn {0.0};
  (void) TempVar;
  for (k = 0; k < ncols; k++) 
  {
    maxPivotColumn = 0.0;
    for (i = k; i < ncols; i++)
    {
      if ((abs(A[k][i]) > maxPivotColumn)) {
	    maxPivotColumn = abs(A[k][i]);
	    p = i;
      }
    }
    // Swap rows k and p
    for (i = k; i < ncols; i++) {
      std::swap(A[i][p], A[i][k]);
      TempVar = B[p];
      B[p] = B[k];
      B[k] = TempVar;
    }
  }
  PrintVector(&B[0], ncols);
  PrintVector2D(&A[0], ncols);
  return 0;
}

/*
Iterate over vector of vectors and for each of the
nested vector print its contents
*/
void PrintVector2D(const vector<double> TempMatrix[], unsigned int nrows)
{
  std::cout << "Displaying the 2D vector:" << endl;
/*
Iterate over vector of vectors and for each of the
nested vector print its contents
*/
// Displaying the 2D vector
  std::cout << std::setprecision(5);
  for (unsigned int i = 0; i < nrows; i++) {
    for (
      auto it = TempMatrix[i].begin();
        it != TempMatrix[i].end(); it++)
        cout << *it << " ";
      cout << endl;
    }
}
void PrintVector(const double TempVector[], unsigned int nrows){
  std::cout << "Displaying the vector: " << endl;
  std::cout << std::setprecision(5);
  for (unsigned int i = 0; i < nrows; i++)
    std::cout << TempVector[i] << "   ";
  std::cout << endl;
}

class Input {
  public:
  unsigned int Input_Begin(void)
  {
    unsigned int i,j, max_iter {1000};
    (void) max_iter;
    double tolerance {0.00001};
    (void) tolerance;
    double eigenvalue {0.0};
    (void) eigenvalue;
    double Scalar {3.0};
    (void) Scalar;
    FILE *myFile;
    myFile = fopen ("A.dat", "r");
    unsigned int dimArray[3];
    unsigned int nrows, number_of_entries_A;
    unsigned int ncols;
    unsigned int i_index, j_index;
    double value, elem {0.0};
    while (myFile == NULL)
    {
	  std::cout << "Error Reading File" << endl;
      exit (0);
    }
    fscanf (myFile, "%*s %*s %*s %*s %*s");
    for (i = 0; i < 3; i++)
    {
    fscanf (myFile, "%u,", &dimArray[i]);
    }
    nrows = dimArray[0];
    ncols = dimArray[1];
    number_of_entries_A = dimArray[2];
    vector <double> array_A;
    for (i = 0; i < number_of_entries_A; i++)
    {
      fscanf (myFile, "%u,", &i_index);
      i_index--;
      fscanf (myFile, "%u,", &j_index);
      j_index--;
      fscanf (myFile, "%lf,", &value);
//    Change program to use the single index array_A
      array_A.push_back(value);
    }
    fclose (myFile);
    FILE *myFile2;
    myFile2 = fopen ("B.dat", "r");
    unsigned int dim_B_Array[3];
    unsigned int number_of_entries_B;
    unsigned int col_B;
    (void) col_B;
    while (myFile2 == NULL)
    {
	  std::cout << "Error Reading File" << endl;
      exit (0);
    }
    fscanf (myFile2, "%*s %*s %*s %*s %*s");
    for (i = 0; i < 3; i++)
    {
      fscanf (myFile2, "%u,", &dim_B_Array[i]);
    }
    col_B = dim_B_Array[1];
    number_of_entries_B = dim_B_Array[2];
    vector <double> B;
    vector <double> X;
    vector <double> MyVector1;
    vector <double> MyVector2;
    vector <double> MyVector3;
    vector <int> P;
    for (i = 0; i < number_of_entries_B; i++)
    {
      fscanf (myFile2, "%u,", &i_index);
      i_index--;
      fscanf (myFile2, "%u,", &j_index);
      j_index--;
      fscanf (myFile2, "%lf,", &value);
      B.push_back(value);
      X.push_back(0.0);
      MyVector1.push_back(10.0+i);
      MyVector2.push_back(10.0-i);
      MyVector3.push_back(0.0);
      P.push_back(0);
    }
    fclose (myFile2);
    // Initializing the vector of vectors
//    vector<vector<double> > A;
    matrix<double> A, AMat, BMat, CMat;
    // Inserting elements into vector
    for (i = 0; i < nrows; i++) {
      // Vector to store column elements
      vector<double> Row1, Row2, Row3, Row4;
        for (j = 0; j < ncols; j++) {
		  elem = array_A[i*nrows+j];
          Row1.push_back(elem);
          Row2.push_back(elem);
          Row3.push_back(elem);
          Row4.push_back(elem);
            if (j == (ncols)) {
				Row1.push_back(elem);
				Row2.push_back(elem);
				Row3.push_back(elem);
				Row4.push_back(elem);				
			}
        }
        // Pushing back above 1D vector
        // to create the 2D vector
        A.push_back(Row1);
        AMat.push_back(Row2);
        BMat.push_back(Row3);
        CMat.push_back(Row4);
    }
    matrix<double> Augmented;
    // Elements to insert in column
    // Inserting elements into vector
    for (i = 0; i < nrows; i++) {
    // Vector to store column elements
      vector<double> Augment;
      for (j = 0; j < ncols+1; j++) {
        Augment.push_back(elem);
      }
      // Pushing back above 1D vector
      // to create the 2D vector
      Augmented.push_back(Augment);
    }
    for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < ncols+1; j++)
	  {
	    if (j != (ncols)) {
	      Augmented[i][j] = A[i][j];
	    }
	    if (j == (ncols)) {
	      Augmented[i][ncols] = B[i];
	    }
      }
    }  
    
//    Full_Pivoting_Solution(&A[0], &X[0], &B[0], nrows, ncols);
//    Rook_Pivoting_Solution(&A[0], &X[0], &B[0], nrows, ncols);
//    Row_Pivoting_Solution(&A[0], &X[0], &B[0], nrows, ncols);
//    Column_Pivoting_Solution(&A[0], &X[0], &B[0], nrows, ncols);
//    Testing_General_Purpose_Parameters (&A[0], &X[0], &B[0], nrows, ncols);
//    Testing_Arithmetic_Matrix_Parameters (&A[0], &X[0], &B[0], &AMat[0], &BMat[0], &CMat[0], &MyVector1[0], &MyVector2[0], &MyVector3[0], Scalar, nrows, ncols); 

    array_A.clear();
    array_A.shrink_to_fit();
    B.clear();
    B.shrink_to_fit();
    X.clear();
    X.shrink_to_fit();
    A.clear();
    A.shrink_to_fit();
    AMat.clear();
    AMat.shrink_to_fit();
    BMat.clear();
    BMat.shrink_to_fit();
    CMat.clear();
    CMat.shrink_to_fit();
    Augmented.clear();
    Augmented.shrink_to_fit();
    MyVector1.clear();
    MyVector1.shrink_to_fit();
    MyVector2.clear();
    MyVector2.shrink_to_fit();
    MyVector3.clear();
    MyVector3.shrink_to_fit();
    P.clear();
    P.shrink_to_fit();
    
//    std::cin.clear(); // reset any error flags
//    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//    ignore any characters in the input buffer until we find an enter character
//    std::cin.get(); // get one more char from the user
    return 0;
  }
};

int main ()
{
  // Using time point and system_clock
  std::chrono::time_point<std::chrono::system_clock> start, end;
 
  start = std::chrono::system_clock::now();	
	
  Input Input_One;
  Input_One.Input_Begin();
  
  end = std::chrono::system_clock::now();
 
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);
 
  std::cout << "finished computation at " << std::ctime(&end_time)
            << "elapsed time: " << elapsed_seconds.count() << "s\n";
            
  return 0;
}
