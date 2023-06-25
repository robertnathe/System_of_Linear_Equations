#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
using namespace std; 
#define a(i,j) a[(i)*nrows+(j)]
int SOR(vector<double> A[], double B[], double X[], unsigned int nrows, unsigned int ncols, double tolerance, unsigned int Max_Iter, bool convergence_check);
int Gauss_Seidel(vector<double> A[], double B[], double X[], unsigned int nrows, unsigned int ncols, double tolerance, unsigned int Max_Iter, bool convergence_check);
int Jacobi(vector<double> A[], double B[], double X[], unsigned int nrows, unsigned int ncols, double tolerance, unsigned int Max_Iter, double convergence_check);
int Output_Data(vector<double> A[], double B[], double X, unsigned int nrows, unsigned int ncols);
void PrintVector2D(const vector<double> A[], unsigned int nrows);
void PrintVector(const double B[], unsigned int nrows);

template<typename T>
using matrix = std::vector< std::vector<T> >;

int SOR (vector<double> A[], double B[], double X[], unsigned int nrows, unsigned int ncols, double tolerance, unsigned int Max_Iter, bool convergence_check)
{  
  std::cout << endl << "Implementing the SOR method: " << endl;
  unsigned int i {0}, j {0}, Iter {1};
  double elem {1.0}, omega {1.25}, distance {0.0}, distance_max {0.0};
  vector<double> XPrevious;
// Change omega to compare the Gauss-Seidel and SOR methods.
//  std::cout << "Please Enter the omega factor for the SOR method." << endl;
//  std::cin >> omega;
  for (j=0; j<ncols; j++) {
    XPrevious.push_back(elem);
  }
  std::cout << std::setprecision(5);
  while(convergence_check != true) {    
    distance_max = 0.0;	
    for (i = 0; i < nrows; i++)
    {
      X[i] = (B[i] / A[i][i]);
      for (j = 0; j < nrows; j++)
      {
        if (j == i) // Without the diagonals
          continue;
        X[i] = X[i] - ((A[i][j] / A[i][i]) * XPrevious[j]);
      }
      // SOR factor omega calculated on the next line
      X[i] = omega*X[i] +((1-omega)*XPrevious[i]);
      distance = abs( X[i] - XPrevious[i] );
	  distance_max = std::max(distance, distance_max);
      XPrevious[i] = X[i];
    }
    if ( (distance_max) < tolerance) {
      std::cout << "The convergence distance is less than the " << tolerance << "." << endl;
	  std::cout << "The number of iterations is " << Iter << "." << endl;
	  convergence_check = true;
      break;
    }
    if (Iter >= Max_Iter) {
	  std::cout << "The distance_max is the following: " << distance_max << endl;
	  std::cout << "The tolerance is the following: " << tolerance << endl;
	  std::cout << "The number of iterations is " << Iter << "." << endl;
      std::cout << "Maximum number of iterations reached before Jacobi_Iterative_Solve could finish." << std::endl;
      convergence_check = true;
      break;
    }    
    Iter++;
  }
  convergence_check = false;
  XPrevious.clear();
  XPrevious.shrink_to_fit();
  return 0;
}

int Gauss_Seidel (vector<double> A[], double B[], double X[], unsigned int nrows, unsigned int ncols, double tolerance, unsigned int Max_Iter, bool convergence_check)
{
  std::cout << endl << "Implementing the Gauss Seidel method: " << endl;  
  unsigned int i {0}, j {0}, Iter {1};
  double elem {0.0}, distance {0.0}, distance_max {0.0};
  vector<double> XPrevious;
  for (unsigned int j=0; j<ncols; j++) {
    XPrevious.push_back(elem);
  }
  std::cout << std::setprecision(5);
  while(convergence_check != true) {    
    distance_max = 0.0;	
    for (i = 0; i < nrows; i++)
    {
      X[i] = (B[i] / A[i][i]);
      for (j = 0; j < nrows; j++)
      {
        if (j == i) // Without the diagonals
          continue;
        X[i] = X[i] - ((A[i][j] / A[i][i]) * XPrevious[j]);
      }
      distance = abs( X[i] - XPrevious[i] );
	  distance_max = std::max(distance, distance_max);
      XPrevious[i] = X[i];
    }
    if ( (distance_max) < tolerance) {
      std::cout << "The convergence distance is less than the " << tolerance << "." << endl;
	  std::cout << "The number of iterations is " << Iter << "." << endl;
	  convergence_check = true;
      break;
    }
    if (Iter >= Max_Iter) {
	  std::cout << "The distance_max is the following: " << distance_max << endl;
	  std::cout << "The tolerance is the following: " << tolerance << endl;
	  std::cout << "The number of iterations is " << Iter << "." << endl;
      std::cout << "Maximum number of iterations reached before Jacobi_Iterative_Solve could finish." << std::endl;
      convergence_check = true;
      break;
    }    
    Iter++;
  }
  convergence_check = false;
  XPrevious.clear();
  XPrevious.shrink_to_fit();
  return 0;
}

int Jacobi(vector<double> A[], double B[], double X[], unsigned int nrows, unsigned int ncols, double tolerance, unsigned int Max_Iter, double convergence_check)
{
  std::cout << endl << "Implementing the Jacobi method: " << endl;
  unsigned int i {0},j {0}, Iter {1};
  double distance {0.0}, distance_max {0.0}, sum {0.0};
  vector<double> XPrevious;
  for (unsigned int j=0; j<ncols; j++) {
    XPrevious.push_back(0.0);
  }
  while(convergence_check != true) {	    
    distance_max = 0.0;
    for (i = 0; i<nrows; i++)
    {
      sum = 0.0;
      for (j=0; j<ncols; j++)
      {
        if (i != j)
        {
          sum = sum + A[i][j] * XPrevious[j];
        }
      }
      X[i] = -sum / (A[i][i]) + B[i] / (A[i][i]);
    }
    for (j=0;j<ncols;j++)
    {
	  distance = abs( X[j] - XPrevious[j] );
	  distance_max = std::max(distance, distance_max);
	  XPrevious[j] = X[j];
    }
    if ( (distance_max) < tolerance) {
      std::cout << "The convergence distance is less than the " << tolerance << "." << endl;
	  std::cout << "The number of iterations is " << Iter << "." << endl;
	  convergence_check = true;
      break;
    }
    if (Iter >= Max_Iter) {
	  std::cout << "The distance_max is the following: " << distance_max << endl;
	  std::cout << "The tolerance is the following: " << tolerance << endl;
	  std::cout << "The number of iterations is " << Iter << "." << endl;
      std::cout << "Maximum number of iterations reached before Jacobi could finish." << std::endl;
      convergence_check = true;
      break;
    }
    Iter++;
  }
  convergence_check = false;
  XPrevious.clear();
  XPrevious.shrink_to_fit();
  return 0;
}

int Output_Data(vector<double> A[], double B[], double X[], unsigned int nrows, unsigned int ncols)
{
  unsigned int i;
  std::cout << std::setprecision (5) << endl;
  std::cout << "******************** Solve Ax = B ********************" << endl;
  // Displaying the 2D vector
  std::cout << "The vector A is the following: " << endl; 
  PrintVector2D(&A[0], nrows);  
  std::cout << "The vector B is the following: " << endl; 
  PrintVector(&B[0], ncols); 
  // Create a new file named "X.dat"
  std::ofstream outputFile("X.dat");
  if (outputFile.is_open()) {
    // Write some text into the file
    outputFile << "%%MatrixMarket_Output_vector_X.dat matrix coordinate pattern general" << endl;
    outputFile << 1 << " " << nrows << " " << nrows;
    outputFile << endl;
    for (i = 0; i < ncols; i++)
      outputFile <<  1 << " " << i+1 << " " << X[i] << endl;
    // Close the file
    outputFile.close();
  } else {
    std::cout << "Error writing to file." << std::endl;
  }
  std::cout << "The vector X is the following: " << endl;
  PrintVector(&X[0], nrows);
  return 0;
}

/*
Iterate over vector of vectors and for each of the 
nested vector print its contents
*/
void PrintVector2D(const vector<double> A[], unsigned int nrows)
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
      auto it = A[i].begin();
        it != A[i].end(); it++)
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
    unsigned int i,j;
    FILE *myFile;
    myFile = fopen ("A.dat", "r");
    unsigned int dimArray[3];
    unsigned int nrows, number_of_entries_A;
    unsigned int ncols, Max_Iter {100};
    bool convergence_check = false;
    unsigned int i_index, j_index;
    double value, elem {0.0}, tolerance {0.00001};
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
    for (i = 0; i < number_of_entries_B; i++)
    {
      fscanf (myFile2, "%u,", &i_index);
      i_index--;
      fscanf (myFile2, "%u,", &j_index);
      j_index--;
      fscanf (myFile2, "%lf,", &value);
      B.push_back(value);
      X.push_back(elem);
    }
    fclose (myFile2);
    // Initializing the vector of vectors
//    vector<vector<double> > A;
    matrix<double> A;
    // Inserting elements into vector
    for (i = 0; i < nrows; i++) {
      // Vector to store column elements
      vector<double> v1;  
        for (j = 0; j < ncols; j++) {
		  elem = array_A[i*nrows+j];
          v1.push_back(elem);           
            if (j == (ncols)) {
				v1.push_back(elem);
			}
        }
        // Pushing back above 1D vector
        // to create the 2D vector
        A.push_back(v1);
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
     
    Jacobi (&A[0], &B[0], &X[0], nrows, ncols, tolerance, Max_Iter, convergence_check);
    Output_Data (&A[0], &B[0], &X[0], nrows, ncols);

    Gauss_Seidel (&A[0], &B[0], &X[0], nrows, ncols, tolerance, Max_Iter, convergence_check);
    Output_Data (&A[0], &B[0], &X[0], nrows, ncols);

    SOR (&A[0], &B[0], &X[0], nrows, ncols, tolerance, Max_Iter, convergence_check);
    Output_Data (&A[0], &B[0], &X[0], nrows, ncols);
     
    array_A.clear();
    array_A.shrink_to_fit();
    B.clear();
    B.shrink_to_fit();
    X.clear();
    X.shrink_to_fit();
    Augmented.clear();
    Augmented.shrink_to_fit();

//    std::cin.clear(); // reset any error flags
//    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); 
//    ignore any characters in the input buffer until we find an enter character
//    std::cin.get(); // get one more char from the user
    return 0;
  }  
};

int main ()
{
  Input Input_One;
  Input_One.Input_Begin();
  return 0;
}
