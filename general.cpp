#include <iostream>
#include <cmath>




double** Allocate_matrix(int Nrows, int Ncols)
{
  double** M;
  M = new double*[Nrows];
  for (int i=0; i<Nrows; i++)
    {
      M[i] = new double[Ncols];
    }
  return M;
}

void Deallocate_matrix(int Nrows, double** M)
{
  for (int i=0; i<Nrows; i++)
    {
      delete[] M[i];
    }
  delete[] M;
}

double* Allocate_vector(int n)
{
  double* vector;
  vector = new double[n];

  return vector;
}

void Deallocate_vector(double* vector)
{
  delete[] vector;
}



void Print_vector(int n,double* vec)
{
    for (int i=0;i<n;i++)
    {
        std::cout << vec[i] << std::endl;
    }
}
