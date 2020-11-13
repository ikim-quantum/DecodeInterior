#include "pauli_product.cpp"
#include <armadillo>
#include <chrono>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace arma;

int main()
{
  int n = 1024;
  int n_itt = 100;

  for (int i=0; i<n_itt; i++)
    {
      cx_vec psi;
      psi.randn(n);
      psi = normalise(psi, 2);

      unsigned x = rand()%n;
      unsigned z = rand()%n;
      
      double m1 = measure_pp(x, z, psi);
      double m2 = measure_pp_slow(x, z, psi);

      cout << m1 - m2 << endl;
    }

  cx_vec psi;
  psi.randn(n);
  psi = normalise(psi, 2);
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  for (int i=0; i<n_itt; i++)
    {
      unsigned x = rand()%n;
      unsigned z = rand()%n;
      
      double m = measure_pp(x, z, psi);
    }
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  std::cout << "Measurement Time (fast) = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/n_itt << "[µs]" << std::endl;
  std::cout << "Measurement Time (fast) = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count()/n_itt << "[ns]" << std::endl;


  psi.randn(n);
  psi = normalise(psi, 2);
  std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
  for (int i=0; i<n_itt; i++)
    { 
      unsigned x = rand()%n;
      unsigned z = rand()%n;
      
      double m = measure_pp_slow(x, z, psi);
    }
  std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();

  std::cout << "Measurement Time (slow) = " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count()/n_itt << "[µs]" << std::endl;
  std::cout << "Measurement Time (slow) = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end2 - begin2).count()/n_itt << "[ns]" << std::endl;

}
