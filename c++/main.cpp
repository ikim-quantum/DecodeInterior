#include <armadillo>
#include <chrono>
#include <stdlib.h>
#include "apply_pauli.cpp"

using namespace std;
using namespace arma;

int main()
{
  int n = 32768;
  cx_vec psi1, psi2;

  psi1.randn(n);
  psi2.randn(n);


  psi2 = psi1;
  //cout << psi1 << endl;
   
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  apply_pauli_fast(0, true, true, psi1);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Fast algorithm" << endl;
  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
  //cout << psi1 << endl;

  
  //  cout << psi2 << endl;
  /* 
  std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
  apply_pauli_slow(0, true, true, psi2);
  std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
  std::cout << "Slow algorithm" << endl;
  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count() << "[µs]" << std::endl;
  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end2 - begin2).count() << "[ns]" << std::endl;
  */
  cout << sum(abs(psi1-psi2)) << endl;
}
