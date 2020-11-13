// Testing PPR
#include "pauli_product.cpp"
#include "decode.cpp"
#include <armadillo>
#include <chrono>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace arma;

int main()
{
  int n = 256;
  cx_vec psi;
  psi.randn(n);
  psi = normalise(psi, 2);

  decode(psi, 6, 0.01);
}
