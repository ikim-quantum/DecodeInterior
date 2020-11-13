#include "rand_state.cpp"
#include "entropy.cpp"
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

int main()
{
  int n_q = 10;
  int k = 3;
  for (int d=1; d<100; d++)
    {
      cx_dvec psi = scrambled_1d(n_q, d);
      cout << "depth=" << d << endl;
      cout << "von Neumann entropy=" << vNentropy(k, psi) << endl;
      cout << "diagonal entropy=" << diagentropy(k, psi) << endl;
    }

  return 0;
}
