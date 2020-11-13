#include "entropy.hpp"
#include <armadillo>
#include <stdlib.h>
using namespace arma;
using namespace std;

// von Neumann entropy of the first k qubits
double vNentropy(unsigned int k, cx_dvec &psi)
{
  int d = psi.size();
  cx_dmat sqrtrho = reshape(psi, 1<<k, d/(1<<k));
  vec s = svd(sqrtrho);
  psi.reshape(size(psi));
  double sum = 0.0;
  for (int i=0; i<(1<<k); i++)
    {
      double ent_sq = norm(s(i));
      sum -= log2(ent_sq) * ent_sq;
    }
  return sum;
}

// Diagonal entropy of the first k qubits
double diagentropy(unsigned int k, cx_dvec &psi)
{
  int d = psi.size();
  double sum = 0.0;
  for (int i=0; i<(1<<k); i++)
    {
      double ent_sq = 0.0;
      for (int j=0; j<(d/(1<<k)); j++)
	{
	  ent_sq += norm(psi(i + (j<<k)));
	}
      sum -= log2(ent_sq) * ent_sq;
    }
  return sum;
}
