#include "rand_state.hpp"
#include "pauli_product.cpp"
#include <armadillo>
#include <cmath>
#include <stdlib.h>

using namespace std;
using namespace arma;

void apply_rand_gate(unsigned int q1, unsigned int q2, cx_dvec &psi)
{
  int x1, x2, z1, z2;
  x1 = rand()%2;
  x2 = rand()%2;
  z1 = rand()%2;
  z2 = rand()%2;

  int x, z;
  x = (x1<<q1) + (x2<<q2);
  z = (z1<<q1) + (z2<<q2);

  double theta = (double)rand() / (double)RAND_MAX * 3.14159265358979;

  apply_ppr(x, z, theta, psi);
}

cx_dvec rand_haar(unsigned int n_q)
{
  cx_vec psi;
  psi.randn(1<<n_q);
  psi = normalise(psi, 2);
  return psi;
}

cx_dvec scrambled_1d(unsigned int n_q, unsigned int depth)
{
  cx_vec psi;
  psi.zeros(1<<n_q);
  psi(0) = 1;
  for (int d=0; d<depth; d++)
    {
      if (d%2==0)
	{
	  for (int i=0; i<n_q; i++)
	    {
	      if (i%2==0)
		{
		  if (i<n_q-1)
		      apply_rand_gate(i, i+1, psi);
		  //		  else
		  //		      apply_rand_gate(i,0, psi);
		}
	    }
	}
      else
	{
	  for (int i=0; i<n_q; i++)
	    {
	      if (i%2==1)
		{
		  if (i<n_q-1)
		      apply_rand_gate(i, i+1, psi);
		  //		  else
		  //		      apply_rand_gate(i,0, psi);
		}
	    }
	}
    }
  return psi;
}
