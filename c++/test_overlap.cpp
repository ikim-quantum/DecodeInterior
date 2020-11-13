// Testing PPR
#include "pauli_product.cpp"
#include "decode.cpp"
#include <armadillo>
#include <chrono>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace arma;

// Checking if the optimized overlap function is working properly.
int main()
{
  int n = 256;
  int n_itt = 100;
  unsigned int x, z;


  cx_vec psi, psi_cpy;
  psi.randn(n);
  psi = normalise(psi, 2);

  double a1, a2, a3, theta;

  for (int x=0;x<n;x++)
    {
      for (int z=0;z<n;z++)
	{
	  a1 = coeff_a1(0, psi);
	  a2 = coeff_a2(0, x, psi);
	  a3 = coeff_a3(0, x, z, psi);	  
	  theta = optimized_theta(0, x, z, psi);
	  apply_ppr(x, z, theta, psi);
	  double my_overlap = optimal_overlap(0, x, z, psi);
	  double actual_overlap = coeff_a1(0, psi);
	  apply_ppr(x, z, -theta, psi);

	  if (abs(actual_overlap-my_overlap)>0.00001)
	    {
	      cout << psi << endl;
	      cout << "x, z=" << x << ", " << z << endl;
	      cout << "a1=" << a1 << endl;
	      cout << "a2=" << a2 << endl;
	      cout << "a3=" << a3 << endl;
	      cout << "theta_opt=" << theta << endl;
	      cout << "optimal_overlap=" << my_overlap << endl;
	      cout << "actual overlap after the update: " << actual_overlap << endl;
	    }
	}
    }
   
}
