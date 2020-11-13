#include "pauli_product.hpp"
#include <armadillo>
#include <cmath>
#include <stdlib.h>


using namespace std;
using namespace arma;

// Applies a Pauli operator i^{xz} X^xZ^z to the k'th qubit to the state psi.
// This algorithm takes at most 2^n complex number multiplications.
void apply_pauli_fast(int k, bool x, bool z, cx_dvec &psi)
{
  int d;
  d = psi.size();
  if (pow(2,k) >= d)
    {
      throw "Qubit index out of range.";
    }
  if (z)
    {
      if (x)
	{
	  for (int i=0; i<d; i++)
	    {
	      if ((i & ( 1 << k )) >> k == 0)
		{
		  cx_double temp1, temp2;
		  temp1 = psi(i) * cx_double(0.0, 1.0);
		  temp2 = psi((1<<k)|i) * cx_double(0.0, -1.0);
		  psi(i) = temp2;
		  psi((1<<k)|i) = temp1;
		}
	    }
	}
      else
	{
	  for (int i=0; i<d; i++)
	    {
	      if ((i & ( 1 << k )) >> k == 1)
		{
		  psi(i) *= -1.0;
		}
	    }
	}
    }
  else
    {
      if (x)
	{
	  for (int i=0; i<d; i++)
	    {
	      if ((i & ( 1 << k )) >> k == 0)
		{
		  cx_double temp1, temp2;
		  temp1 = psi(i);
		  temp2 = psi((1<<k)|i);
		  psi(i) = temp2;
		  psi((1<<k)|i) = temp1;
		}
	    }
	}
    }

}

// Applies a Pauli operator i^{xz} X^xZ^z to the k'th qubit to the state psi.
// This algorithm takes O(4^n)-time, quadratically worse than the apply_pauli_fast.
// This function is primarily used as a bookkeeping device, instead of performance.
void apply_pauli_slow(int k, bool x, bool z, cx_dvec &psi)
{
  int d, d_1, d_2;
  cx_mat PI(2, 2, fill::zeros), PX(2, 2, fill::zeros), PY(2, 2, fill::zeros), PZ(2, 2, fill::zeros);
  cx_mat p;

  PI(0, 0) = cx_double(1.0, 0.0); PI(1, 1) = cx_double(1.0, 0.0);
  PX(0, 1) = cx_double(1.0, 0.0); PI(1, 0) = cx_double(1.0, 0.0);
  PY(0, 1) = cx_double(0.0, -1.0); PY(1, 0) = cx_double(0.0, 1.0);
  PZ(0, 0) = cx_double(1.0, 0.0); PZ(1, 1) = cx_double(-1.0, 0.0);
  
  d = psi.size();
  d_1 = pow(2, k);
  d_2 = d / 2 / d_1;

  if (x)
    {
      if (z)
	{
	  p = PY;
	}
      else
	{
	  p = PX;
	}
       
    }
  else
    {
      if (z)
	{
	  p = PZ;
	}
      else
	{
	  p = PI;
	}
    }
  

  cx_mat I1(d_1, d_1,fill::eye), I2(d_2, d_2, fill::eye);

  cx_mat Pauli;
  Pauli = kron(kron(I2, p), I1);
  psi = Pauli * psi;
      
}

// Apply exp(iPθ) to psi, where P = i^{x·z}(⨂_k X_k^(x_k)) (⨂_k Z_k^(z_k))
void apply_ppr(unsigned int x, unsigned int z, double theta, cx_dvec &psi)
{
  int d;
  d = psi.size();
  cx_double c, s;
  c = cos(theta);
  s = sin(theta);
  //  cx_double csc;
  //  csc = 1.0 / cos(theta);
  bool not_changed[2048];
  fill_n(not_changed, 2048, true);

  for (int y=0; y<d; y++)
    {
      if (not_changed[y])
	{
	  
	    unsigned int xz = __builtin_popcount(x&z);
	    unsigned int yz = __builtin_popcount(y&z);
	    
	    cx_double phase1, phase2;

	    switch((1 + xz + 2*(yz+xz))%4)
	      {
	      case 0:
		phase1 = cx_double(1.0, 0.0);
		break;
	      case 1:
		phase1 = cx_double(0.0, 1.0);
		break;
	      case 2:
		phase1 = cx_double(-1.0, 0.0);
		break;
	      case 3:
		phase1 = cx_double(0.0, -1.0);
		break;
	      }
	    
	    switch((1 + xz + 2 * yz)%4)
	      {
	      case 0:
		phase2 = cx_double(1.0, 0.0);
		break;
	      case 1:
		phase2 = cx_double(0.0, 1.0);
		break;
	      case 2:
		phase2 = cx_double(-1.0, 0.0);
		break;
	      case 3:
		phase2 = cx_double(0.0, -1.0);
		break;
	      }
	    cx_double temp1, temp2;
	    temp1 = psi(y) * c + psi(y^x) * s * phase1;
	    temp2 = psi(y^x) * c + psi(y) * s * phase2;
	    psi(y) = temp1;
	    psi(y^x) = temp2;

	    //psi(y) *= c; psi(y) += s * phase1 * psi(y^x);
	    //psi(y^x) += s * phase2 * psi(y);
	    //psi(y^x) *= csc;

	    not_changed[y] = false; not_changed[y^x] = false;
	}
      
    }
}

void apply_ppr_slow(unsigned int x, unsigned int z, double theta, cx_dvec &psi)
{
    int d;
    d = psi.size();
    unsigned int bits, n_temp = (d < 0) ? -d : d;
    for(bits = 0; n_temp != 0; ++bits) n_temp >>= 1;
    unsigned int n_q = bits - 1;
    
    cx_double c, s;
    c = cos(theta);
    s = sin(theta);

    cx_dmat PX(2,2, fill::zeros), PY(2,2, fill::zeros), PZ(2,2, fill::zeros), PI(2,2, fill::zeros), P, PString(1,1, fill::ones);
    PX(0, 1) = 1.0; PX(1, 0) = 1.0;
    PZ(0, 0) = 1.0; PZ(1, 1) = -1.0;
    PY(1, 0) = cx_double(0.0, 1.0); PY(0, 1) = cx_double(0.0, -1.0);
    PI(0, 0) = 1.0; PI(1, 1) = 1.0;
		   
    for (int i=0; i<n_q; i++)
      {
	
	if (x & (1<< i))
	  {
	    if (z & (1<< i))
		P = PY;  
	    else
		P = PX;
	  }
	else
	  {
	    if (z & (1<< i))
		P = PZ;
	    else
		P = PI;
	  }
	PString = kron(P, PString);
      }

    psi = c * psi + PString * psi * s * cx_double(0.0, 1.0);

}

// Measure Pauli Product operator (fast)
double measure_pp(unsigned int x, unsigned int z, cx_dvec &psi)
{
  int d;
  d = psi.size();
  unsigned int bits, n_temp = (d < 0) ? -d : d;
  for(bits = 0; n_temp != 0; ++bits) n_temp >>= 1;
  unsigned int n_q = bits - 1;

  cx_double mysum = 0.0;
  for (int y=0; y<d; y++)
    {
      unsigned int xz = __builtin_popcount(x&z);
      unsigned int yz = __builtin_popcount(y&z);
      cx_double phase;
      switch((xz + 2*yz)%4)
	{
	case 0:
	  phase = cx_double(1.0, 0.0);
	  break;
	case 1:
	  phase = cx_double(0.0, 1.0);
	  break;
	case 2:
	  phase = cx_double(-1.0, 0.0);
	  break;
	case 3:
	  phase = cx_double(0.0, -1.0);
	  break;
	}
      mysum += psi(y) * conj(psi(y^x)) * phase;
    }
  return real(mysum);
}


// Meausre Pauli Product operator (slow)
double measure_pp_slow(unsigned int x, unsigned int z, cx_dvec &psi)
{
  int d;
  d = psi.size();
  unsigned int bits, n_temp = (d < 0) ? -d : d;
  for(bits = 0; n_temp != 0; ++bits) n_temp >>= 1;
  unsigned int n_q = bits - 1;

  cx_dmat PX(2,2, fill::zeros), PY(2,2, fill::zeros), PZ(2,2, fill::zeros), PI(2,2, fill::zeros), P, PString(1,1, fill::ones);
  PX(0, 1) = 1.0; PX(1, 0) = 1.0;
  PZ(0, 0) = 1.0; PZ(1, 1) = -1.0;
  PY(1, 0) = cx_double(0.0, 1.0); PY(0, 1) = cx_double(0.0, -1.0);
  PI(0, 0) = 1.0; PI(1, 1) = 1.0;
  
  for (int i=0; i<n_q; i++)
    {
      
      if (x & (1<< i))
	{
	  if (z & (1<< i))
	    P = PY;  
	  else
	    P = PX;
	}
      else
	{
	  if (z & (1<< i))
	    P = PZ;
	  else
	    P = PI;
	}
      PString = kron(P, PString);
    }
  return real(cdot(psi, PString * psi));
}
