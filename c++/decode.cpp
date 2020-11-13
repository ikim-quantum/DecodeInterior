#include "decode.hpp"
#include "pauli_product.hpp"
#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

double coeff_a1(unsigned int q_idx, cx_dvec &psi)
{
  int d = psi.size();
  double mysum = 0.0;
  for (int y=0; y<d; y++)
    {
      if (!(y & (1<<q_idx)))
	{
	   mysum += std::norm(psi(y));
	}
    }
  return mysum;
}

double coeff_a2(unsigned int q_idx, unsigned int x, cx_dvec &psi)
{
  int d = psi.size();
  double mysum = 0.0;
  for (int y=0; y<d ; y++)
    {
      int temp = 1<< q_idx;
      int temp2 = y^x;
      if (!((y^x) & (1<<q_idx)))
	{
	  mysum += std::norm(psi(y));
	}
      else
	{
	}
    }
  return mysum;
  
}

double coeff_a3(unsigned int q_idx, unsigned int x, unsigned int z, cx_dvec &psi)
{
  int d = psi.size();
  double mysum = 0.0;
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
      if (!((y^x) & (1<<q_idx)))
	{
	  //	  cout << y << endl;
	  //	  cout << (y^x) << endl;
	  mysum += (-2) * std::imag(std::conj(psi(y^x)) * psi(y) * phase);
	}
    }
  return mysum;
}

// Given a Pauli string, find the rotation angle with respect to the corresponding PPR that gives the optimal overlap with the |0>_{q_idx} state.
double optimized_theta(unsigned int q_idx, unsigned int x, unsigned int z, cx_dvec &psi)
{
  double a1 = coeff_a1(q_idx, psi);
  double a2 = coeff_a2(q_idx, x, psi);
  double a3 = coeff_a3(q_idx, x, z, psi);

  double theta1, theta2;
  if (a1==a2)
    theta1= 3.14159265358979/2.0;
  else
    theta1 = atan(a3/(a1 - a2))/2.0;
  
  theta2 = theta1 + 3.14159265358979/2.0;

  double f1 = ((a1-a2) * cos(2*theta1) + a3 * sin(2*theta1))/2.0;
  double f2 = ((a1-a2) * cos(2*theta2) + a3 * sin(2*theta2))/2.0;

  f1 += (a1+a2)/2.0;
  f2 += (a1+a2)/2.0;

  if (f1>f2)
    return theta1;
  else
    return theta2;
}

// Given a Pauli string, find the optimal overlap with the |0>_{q_idx} using
// a PPR.
double optimal_overlap(unsigned int q_idx, unsigned int x, unsigned int z, cx_dvec &psi)
{
  double a1 = coeff_a1(q_idx, psi);
  double a2 = coeff_a2(q_idx, x, psi);
  double a3 = coeff_a3(q_idx, x, z, psi);

  double theta1, theta2;
  if (a1==a2)
    theta1= 3.14159265358979/2.0;
  else
    theta1 = atan(a3/(a1 - a2))/2.0;
  
  theta2 = theta1 + 3.14159265358979/2.0;
  double f1 = ((a1-a2) * cos(2*theta1) + a3 * sin(2*theta1))/2.0;
  double f2 = ((a1-a2) * cos(2*theta2) + a3 * sin(2*theta2))/2.0;

  f1 += (a1+a2)/2.0;
  f2 += (a1+a2)/2.0;

  if (f1>f2)
    return f1;
  else
    return f2;
}

// Find the optimal Pauli string that maximizes the overlap with |0>_q
// and apply the corresponding PPR to psi. Here the PPR is assumed to be restricted
// to a set of first 'bits' qubits
void optimal_update(unsigned int q_idx, unsigned int setbits, cx_dvec &psi)
{
  int i_best, j_best;
  double f = 0.0;
  double theta = 0.0;
  for (int i=0; i<(1<<setbits); i++)
    {
      for (int j=0; j<(1<<setbits); j++)
	{
	  /*	  cout << "State=" << endl;
	  cout << psi << endl;
	  cout << "i, j=" << i << ", " << j << endl;
	  double a1 = coeff_a1(0, psi);
	  double a2 = coeff_a2(0, i, psi);
	  double a3 = coeff_a3(0, i, j, psi);
	  double f_opt = optimal_overlap(0, i, j, psi);
	  cout << "a1=" << a1 << endl;
	  cout << "a2=" << a2 << endl;
	  cout << "a3=" << a3 << endl;
	  cout << "f_opt=" << f_opt << endl;
	  if ((a1==a2) && (abs(a3)>0.0001))
	  exit (EXIT_FAILURE);*/

	  // Apply optimal update only if the PPR acts nontrivially on q_idx.
	  if ((i & (1<<q_idx)) || (j & (1<<q_idx)))
	    {
	      if (optimal_overlap(q_idx, i, j, psi) > f)
		{
		  f = optimal_overlap(q_idx, i, j, psi);
		  //	      cout << "New record: " << f << endl;
		  i_best = i;
		  j_best = j;
		}
	    }
	  
	}
    }
  theta = optimized_theta(0, i_best, j_best, psi);
  apply_ppr(i_best, j_best, theta, psi);
}

// Skip the given choice of i and j with probability 1-ratio.
void randomized_update(unsigned int q_idx, unsigned int setbits, double ratio, cx_dvec &psi)
{
  int i_best, j_best;
  double f = 0.0;
  double theta = 0.0;
  for (int i=0; i<(1<<setbits); i++)
    {
      for (int j=0; j<(1<<setbits); j++)
	{
	  if ((i & (1<<q_idx)) || (j & (1<<q_idx)))
	    {
	      if ((double)rand() / (double)RAND_MAX < ratio)
		{
		  if (optimal_overlap(q_idx, i, j, psi) > f)
		    {
		      f = optimal_overlap(q_idx, i, j, psi);
		      //	      cout << "New record: " << f << endl;
		      i_best = i;
		      j_best = j;
		    }
		}
	    }
	  
	}
    }
  theta = optimized_theta(0, i_best, j_best, psi);
  apply_ppr(i_best, j_best, theta, psi);
}

// Decode and print the result.
void decode(unsigned q_idx, cx_dvec &psi, unsigned int setbits, double eps)
{
  double f = 0.0;
  int itt = 0;
  while (f<1-eps)
    {
      itt++;
      optimal_update(q_idx, setbits, psi);
      f = coeff_a1(0, psi);
      cout << itt << ":" << f << endl;
    }
}

void randomized_decode(unsigned q_idx, cx_dvec &psi, unsigned int setbits, double ratio, double eps)
{
  double f = 0.0;
  int itt = 0;
  while (f<1-eps)
    {
      itt++;
      randomized_update(q_idx, setbits, ratio, psi);
      f = coeff_a1(0, psi);
      cout << itt << ":" << f << endl;
    }
}

double decode_cost(unsigned q_idx, cx_dvec &psi, unsigned int setbits, double eps)
{
  double f = 0.0;
  double cost = 0.0;
  int itt = 0;
  while (f<1-eps)
    {
      itt++;
      optimal_update(q_idx, setbits, psi);
      f = coeff_a1(0, psi);
      //      cout << itt << ":" << f << endl;
      cost += 1.0;
    }
  return cost;
}

double randomized_decode_cost(unsigned q_idx, cx_dvec &psi, unsigned int setbits, double ratio, double eps)
{
  double f = 0.0;
  double cost = 0.0;
  int itt = 0;
  while (f<1-eps)
    {
      itt++;
      randomized_update(q_idx, setbits, ratio, psi);
      f = coeff_a1(0, psi);
      //      cout << itt << ":" << f << endl;
      cost += 1.0;
    }
  return cost;
}
