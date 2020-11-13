#include "complexity.hpp"
#include "rand_state.cpp"
#include "decode.cpp"


double complexity_rand_1d(int n_q, int d, int k, double eps, int samples)
{
  double cost_total = 0.0;
  for (int i=0; i<samples; i++)
    {
      //      cout << "sample=" << i << endl;
      cx_dvec psi=scrambled_1d(n_q, d);
      //decode(psi, k, eps);
      cost_total += decode_cost(0, psi, k, eps);
    }
  return cost_total/samples;
}

double complexity_haar(int n_q, int k, double eps, int samples)
{
  double cost_total = 0.0;
  for (int i=0; i<samples; i++)
    {
      cx_dvec psi= rand_haar(n_q);
      cost_total += decode_cost(0, psi, k, eps);
    }
  return cost_total/samples;
}
