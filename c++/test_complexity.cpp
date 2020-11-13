#include "complexity.cpp"

int main()
{
  int n_q = 4;
  int d_max= 100;
  int k = 4;
  int samples = 100;
  double eps = 0.1;

  for (int d=1; d<d_max; d++)
    {
      //      cout << "Random circuit" << endl;
      cout <<complexity_rand_1d(n_q, d, k, eps, samples) << ", " <<endl;
      //      cout << "Haar" << endl;
      //      cout <<complexity_haar(n_q, k, eps, samples) << ", " <<endl;
    }
  return 0;
}
