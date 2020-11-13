#ifndef RAND_STATE_HPP_
#define RAND_STATE_HPP_

#include <armadillo>

// Random number generator seed needed?
void apply_rand_gate(unsigned int q1, unsigned int q2, arma::cx_dvec &psi);
arma::cx_dvec rand_haar(unsigned int n_q);
arma::cx_dvec scrambled_1d(unsigned int n_q, unsigned int depth);
arma::cx_dvec scrambled_prosen(unsigned int n_q);

#endif // RAND_STATE_HPP_
