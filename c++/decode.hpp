#ifndef QSIM_DECODE_HPP_
#define QSIM_DECODE_HPP_
#include <armadillo>

double coeff_a1(unsigned int q_idx, arma::cx_dvec &psi);
double coeff_a2(unsigned int q_idx, unsigned int x, arma::cx_dvec &psi);
double coeff_a3(unsigned int q_idx, unsigned int x, unsigned int z, arma::cx_dvec &psi);
double optimized_theta(unsigned int q_idx, unsigned int x, unsigned int z, arma::cx_dvec &psi);
double optimal_overlap(unsigned int q_idx, unsigned int x, unsigned int z, arma::cx_dvec &psi);
void optimal_update(unsigned int setbits, arma::cx_dvec &psi);
void decode(arma::cx_dvec &psi, unsigned int setbits, double eps);
double decode_cost(arma::cx_dvec &psi, unsigned int setbits, double eps);

#endif // QSIM_DECODE_HPP_
