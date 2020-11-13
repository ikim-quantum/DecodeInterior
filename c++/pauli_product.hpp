#ifndef QSIM_APPLY_PAULI_HPP_
#define QSIM_APPLY_PAULI_HPP_

#include <armadillo>

void apply_pauli(int k, bool x, bool z, arma::cx_dvec &psi);
void apply_pauli_slow(int k, bool x, bool z, arma::cx_dvec &psi);
void apply_ppr(unsigned int x, unsigned int z, double theta, arma::cx_dvec &psi);
void apply_ppr_slow(unsigned int x, unsigned int z, double theta, arma::cx_dvec &psi);
double measure_pp(unsigned int x, unsigned int z, arma::cx_dvec &psi);
double measure_pp_slow(unsigned int x, unsigned int z, arma::cx_dvec & psi);

#endif // QSIM_APPLY_PAULI_HPP_
