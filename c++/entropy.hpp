#ifndef ENTROPY_HPP_
#define ENTROPY_HPP_

#include <armadillo>

double vNentropy(unsigned int k, arma::cx_dvec &psi);
double diagentropy(unsigned int k, arma::cx_dvec &psi);

#endif // ENTROPY_HPP_
