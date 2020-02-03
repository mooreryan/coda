#ifndef CODA_RSVD_H
#define CODA_RSVD_H

#include <rsvd/Prelude.hpp>
#include <Eigen/Dense>
#include <random>

Rsvd::RandomizedSvd<Eigen::MatrixXf, std::mt19937_64, Rsvd::SubspaceIterationConditioner::Lu>
randomized_svd(std::mt19937_64& random_engine,
               const Eigen::MatrixXf& m,
               const int rank);

#endif //CODA_RSVD_H
