#include "rsvd.h"

Rsvd::RandomizedSvd<Eigen::MatrixXf, std::mt19937_64, Rsvd::SubspaceIterationConditioner::Lu>
randomized_svd(std::mt19937_64& random_engine,
               const Eigen::MatrixXf& m,
               const int rank)
{
  Rsvd::RandomizedSvd<Eigen::MatrixXf, std::mt19937_64, Rsvd::SubspaceIterationConditioner::Lu> rsvd(random_engine);

  rsvd.compute(m, rank);

  return rsvd;
}
