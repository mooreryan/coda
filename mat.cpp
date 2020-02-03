#include "mat.h"

/// Gives a 1D column major index given the row index and column index.
size_t
col_major_index(const size_t nrows, const size_t ridx, const size_t cidx)
{
  return cidx * nrows + ridx;
}

/// In-place CLR transformation
void
clr_in_place(Eigen::MatrixXf& otu_table)
{
  for (int j = 0; j < otu_table.cols(); ++j) {
    Eigen::ArrayXf log_col      = otu_table.col(j).array().log();
    float   mean_log_col = log_col.sum() / log_col.size();

    for (int i = 0; i < otu_table.rows(); ++i) {
      otu_table(i, j) = logf(otu_table(i, j)) - mean_log_col;
    }
  }
}

/// Euclidean distance between 2 vectors
float
distance(const Eigen::VectorXf& v1, const Eigen::VectorXf& v2)
{
  return sqrtf((v1 - v2).array().pow(2).sum());
}

/// @brief All v. all symmetric distance by columns
/// @param m matrix
/// @return a m.cols() by m.cols() matrix with distances
Eigen::MatrixXf
colwise_distance(const Eigen::MatrixXf& m)
{
  Eigen::MatrixXf d(m.cols(), m.cols());

  // All distances default to zero.
  d.setZero();

  for (int i = 0; i < m.cols() - 1; ++i) {
    for (int j = i + 1; j < m.cols(); ++j) {
      float dist = distance(m.col(i), m.col(j));

      d(i, j) = dist;
      d(j, i) = dist;
    }
  }

  return d;
}
