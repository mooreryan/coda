#ifndef CODA_MAT_H
#define CODA_MAT_H

#include <Eigen/Dense>

/// Output format for vector
#define coeff_separator "\t"
#define row_separator "\n"

const Eigen::IOFormat TSVFormat(Eigen::StreamPrecision,
                                Eigen::DontAlignCols,
                                coeff_separator,
                                row_separator);


void
clr_in_place(Eigen::MatrixXf& otu_table);

float
distance(const Eigen::VectorXf& v1, const Eigen::VectorXf& v2);

/// Gives a 1D column major index given the row index and column index.
size_t
col_major_index(const size_t nrows, const size_t ridx, const size_t cidx);

Eigen::MatrixXf
colwise_distance(const Eigen::MatrixXf& m);

#endif //CODA_MAT_H
