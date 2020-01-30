#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <Eigen/Dense>
#include <rsvd/Prelude.hpp>
#include <assert.h>

#define VERSION "0.1.0"

using namespace std;
using namespace Eigen;

using Rsvd::RandomizedSvd;
using Rsvd::SubspaceIterationConditioner;

void
log_msg(const char* msg)
{
  fprintf(stderr, "LOG -- %s\n", msg);
}

RandomizedSvd<MatrixXf, mt19937_64, SubspaceIterationConditioner::Lu>
randomized_svd(mt19937_64& random_engine,
               const MatrixXf& m,
               const int rank)
{
  RandomizedSvd<MatrixXf, mt19937_64, SubspaceIterationConditioner::Lu> rsvd(random_engine);

  rsvd.compute(m, rank);

  return rsvd;
}

/// Output format for vector
#define coeff_separator "\t"
#define row_separator "\n"
const

static IOFormat TSVFormat(StreamPrecision,
                          DontAlignCols,
                          coeff_separator,
                          row_separator);

/// What DivNet uses
const float zero_replacement = 0.05;

/// In-place CLR transformation
void
clr_in_place(MatrixXf& otu_table)
{
  for (int j = 0; j < otu_table.cols(); ++j) {
    ArrayXf log_col      = otu_table.col(j).array().log();
    float   mean_log_col = log_col.sum() / log_col.size();

    for (int i = 0; i < otu_table.rows(); ++i) {
      otu_table(i, j) = logf(otu_table(i, j)) - mean_log_col;
    }
  }
}

/// Euclidean distance between 2 vectors
float
distance(const VectorXf& v1, const VectorXf& v2)
{
  return sqrtf((v1 - v2).array().pow(2).sum());
}

/// @brief All v. all symmetric distance by columns
/// @param m matrix
/// @return a m.cols() by m.cols() matrix with distances
MatrixXf
colwise_distance(const MatrixXf& m)
{
  MatrixXf d(m.cols(), m.cols());

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

vector<string>
split(const string& s, char delimiter)
{
  vector<string> tokens;
  string         token;
  istringstream  token_stream(s);

  while (getline(token_stream, token, delimiter)) {
    tokens.push_back(token);
  }

  return tokens;
}

/// Will throw if there are more tokens than will fit in `tokens`.
void
split(array<string, 4>& tokens, const string& s, char delimiter)
{
  string        token;
  istringstream token_stream(s);

  int i    = 0;
  int size = tokens.size();
  while (getline(token_stream, token, delimiter)) {
    if (i >= size) {
      throw runtime_error("Too many tokens; won't fit!");
    }

    tokens[i++] = token;
  }
}

typedef struct Element
{
  string sample;
  string otu;
  size_t sample_index;
  size_t otu_index;
  float  otu_len;
  float  count;
}           Element;

void
element_print(const Element& element)
{
  fprintf(stdout,
          "sample: %s, sample_index: %lu, otu: %s, otu_index: %lu, otu_len: %f, count: %f\n",
          element.sample.c_str(),
          element.sample_index,
          element.otu.c_str(),
          element.otu_index,
          element.otu_len,
          element.count);
}


// TODO does this copy the strings or just refs?
Element
parse_data_line(const array<string, 4>& tokens)
{
  Element element;

  element.sample       = tokens[0];
  element.otu          = tokens[1];
  element.sample_index = 0;
  element.otu_index    = 0;
  element.otu_len      = stof(tokens[2]);
  element.count        = 0.0f;

  int count = stoi(tokens[3]);

  // Replace zero counts with a small non-zero number.
  if (count == 0) {
    element.count = 0.05f;
  }
  else {
    element.count = (float)count;
  }

  return element;
}

/// Counts/hits per 100 AAs/bases
float
element_normalize_count(const Element& element)
{
  return element.count / element.otu_len * 100;
}


/// Gives a 1D column major index given the row index and column index.
size_t
col_major_index(const size_t nrows, const size_t ridx, const size_t cidx)
{
  return cidx * nrows + ridx;
}

int main(int argc, char* argv[])
{
  if (argc != 6) {
    fprintf(stderr,
            "VERSION: %s"
            "\n\n"
            "Usage: %s <seed> <counts> <clr_out.tsv> <aitchison_dist_out.tsv> <sample_projection_out.tsv>\n\n"
            "Rows are OTUs, columns are samples.\n\n"
            "IMPORTANT:  assumes more OTUs than samples, but doesn't check for it!\n\n",
            VERSION,
            argv[0]);

    return 1;
  }

  char* seed_str                = argv[1];
  char* in_fname                = argv[2];
  char* clr_out_fname           = argv[3];
  char* ait_out_fname           = argv[4];
  char* sample_projection_fname = argv[5];

  int seed;
  try {
    seed = stoi(seed_str);
  }
  catch (const invalid_argument& ia) {
    cerr << "ERROR -- could not convert the seed!" << endl;
    return 1;
  }
  catch (const out_of_range& oof) {
    cerr << "ERROR -- the seed you entered is out of range!" << endl;
    return 1;
  }

  ofstream clr_out;
  ofstream ait_out;
  ofstream sample_projection_out;

  clr_out.open(clr_out_fname);
  ait_out.open(ait_out_fname);
  sample_projection_out.open(sample_projection_fname);

  if (!clr_out.is_open()) {
    fprintf(stderr,
            "ERROR -- Couldn't open '%s' for writing.",
            clr_out_fname);
  }

  if (!ait_out.is_open()) {
    fprintf(stderr,
            "ERROR -- Couldn't open '%s' for writing.",
            ait_out_fname);
  }

  if (!sample_projection_out.is_open()) {
    fprintf(stderr,
            "ERROR -- Couldn't open '%s' for writing.",
            sample_projection_fname);
  }

  array<string, 4> tokens;
  vector<Element>  elements;

  // Map samples to columns.  Will be columns.
  unordered_map<string, size_t> sample_index;

  // Map otu to rows.  Will be rows.
  unordered_map<string, size_t> otu_index;

  string   line;
  ifstream in_file(in_fname);

  if (in_file.is_open()) {
    while (getline(in_file, line)) {
      split(tokens, line, '\t');
      Element element = parse_data_line(tokens);

      sample_index.try_emplace(element.sample, sample_index.size());
      otu_index.try_emplace(element.otu, otu_index.size());

      element.sample_index = sample_index.at(element.sample);
      element.otu_index    = otu_index.at(element.otu);

      elements.push_back(element);
    }
    in_file.close();
  }
  else {
    cerr << "ERROR -- Couldn't open the file!";
    return 1;
  }

  // Need a map from index to sample name.
  unordered_map<size_t, string> index_to_sample;
  for (const auto& pair : sample_index) {
    index_to_sample.try_emplace(pair.second, pair.first);
  }

  // Need a map from index to otu_name name.
  unordered_map<size_t, string> index_to_otu;

  for (const auto& pair : otu_index) {
    index_to_otu.try_emplace(pair.second, pair.first);
  }

  // nrows = nOTUs, ncols = nsamples
  MatrixXf otu_table(otu_index.size(), sample_index.size());
  // any missing values will be set to zero
  otu_table.setZero();

  for (const auto& element : elements) {
    float norm_count = element_normalize_count(element);

    otu_table(element.otu_index, element.sample_index) = norm_count;
  }

  log_msg("Calculating CLR transformation");

  clr_in_place(otu_table);

  log_msg("Printing CLR matrix");

  // Print CLR matrix
  clr_out << "otu";
  for (size_t i = 0; i < index_to_sample.size(); ++i) {
    clr_out << "\t" << index_to_sample.at(i);
  }
  clr_out << endl;

  for (Index i = 0; i < otu_table.rows(); ++i) {
    string otu = index_to_otu.at(i);
    clr_out << otu
            << "\t"
            << otu_table.row(i).format(TSVFormat)
            << endl;
  }

  log_msg("Calculating Aitchison distance");

  MatrixXf aitchison_distance = colwise_distance(otu_table);

  log_msg("Printing Aitchison distance matrix");

  // Print Aitchison distance
  ait_out << "sample";
  for (size_t i = 0; i < index_to_sample.size(); ++i) {
    ait_out << "\t" << index_to_sample.at(i);
  }
  ait_out << endl;

  for (int i = 0; i < aitchison_distance.rows(); ++i) {
    string sample = index_to_sample.at(i);
    ait_out << sample
            << "\t"
            << aitchison_distance.row(i).format(TSVFormat)
            << endl;
  }

  std::mt19937_64 random_engine{};
  random_engine.seed(seed);

  // Do the SVD
  RandomizedSvd<MatrixXf, std::mt19937_64, SubspaceIterationConditioner::Lu> rsvd(random_engine);

  rsvd.compute(otu_table, otu_table.cols());

  log_msg("Calculating SVD");

  const MatrixXf sample_projection = rsvd.matrixV() * rsvd.singularValues().asDiagonal();

  assert(sample_projection.cols() == sample_projection.rows());

  log_msg("Writing sample projection");

  sample_projection_out << "sample";
  for (size_t i = 0; i < index_to_sample.size(); ++i) {
    sample_projection_out << "\tPC" << (i + 1);
  }
  sample_projection_out << endl;

  for (Index i = 0; i < sample_projection.rows(); ++i) {
    string sample = index_to_sample.at(i);
    sample_projection_out << sample
                          << "\t"
                          << sample_projection.row(i).format(TSVFormat)
                          << endl;
  }

  return 0;
}
