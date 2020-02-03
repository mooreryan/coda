#include <stdlib.h> // exit, EXIT_FAILURE

#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>

#include "utils.h"
#include "mat.h"
#include "rsvd.h"
#include "element.h"

#define VERSION "0.2.0"

using namespace std;


int main(int argc, char* argv[])
{
  if (argc != 3) {
    fprintf(stderr,
            "VERSION: %s"
            "\n\n"
            "Usage: %s <seed> <counts>\n\n"
            "Rows are OTUs, columns are samples.\n\n"
            "IMPORTANT:  assumes more OTUs than samples, but doesn't check for it!\n\n",
            VERSION,
            argv[0]);

    return 1;
  }

  char* seed_str = argv[1];
  char* in_fname = argv[2];

  unsigned long zero_counts = 0;

  array<string, 4> tokens;
  vector <Element> elements;

  // Map samples to columns.  Will be columns.
  unordered_map <string, size_t> sample_index;

  // Map otu to rows.  Will be rows.
  unordered_map <string, size_t> otu_index;

  // Need a map from index to sample name.
  unordered_map <size_t, string> index_to_sample;

  // Need a map from index to otu_name name.
  unordered_map <size_t, string> index_to_otu;

  string   line;
  ifstream in_file(in_fname);

  // Parse seed
  int seed;
  try {
    seed = stoi(seed_str);
  }
  catch (const invalid_argument& ia) {
    cerr << "ERROR -- could not convert the seed!" << endl;
    return EXIT_FAILURE;
  }
  catch (const out_of_range& oof) {
    cerr << "ERROR -- the seed you entered is out of range!" << endl;
    return EXIT_FAILURE;
  }

  fprintf(stderr, "INFO -- seed: %d\n", seed);

  log_msg("Reading data");

  // Read input file
  if (in_file.is_open()) {
    while (getline(in_file, line)) {
      split(tokens, line, '\t');
      Element element = parse_data_line(tokens, true);

      // Track "zero" counts
      if (element.count < 1) {
        zero_counts++;
      }

      // map sample to index
      sample_index.try_emplace(element.sample, sample_index.size());

      // map otu to index
      otu_index.try_emplace(element.otu, otu_index.size());

      // track the indices of this element
      element.sample_index = sample_index.at(element.sample);
      element.otu_index    = otu_index.at(element.otu);

      elements.push_back(element);
    }
    in_file.close();
  }
  else {
    cerr << "ERROR -- Couldn't open the file!";
    return EXIT_FAILURE;
  }

  log_msg("Creating sample lookup table");
  // Need to look up samples and OTUs by index

  for (const auto& pair : sample_index) {
    index_to_sample.try_emplace(pair.second, pair.first);
  }

  log_msg("Creating OTU lookup table");
  for (const auto& pair : otu_index) {
    index_to_otu.try_emplace(pair.second, pair.first);
  }

  log_msg("Normalizing counts");

  // nrows = nOTUs, ncols = nsamples
  Eigen::MatrixXf otu_table = Eigen::MatrixXf::Constant(otu_index.size(), sample_index.size(), zero_replacement);

  fprintf(stderr,
          "INFO -- Table %% zeros: %.3f\n",
          zero_counts / (double)otu_table.size() * 100);

  // Convert input to normalized counts and into an Eigen matrix
  for (const auto& element : elements) {
    float norm_count = element_normalize_count(element);

    otu_table(element.otu_index, element.sample_index) = norm_count;
  }

  // Empty and free elements.
  elements.clear();
  vector<Element>().swap(elements);

  log_msg("Calculating CLR transformation");

  clr_in_place(otu_table);

  log_msg("Printing CLR matrix");

  // Print CLR matrix
  ofstream clr_out;
  clr_out.open("coda__clr.tsv");

  if (!clr_out.is_open()) {
    fprintf(stderr,
            "ERROR -- Couldn't open clr_outfile for writing.");
  }

  clr_out << "otu";
  for (size_t i = 0; i < index_to_sample.size(); ++i) {
    clr_out << "\t" << index_to_sample.at(i);
  }
  clr_out << endl;

  for (Eigen::Index i = 0; i < otu_table.rows(); ++i) {
    string otu = index_to_otu.at(i);
    clr_out << otu
            << "\t"
            << otu_table.row(i).format(TSVFormat)
            << endl;
  }

  clr_out.close();

  log_msg("Calculating Aitchison distance");

  Eigen::MatrixXf aitchison_distance = colwise_distance(otu_table);

  ofstream ait_out;
  ait_out.open("coda__ait.tsv");

  if (!ait_out.is_open()) {
    fprintf(stderr,
            "ERROR -- Couldn't open clr_outfile for writing.");
  }

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

  ait_out.close();

  // Do the SVD
  std::mt19937_64 random_engine{};
  random_engine.seed(seed);

  log_msg("Calculating SVD");

  Rsvd::RandomizedSvd<Eigen::MatrixXf, std::mt19937_64, Rsvd::SubspaceIterationConditioner::Lu> rsvd = randomized_svd(random_engine, otu_table, otu_table.cols());

  rsvd.compute(otu_table, otu_table.cols());

  const Eigen::MatrixXf sample_projection = rsvd.matrixV() * rsvd.singularValues().asDiagonal();

  assert(sample_projection.cols() == sample_projection.rows());

  // Write the SVD
  log_msg("Writing sample projection");

  ofstream sample_projection_out;
  sample_projection_out.open("coda__projection.tsv");

  if (!sample_projection_out.is_open()) {
    fprintf(stderr,
            "ERROR -- Couldn't open clr_outfile for writing.");
  }

  sample_projection_out << "sample";
  for (size_t i = 0; i < index_to_sample.size(); ++i) {
    sample_projection_out << "\tPC" << (i + 1);
  }
  sample_projection_out << endl;

  for (Eigen::Index i = 0; i < sample_projection.rows(); ++i) {
    string sample = index_to_sample.at(i);
    sample_projection_out << sample
                          << "\t"
                          << sample_projection.row(i).format(TSVFormat)
                          << endl;
  }

  return 0;
}