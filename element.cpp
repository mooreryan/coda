#include "element.h"

using namespace std;

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
parse_data_line(const std::array<string, 4>& tokens, bool replace_zeros)
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
  if (replace_zeros && count == 0) {
    element.count = zero_replacement;
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
