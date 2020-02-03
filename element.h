#ifndef CODA_ELEMENT_H
#define CODA_ELEMENT_H

#include <string>
#include <array>

#define zero_replacement 0.05f

typedef struct Element
{
  std::string sample;
  std::string otu;
  size_t sample_index;
  size_t otu_index;
  float  otu_len;
  float  count;
} Element;

void
element_print(const Element& element);

Element
parse_data_line(const std::array<std::string, 4>& tokens, bool replace_zeros);

/// Counts/hits per 100 AAs/bases
float
element_normalize_count(const Element& element);

#endif //CODA_ELEMENT_H
