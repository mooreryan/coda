#include "utils.h"

using namespace std;

void
log_msg(const char* msg)
{
  fprintf(stderr, "LOG -- %s\n", msg);
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
