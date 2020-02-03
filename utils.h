#ifndef CODA_UTILS_H
#define CODA_UTILS_H

#include <array>
#include <string>
#include <iostream>
#include <Eigen/Dense>

void
log_msg(const char* msg);


std::vector<std::string>
split(const std::string& s, char delimiter);

void
split(std::array<std::string, 4>& tokens, const std::string& s, char delimiter);


#endif //CODA_UTILS_H
