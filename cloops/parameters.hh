#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <random>

#include "vector.hh"

struct Parameters {
    double ls;  // scattering length
    double la;  // absorption length

    Parameters(double scattering_length, double absorption_length) : ls(scattering_length), la(absorption_length) {}

    double sample_scattering_length(std::mt19937&);
    double sample_absorption_length(std::mt19937&);
};

double sample_length(double, std::mt19937&);
Vector sample_direction(std::mt19937&);

#endif // PARAMETERS_H_
