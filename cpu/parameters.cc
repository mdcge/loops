#include "parameters.hh"

#include <random>
#include <cmath>

double Parameters::sample_scattering_length(std::mt19937& rng)  {
    return sample_length(ls, rng);
}

double Parameters::sample_absorption_length(std::mt19937& rng)  {
    return sample_length(la, rng);
}

double sample_length(double lambda, std::mt19937& rng)  {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double u = uniform(rng);
    return -lambda * std::log(u);
}
