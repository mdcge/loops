#include "parameters.hh"

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

// Sample direction isotropically
Vector sample_direction(std::mt19937& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    double cos_theta = 2.0 * uniform(rng) - 1.0;  // uniform in [-1, 1]
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = 2.0 * M_PI * uniform(rng);  // uniform in [0, 2pi]

    return Vector(
        sin_theta * std::cos(phi),
        sin_theta * std::sin(phi),
        cos_theta
    );
}
