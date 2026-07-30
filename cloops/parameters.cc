#include "parameters.hh"

#include <cmath>

// === OpticalProperty ===

double OpticalProperty::at(double wavelength) const {
    // Pointer to first wavelength >= value being searched
    std::vector<double>::const_iterator it = std::lower_bound(wavelengths.begin(), wavelengths.end(), wavelength);

    if (it == wavelengths.begin()) return values.front();  // if wavelength < min(wavelengths), return `values[0]`
    if (it == wavelengths.end())   return values.back();  // if wavelength > max(wavelengths), return `values[-1]`

    // If exact match, just return corresponding value
    size_t i = std::distance(wavelengths.begin(), it);
    if (wavelengths[i] == wavelength) return values[i];

    // Otherwise, interpolate
    double x0 = wavelengths[i-1], x1 = wavelengths[i];
    double y0 = values[i-1],     y1 = values[i];
    return y0 + (wavelength - x0) * (y1 - y0) / (x1 - x0);
}

// === Parameters ===

double Parameters::sample_scattering_length(double wavelength, std::mt19937& rng)  {
    double scat_len = scattering_length.at(wavelength);
    return sample_length(scat_len, rng);
}

double Parameters::sample_absorption_length(double wavelength, std::mt19937& rng)  {
    double abs_len = absorption_length.at(wavelength);
    return sample_length(abs_len, rng);
}

double Parameters::lookup_refractive_index(double wavelength) {
    return refractive_index.at(wavelength);
}

// === General ===

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
