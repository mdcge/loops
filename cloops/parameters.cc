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

// === General ===

double sample_length(double lambda, xoshiro256pp& rng)  {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double u = uniform(rng);
    return -lambda * std::log(u);
}

// Sample direction isotropically
Vector sample_direction(xoshiro256pp& rng) {
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    double x1, x2, s;
    do {
        x1 = uniform(rng);
        x2 = uniform(rng);
        s = x1*x1 + x2*x2;
    } while (s >= 1.0);   // rejection: keep points inside the unit disk (~78.5% acceptance)

    double factor = 2.0 * std::sqrt(1.0 - s);
    return Vector(
        x1 * factor,
        x2 * factor,
        1.0 - 2.0*s
    );
}
