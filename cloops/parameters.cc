#include "parameters.hh"

#include <cmath>

double Parameters::sample_scattering_length(double wavelength, std::mt19937& rng)  {
    double scat_len = lookup_value(scattering_length, wavelength);
    return sample_length(scat_len, rng);
}

double Parameters::sample_absorption_length(double wavelength, std::mt19937& rng)  {
    double abs_len = lookup_value(absorption_length, wavelength);
    return sample_length(abs_len, rng);
}

double Parameters::lookup_refractive_index(double wavelength) {
    return lookup_value(refractive_index, wavelength);
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

// Extract value from property (xs and ys), given input x
// Interpolate if needed
// If x not in [min(xs), max(xs)], just take the boundary value of ys
double lookup_value(const OpticalProperty& property, double x) {
    const std::vector<double>& xs = property.wavelengths;
    const std::vector<double>& ys = property.values;
    
    // Pointer to first element >= value being searched
    std::vector<double>::const_iterator it = std::lower_bound(xs.begin(), xs.end(), x);

    // If x < min(xs), return ys[min(xs)]
    if (it == xs.begin())
        return ys.front();

    // If x > max(xs), return ys[max(xs)]
    if (it == xs.end())
        return ys.back();

    // Index of desired value
    size_t i = std::distance(xs.begin(), it);

    // If match is exact, just return corresponding y
    if (xs[i] == x)
        return ys[i];

    // Otherwise, interpolate
    double x0 = xs[i - 1], x1 = xs[i];
    double y0 = ys[i - 1], y1 = ys[i];

    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}
