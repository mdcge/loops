#include "photon.hh"

#include "parameters.hh"

const double C = 299.8;  // mm/ns
const double C_INV = 1.0 / C;  // ns/mm

// Propagate photon by distance `d` [mm] in refractive index `n`
void Photon::propagate(double d, double n) {
    r += d * p;
    t += d * n * C_INV;
}

// Track a single photon, from creation to absorption
void Photon::track(Parameters& parameters, std::mt19937& rng, int max_steps) {
    // Set random direction
    p = sample_direction(rng);
    
    // Sample absorption length
    abs_dist = sample_length(parameters.absorption_length.at(l), rng);
    double refractive_index = parameters.refractive_index.at(l);

    // Loop until max_steps reached
    for (int i=0; i<max_steps; i++) {
        // Sample scattering length for this step
        double scattering_length = sample_length(parameters.scattering_length.at(l), rng);

        if (abs_dist < scattering_length) {  // if absorption happens first, break out of the loop
            propagate(abs_dist, refractive_index);
            break;
        } else {  // if scattering happens first, propagate the photon and resample direction
            abs_dist -= scattering_length;
            propagate(scattering_length, refractive_index);
            p = sample_direction(rng);
        }
    }
}
