#include "photon.hh"

const double C = 299.8;  // mm/ns
const double C_INV = 1.0 / C;  // ns/mm

// Propagate photon by distance `d` [mm]
void Photon::propagate(double d) {
    r += d * p;
    t += d * C_INV;
}

// Track a single photon, from creation to absorption
void Photon::track(Parameters& parameters, std::mt19937& rng, int max_steps) {
    // Set random direction
    p = sample_direction(rng);
    
    // Sample absorption length
    abs_dist = parameters.sample_absorption_length(rng);

    // Loop until max_steps reached
    for (int i=0; i<max_steps; i++) {
        // Sample scattering length for this step
        double scattering_length = parameters.sample_scattering_length(rng);

        if (abs_dist < scattering_length) {  // if absorption happens first, break out of the loop
            propagate(abs_dist);
            break;
        } else {  // if scattering happens first, propagate the photon and resample direction
            abs_dist -= scattering_length;
            propagate(scattering_length);
            p = sample_direction(rng);
        }
    }
}
