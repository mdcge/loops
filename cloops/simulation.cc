#include "simulation.hh"

const double C = 299.8;  // mm/ns
const double C_INV = 1.0 / C;  // ns/mm

Photon Simulation::generate_photon(Vector origin) {
    Photon photon;
    photon.r = origin;
    photon.p = sample_direction(rng);
    photon.wl = params.emission_spectrum.wavelengths[emission_spectrum_distribution(rng)];
    
    int time_component = scintillation_component_distribution(rng);
    std::exponential_distribution<double> rise(1.0 / params.scint_rise_time);
    std::exponential_distribution<double> decay(-1.0 / params.scint_time_cnsts[time_component]);
    photon.t = rise(rng) + decay(rng);
    
    return photon;
}

// Propagate photon by distance `d` [mm] in refractive index `n`
void Simulation::propagate_photon(Photon& photon, double d, double n) {
    photon.r += d * photon.p;
    photon.t += d * n * C_INV;
}

// Track a single photon, from creation to absorption
void Simulation::track_photon(Photon& photon, int max_steps) {
    // Set random direction
    photon.p = sample_direction(rng);
    
    // Sample absorption length
    photon.abs_dist = sample_length(params.absorption_length.at(photon.wl), rng);
    double refractive_index = params.refractive_index.at(photon.wl);

    // Loop until max_steps reached
    for (int i=0; i<max_steps; i++) {
        // Sample scattering length for this step
        double scattering_length = sample_length(params.scattering_length.at(photon.wl), rng);

        if (photon.abs_dist < scattering_length) {  // if absorption happens first, break out of the loop
            propagate_photon(photon, photon.abs_dist, refractive_index);
            break;
        } else {  // if scattering happens first, propagate the photon and resample direction
            photon.abs_dist -= scattering_length;
            propagate_photon(photon, scattering_length, refractive_index);
            photon.p = sample_direction(rng);
        }
    }
}
