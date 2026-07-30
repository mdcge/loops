#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <random>

#include "vector.hh"
#include "photon.hh"

struct Simulation {
    Parameters params;  // simulation parameters
    std::mt19937 rng;  // random number generator
    std::discrete_distribution<> emission_spectrum_distribution;  // emission spectrum distribution for sampling photon wavelengths
    std::discrete_distribution<> scintillation_component_distribution;  // weighted distribution of scintillation time components

    Simulation(
        Parameters parameters,
        unsigned seed
    ) : params(parameters),
        rng(seed),
        emission_spectrum_distribution(params.emission_spectrum.values.begin(), params.emission_spectrum.values.end()),
        scintillation_component_distribution(params.scint_time_amps.begin(), params.scint_time_amps.end())
        {}

    Photon generate_photon(Vector);
    void track_photon(Photon&, int);
    void propagate_photon(Photon&, double, double);
};

#endif // SIMULATION_H_
