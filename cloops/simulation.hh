#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <random>

#include "vector.hh"
#include "photon.hh"
#include "geometry.hh"
#include "fibres.hh"

struct Simulation {
    Parameters params;  // simulation parameters
    Detector detector;  // detector geometry
    Fibres fibres;  // fibres
    std::mt19937 rng;  // random number generator
    std::discrete_distribution<> emission_spectrum_distribution;  // emission spectrum distribution for sampling photon wavelengths
    std::discrete_distribution<> scintillation_component_distribution;  // weighted distribution of scintillation time components

    Simulation(
        Parameters parameters,
        Detector detector,
        Fibres fibres,
        unsigned seed
    ) : params(std::move(parameters)),
        detector(std::move(detector)),
        fibres(fibres),
        rng(seed),
        emission_spectrum_distribution(params.emission_spectrum.values.begin(), params.emission_spectrum.values.end()),
        scintillation_component_distribution(params.scint_time_amps.begin(), params.scint_time_amps.end())
        {}

    Photon generate_photon(Vector);
    void track_photon(Photon&, int);
    void propagate_photon(Photon&, double, double);
};

#endif // SIMULATION_H_
