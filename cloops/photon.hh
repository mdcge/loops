#ifndef PHOTON_H_
#define PHOTON_H_

#include "vector.hh"

#include "parameters.hh"

#include <limits>
#include <random>

struct Photon {
    Vector r;  // position [mm]
    Vector p;  // direction (normalised)
    double t;  // time [ns]
    double l;  // wavelength [nm]
    double abs_dist;  // distance to absorption [mm]

    Photon(Vector r, Vector p, double t, double l) : r(r), p(p), t(t), l(l), abs_dist(std::numeric_limits<double>::infinity()) {}
    Photon(Vector r, double t, double l) : r(r), p(Vector(0.0, 0.0, 0.0)), t(t), l(l), abs_dist(std::numeric_limits<double>::infinity()) {}
    Photon(Vector r, Parameters& parameters, std::mt19937& rng) : r(r), p(Vector(0.0, 0.0, 0.0)), abs_dist(std::numeric_limits<double>::infinity()) {
        std::discrete_distribution<int> component(parameters.scint_time_amps.begin(), parameters.scint_time_amps.end());
        int i = component(rng);

        std::exponential_distribution<double> rise(1.0 / parameters.scint_rise_time);
        std::exponential_distribution<double> decay(-1.0 / parameters.scint_time_cnsts[i]);
        
        t = rise(rng) + decay(rng);

        std::discrete_distribution<> dist(parameters.emission_spectrum.values.begin(), parameters.emission_spectrum.values.end());
        l = parameters.emission_spectrum.wavelengths[dist(rng)];
    }

    void propagate(double, double);
    
    void track(Parameters&, std::mt19937&, int);
};

#endif // PHOTON_H_
