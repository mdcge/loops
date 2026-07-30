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
    double wl;  // wavelength [nm]
    
    void propagate(double, double);    
    void track(Parameters&, std::mt19937&, int);
};

#endif // PHOTON_H_
