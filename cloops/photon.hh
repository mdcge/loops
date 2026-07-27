#ifndef PHOTON_H_
#define PHOTON_H_

#include "vector.hh"

#include "parameters.hh"

#include <limits>

struct Photon {
    Vector r;  // position [mm]
    Vector p;  // direction (normalised)
    double t;  // time [ns]
    double abs_dist;  // distance to absorption [mm]

    Photon(Vector r, Vector p, double t) : r(r), p(p), t(t), abs_dist(std::numeric_limits<double>::infinity()) {}
    Photon(Vector r, double t) : r(r), p(Vector(0.0, 0.0, 0.0)), t(t), abs_dist(std::numeric_limits<double>::infinity()) {}

    void propagate(double, double);
    
    void track(Parameters&, double, std::mt19937&, int);
};

#endif // PHOTON_H_
