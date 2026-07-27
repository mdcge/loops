#ifndef PHOTON_H_
#define PHOTON_H_

#include "vector.hh"

struct Photon {
    Vector r;  // position [mm]
    Vector p;  // direction (normalised)
    double t;  // time [ns]

    Photon(Vector r, Vector p, double t) : r(r), p(p), t(t) {}

    void propagate(double);
};

#endif // PHOTON_H_
