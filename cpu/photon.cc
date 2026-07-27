#include "photon.hh"

const double C = 299.8;  // mm/ns

// Propagate photon by distance `d` [mm]
void Photon::propagate(double d) {
    r += d * p;
    t += d / C;
}
