#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <variant>

#include "vector.hh"

struct Cylinder {
    Vector O;  // detector origin [mm]
    double r;  // cylinder radius [mm]
    double z;  // cylinder half-height [mm]

    double intersection_distance(const Vector&, const Vector&) const;
};

struct Cuboid {
    Vector O;  // detector origin [mm]
    Vector s;  // cuboid half-sizes [mm]

    double intersection_distance(const Vector&, const Vector&) const;
};

struct Detector {
    std::variant<Cylinder, Cuboid> shape;  // detector shape

    Detector(Vector origin, double radius, double half_height) : shape(Cylinder {origin, radius, half_height}) {}
    Detector(Vector origin, Vector half_sizes) : shape(Cuboid {origin, half_sizes}) {}

    double intersection_distance(const Vector&, const Vector&) const;
};

#endif // GEOMETRY_H_
