#include "geometry.hh"

#include <cmath>
#include <algorithm>

#include "vector.hh"

// Calculate distance from point `x` to cylinder along direction `d`
double Cylinder::intersection_distance(const Vector& x, const Vector& d) const {
    Vector point = x - O;  // "local" point accounting for detector origin

    // Side wall: solve 2D case
    //     |point.xy + t * d.xy|^2 = r^2
    // =>  point.xy^2 + d.xy^2 * t^2 + 2 * point.xy * d.xy * t = r^2
    // =>  d.xy^2 * t^2  +  2 * point.xy * d.xy * t  +  point.xy^2 - r^2  =  0
    // =>     a   * t^2  +            b         * t  +            c       =  0
    // =>  t = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
    // for t
    double t_side = std::numeric_limits<double>::infinity();
    double a = 1.0 - d.z*d.z;  // = d.x*d.x + d.y*d.y because d is unit
    if (a > 0.0) { // if a == 0, it's a purely axial photon: skip side wall check
        double b = 2.0 * (point.x*d.x + point.y*d.y);
        double c = point.x*point.x + point.y*point.y - r*r;
        double discriminator = b*b - 4*a*c;
        if (discriminator >= 0.0) { // if discriminator gives valid solutions
            t_side = (-b + std::sqrt(discriminator)) / (2*a);
        }
    }

    // Caps: solve
    // point.z + t * d.z = +-z
    double t_cap = std::numeric_limits<double>::infinity();
    if (d.z != 0.0) { // if d.z == 0, it's a purely radial photon: skip cap check
        double cap_z = std::copysign(z, d.z);  // +z if d -> +z ,  -z if d -> -z
        t_cap = (cap_z - point.z) / d.z;  // distance to cap
    }

    return std::min(t_side, t_cap);
}

// Calculate distance from point `x` to cuboid along direction `d`
double Cuboid::intersection_distance(const Vector& x, const Vector& d) const {
    Vector point = x - O;  // "local" point accounting for detector origin

    // Closure for calculating distance to each each "slab", i.e. axis plane
    auto slab = [](double point_i, double d_i, double s_i) { // for each axis `i`
        if (d_i == 0.0) { // direction parallel to this axis plane
            return std::numeric_limits<double>::infinity();
        }
        double face = std::copysign(s_i, d_i);  // determine which (+/-) slab direction is pointing at
        return (face - point_i) / d_i;
    };

    return std::min({slab(point.x, d.x, s.x), slab(point.y, d.y, s.y), slab(point.z, d.z, s.z)});
}

// Dispatch intersection calculation to relevant shape
double Detector::intersection_distance(const Vector& x, const Vector& d) const {
    if (auto* cylinder = std::get_if<Cylinder>(&shape)) {
        return cylinder -> intersection_distance(x, d);
    } else {
        const Cuboid& cuboid = std::get<Cuboid>(shape);
        return cuboid.intersection_distance(x, d);
    }
}
