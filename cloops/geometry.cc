#include "geometry.hh"

#include <cmath>
#include <algorithm>
#include <limits>

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
        double discriminant = b*b - 4*a*c;
        if (discriminant >= 0.0) { // if discriminant gives valid solutions
            t_side = (-b + std::sqrt(discriminant)) / (2*a);
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

// Get detector origin
Vector Detector::O() const {
    return std::visit([](const auto& s) -> const Vector& {
        return s.O;
    }, shape);
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

// Calculate intersection distance from `r` to circle of centre `C` and radius `cr` along `p`
double circle_intersection_distance(const Vector2D& r, const Vector2D& p, const Vector2D& C, double cr) {
    Vector2D point = r - C;  // "local" point accounting for circle centre
    // Solve:
    //     |r + t * p|^2 = cr^2
    // =>  r^2 + p^2 * t^2 + 2 * r * p * t = cr^2
    // =>  p^2 * t^2  +  2 * r * p * t  +  r^2 - cr^2  =  0
    // =>   a  * t^2  +      b     * t  +      c       =  0
    // =>  t = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
    // for t
    double a = p.x * p.x + p.y * p.y;
    double b = 2 * (point.x * p.x + point.y * p.y);
    double c = (point.x * point.x + point.y * point.y) - cr * cr;
    double discriminant = b*b - 4*a*c;
    if (discriminant < 0.0) { // if discriminant not valid, return infinity
        return std::numeric_limits<double>::infinity();
    } else { // if discriminant is valid, pick `t` value corresponding to entering distance (assume r outside circle)
        double t = (-b - std::sqrt(discriminant)) / (2 * a);
        if (t < 0.0) {
            return std::numeric_limits<double>::infinity();
        }
        return t;
    }
}
