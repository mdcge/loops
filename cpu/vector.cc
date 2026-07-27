#include "vector.hh"

double Vector::mag() const {
    return std::sqrt(x * x + y * y + z * z);
}

Vector Vector::norm() const {
    double m = mag();
    if (m == 0.0) { return Vector(0.0, 0.0, 0.0); }
    double mag_inv = 1.0 / m;
    return Vector(x * mag_inv, y * mag_inv, z * mag_inv);
}

double dot(const Vector& v1, const Vector& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
