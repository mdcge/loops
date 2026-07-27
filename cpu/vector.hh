#ifndef VECTOR_H_
#define VECTOR_H_

#include <cmath>

struct Vector {
    double x;
    double y;
    double z;

    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    double mag() const;
    Vector norm() const;
};

double dot(const Vector&, const Vector&);

// Operator overloading
inline Vector operator+(const Vector& v1, const Vector& v2) {
    return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector operator-(const Vector& v1, const Vector& v2) {
    return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Vector operator-(const Vector& v) {
    return Vector(-v.x, -v.y, -v.z);
}

inline Vector operator*(double s, const Vector& v) {
    return Vector(s * v.x, s * v.y, s * v.z);
}

inline Vector operator*(const Vector& v, double s) {
    return Vector(s * v.x, s * v.y, s * v.z);
}

#endif // VECTOR_H_
