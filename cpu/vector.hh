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

double dot(const Vector& v1, const Vector& v2);

#endif // VECTOR_H_
