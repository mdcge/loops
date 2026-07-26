#ifndef VECTOR_H_
#define VECTOR_H_

struct Vector {
    double x;
    double y;
    double z;

    Vector(double x, double y, double z) : x(x), y(y), z(z) {}
};

double dot(const Vector& v1, const Vector& v2);

#endif // VECTOR_H_
