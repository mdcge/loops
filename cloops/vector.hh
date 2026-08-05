#ifndef VECTOR_H_
#define VECTOR_H_

#include <cmath>

// 3D vector
struct Vector {
    double x;
    double y;
    double z;

    Vector() : x(0.0), y(0.0), z(0.0) {}
    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    double mag() const;
    Vector norm() const;
};

double dot(const Vector&, const Vector&);


// Operator overloading
inline Vector operator+(const Vector& v1, const Vector& v2) {
    return Vector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector& operator+=(Vector& v1, const Vector& v2) {
    v1.x += v2.x; v1.y += v2.y; v1.z += v2.z;
    return v1;
}

inline Vector operator-(const Vector& v1, const Vector& v2) {
    return Vector(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Vector& operator-=(Vector& v1, const Vector& v2) {
    v1.x -= v2.x; v1.y -= v2.y; v1.z -= v2.z;
    return v1;
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

inline Vector& operator*=(Vector& v, double s) {
    v.x *= s; v.y *= s; v.z *= s;
    return v;
}


// 2D vector
struct Vector2D {
    double x;
    double y;

    Vector2D() : x(0.0), y(0.0) {}
    Vector2D(double x, double y) : x(x), y(y) {}

    double mag() const;
    Vector2D norm() const;
};

double dot(const Vector2D&, const Vector2D&);


// Operator overloading
inline Vector2D operator+(const Vector2D& v1, const Vector2D& v2) {
    return Vector2D(v1.x + v2.x, v1.y + v2.y);
}

inline Vector2D& operator+=(Vector2D& v1, const Vector2D& v2) {
    v1.x += v2.x; v1.y += v2.y;
    return v1;
}

inline Vector2D operator-(const Vector2D& v1, const Vector2D& v2) {
    return Vector2D(v1.x - v2.x, v1.y - v2.y);
}

inline Vector2D& operator-=(Vector2D& v1, const Vector2D& v2) {
    v1.x -= v2.x; v1.y -= v2.y;
    return v1;
}

inline Vector2D operator-(const Vector2D& v) {
    return Vector2D(-v.x, -v.y);
}

inline Vector2D operator*(double s, const Vector2D& v) {
    return Vector2D(s * v.x, s * v.y);
}

inline Vector2D operator*(const Vector2D& v, double s) {
    return Vector2D(s * v.x, s * v.y);
}

inline Vector2D& operator*=(Vector2D& v, double s) {
    v.x *= s; v.y *= s;
    return v;
}

#endif // VECTOR_H_
