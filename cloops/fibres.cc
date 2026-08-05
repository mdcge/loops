#include "fibres.hh"

Vector2D Fibres::f00(Vector2D O) const {
    double min_dist2 = std::numeric_limits<double>::infinity();
    Vector2D foo;
    for (const Vector2D& f : fs) {
        double dist2 = (f-O).x * (f-O).x + (f-O).y * (f-O).y;
        if (dist2 < min_dist2) {
            min_dist2 = dist2;
            foo = f;
        }
    }
    return foo;
}
