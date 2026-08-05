#include "fibres.hh"

#include <algorithm>

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

// Calculate (sx, sy) spacings of rectangular grid
Vector2D Fibres::spacing_rectangle() {
    double epsilon = 0.0001; // to avoid comparison with 0.0 problems

    // Sort the x & y positions, then take minimum distance between consecutive elements
    std::vector<double> fxs; std::vector<double> fys;
    fxs.reserve(fs.size()); fys.reserve(fs.size());
    for (const Vector2D& f : fs) {
        fxs.push_back(f.x);
        fys.push_back(f.y);
    }
    std::sort(fxs.begin(), fxs.end());
    std::sort(fys.begin(), fys.end());

    double min_spacing_x = std::numeric_limits<double>::infinity();
    double min_spacing_y = std::numeric_limits<double>::infinity();
    for (size_t i=1; i<fs.size(); ++i) {
        double spacing_x = fxs[i] - fxs[i-1];
        double spacing_y = fys[i] - fys[i-1];
        if (spacing_x > epsilon && spacing_x < min_spacing_x) {
            min_spacing_x = spacing_x;
        }
        if (spacing_y > epsilon && spacing_y < min_spacing_y) {
            min_spacing_y = spacing_y;
        }
    }
    return Vector2D(min_spacing_x, min_spacing_y);
}

// Calculate s spacing of hexagonal grid
double Fibres::spacing_hexagon() {
    double epsilon = 0.0001; // to avoid comparison with 0.0 problems

    // Compare distances between each pair of points, take the minimum
    double min_dist2 = std::numeric_limits<double>::infinity();
    for (size_t i=0; i<fs.size(); i++) {
        for (size_t j=0; j<fs.size(); j++) {
            Vector2D df = fs[i] - fs[j];
            double dist2 = df.x * df.x + df.y * df.y;
            if (dist2 > epsilon && dist2 < min_dist2) {
                min_dist2 = dist2;
            }
        }
    }
    return std::sqrt(min_dist2);
}
