#include "fibres.hh"

#include <algorithm>
#include <cmath>

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
Vector2D Fibres::spacing_rectangle() const {
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
double Fibres::spacing_hexagon() const {
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

// Calculates the rectangular cell index (i, j) in which point `r` lies in a detector with origin `O` with grid spacing (`sx`, `sy`), given a closest fibre to origin `foo`
std::pair<int, int> cell_index_rectangle(const Vector2D& r, const Vector2D& O, const Vector2D& s, const Vector2D& foo) {
    Vector2D alpha((O.x - foo.x) / s.x, (O.y - foo.y) / s.y);
    int i = std::floor((r.x - O.x) / s.x + alpha.x + 0.5);
    int j = std::floor((r.y - O.y) / s.y + alpha.y + 0.5);
    return std::make_pair(i, j);
}

// Calculates the hexagonal cell index (i, j) in which point `r` lies in a detector with origin `O` with grid spacing `s`, given a closest fibre to origin `foo`
std::pair<int, int> cell_index_hexagon(const Vector2D& r, const Vector2D& O, double s, const Vector2D& foo) {
    double s_inv = 1/s;
    double sqrt3_inv = 1/std::sqrt(3);
    
    Vector2D alpha(O.x - foo.x, O.y - foo.y);
    // Calculate cube coordinates lying on x + y + z = 0
    double xc = (r.x - O.x + alpha.x) * s_inv - (r.y - O.y + alpha.y) * s_inv * sqrt3_inv;
    double yc = 2 * (r.y - O.y + alpha.y) * s_inv * sqrt3_inv;
    double zc = -(xc + yc);

    // Rounded cube coordinates
    int qc = std::round(xc);
    int rc = std::round(yc);
    int sc = std::round(zc);

    // Round amounts for each coordinate: largest round amount is closest to hexagon boundary
    double dxc = std::abs(qc - xc);
    double dyc = std::abs(rc - yc);
    double dzc = std::abs(sc - zc);

    // Bump coordinate closest to boundary back on the x+y+z=0 plane
    if (dxc > dyc && dxc > dzc) {
        qc = -(rc + sc);
    } else if (dyc > dxc && dyc > dzc) {
        rc = -(qc + sc);
    } // no need for `sc = -(qc + rc)` as it's never used
    return std::make_pair(qc, rc);
}
