#include "fibres.hh"

#include <algorithm>
#include <cmath>
#include<limits>
#include <stdexcept>

// Infer the type of lattice given by the fibre positions: rectangle or hexagon
LatticeType infer_lattice_type(const std::vector<Vector2D>& fs) {
    double epsilon = 0.001;

    // Smallest vector with only horizontal component 
    double hmin = std::numeric_limits<double>::infinity();
    Vector2D hvec;
    // Smallest vector with non-zero vertical component 
    double vmin = std::numeric_limits<double>::infinity();
    Vector2D vvec;

    for (size_t i=0; i<fs.size(); ++i) {
        for (size_t j=0; j<fs.size(); ++j) {
            if (i == j) { continue; }

            Vector2D d = fs[i] - fs[j];
            double dist = d.mag();
            if (std::abs(d.y) < epsilon) { // is vector horizontal-only?
                if (dist < hmin) {
                    hmin = dist;
                    hvec = d;
                }
            } else { // does vector have some vertical component?
                if (dist < vmin) {
                    vmin = dist;
                    vvec = d;
                }
            }
        }
    }

    // Angle between shortest horizontal and part-vertical vectors: 0 if rectangular, >0 otherwise
    double costheta = dot(hvec, vvec) / (hmin * vmin);
    return std::abs(costheta) < 0.1 ? LatticeType::Rectangular : LatticeType::Hexagonal;
}

// Closest fibre to origin
Vector2D Fibres::f00(const Vector2D& O) const {
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

// Location of fibre in cell (i, j) in a rectangular grid
Vector2D Fibres::fij_rectangle(const Vector2D& foo, int i, int j, const Vector2D& s) const {
    return Vector2D(foo.x + i * s.x, foo.y + j * s.y);
}

// Location of fibre in cell (i, j) in a hexagonal grid
Vector2D Fibres::fij_hexagon(const Vector2D& foo, int i, int j, double s) const {
    return Vector2D(foo.x + (i + j/2.0) * s, foo.y + (j * std::sqrt(3)/2.0) * s);
}

Vector2D Fibres::fij(const Vector2D& foo, int i, int j) {
    if (lattice_type == LatticeType::Rectangular) {
        return fij_rectangle(foo, i, j, std::get<Vector2D>(s));
    } else {
        return fij_hexagon(foo, i, j, std::get<double>(s));
    }
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

// Populate the fibre cell set with (i, j) cells which contain fibres
void Fibres::build_fibre_cells_rectangle(const Vector2D& foo, const Vector2D& s) {
    fibre_cells.clear();
    for (const Vector2D& f : fs) {
        fibre_cells.insert(cell_index_rectangle(f, s, foo));
    }
}

// Populate the fibre cell set with (i, j) cells which contain fibres
void Fibres::build_fibre_cells_hexagon(const Vector2D& foo, double s) {
    fibre_cells.clear();
    for (const Vector2D& f : fs) {
        fibre_cells.insert(cell_index_hexagon(f, s, foo));
    }
}

void Fibres::build_fibre_cells() {
    if (lattice_type == LatticeType::Rectangular) {
        build_fibre_cells_rectangle(foo, std::get<Vector2D>(s));
    } else {
        build_fibre_cells_hexagon(foo, std::get<double>(s));
    }
}

// Given a cell (i, j) in a rectangular grid, calculate the tx and ty distances along the photon direction `php` from the photon position `phr` needed to reach the x & y cell boundaries
std::array<double, 3> cell_ts_rectangle(const std::pair<int, int> ij, const Vector2D& phr, const Vector2D& php, const Vector2D& foo, const Vector2D& s) {
    int i = ij.first;
    int j = ij.second;
    double fx = foo.x + i * s.x;
    double fy = foo.y + j * s.y;

    double tx;
    if (php.x > 0) {
        tx = (fx + s.x/2 - phr.x) / php.x;
    } else if (php.x < 0) {
        tx = (fx - s.x/2 - phr.x) / php.x;
    } else {
        tx = std::numeric_limits<double>::infinity();
    }
    
    double ty;
    if (php.y > 0) {
        ty = (fy + s.y/2 - phr.y) / php.y;
    } else if (php.y < 0) {
        ty = (fy - s.y/2 - phr.y) / php.y;
    } else {
        ty = std::numeric_limits<double>::infinity();
    }

    return {tx, ty, std::numeric_limits<double>::infinity()};
}

// Given a cell (i, j) in a hexagonal grid, calculate the tx and ty distances along the photon direction `php` from the photon position `phr` needed to reach the x & y cell boundaries
std::array<double, 3> cell_ts_hexagon(const std::pair<int, int> ij, const Vector2D& phr, const Vector2D& php, const Vector2D& foo, double s) {
    static const double sqrt3 = std::sqrt(3);
    static const Vector2D n1(1.0, 0.0);
    static const Vector2D n2(0.5, sqrt3/2.0);
    static const Vector2D n3(-0.5, sqrt3/2.0);
    int i = ij.first; int j = ij.second;
    Vector2D f(foo.x + (i + j/2.0) * s, foo.y + (j * sqrt3/2.0) * s);

    double t1;
    if (dot(n1, php) > 0) {
        t1 = (dot(n1, f) + s/2 - dot(n1, phr)) / dot(n1, php);
    } else if (dot(n1, php) < 0) {
        t1 = (dot(n1, f) - s/2 - dot(n1, phr)) / dot(n1, php);
    } else {
        t1 = std::numeric_limits<double>::infinity();
    }
    
    double t2;
    if (dot(n2, php) > 0) {
        t2 = (dot(n2, f) + s/2 - dot(n2, phr)) / dot(n2, php);
    } else if (dot(n2, php) < 0) {
        t2 = (dot(n2, f) - s/2 - dot(n2, phr)) / dot(n2, php);
    } else {
        t2 = std::numeric_limits<double>::infinity();
    }
    
    double t3;
    if (dot(n3, php) > 0) {
        t3 = (dot(n3, f) + s/2 - dot(n3, phr)) / dot(n3, php);
    } else if (dot(n3, php) < 0) {
        t3 = (dot(n3, f) - s/2 - dot(n3, phr)) / dot(n3, php);
    } else {
        t3 = std::numeric_limits<double>::infinity();
    }

    return {t1, t2, t3};
}

std::array<double, 3> cell_ts(const std::pair<int, int> ij, const Vector2D& phr, const Vector2D& php, const Vector2D& foo, const Spacing& s, LatticeType lattice_type) {
    if (lattice_type == LatticeType::Rectangular) {
        return cell_ts_rectangle(ij, phr, php, foo, std::get<Vector2D>(s));
    } else {
        return cell_ts_hexagon(ij, phr, php, foo, std::get<double>(s));
    }
}

// Calculate the new cell index in a rectangular grid based on the photon position and direction using the pre-computed (tx, ty) values
std::pair<int, int> new_cell_index_rectangle(const std::pair<int, int> ij, const std::array<double, 3>& ts, const Vector2D& php) {
    int i = ij.first; int j = ij.second;
    double tx = ts[0]; double ty = ts[1];

    if (tx <= ty) {
        return std::make_pair(i + std::copysign(1, php.x), j);
    } else {
        return std::make_pair(i, j + std::copysign(1, php.y));
    }
}

// Calculate the new cell index in a hexagonal grid based on the photon position and direction using the pre-computed (tx, ty) values
std::pair<int, int> new_cell_index_hexagon(const std::pair<int, int> ij, const std::array<double, 3>& ts, const Vector2D& php) {
    static const double sqrt3 = std::sqrt(3);
    static const Vector2D n1(1.0, 0.0);
    static const Vector2D n2(0.5, sqrt3/2.0);
    static const Vector2D n3(-0.5, sqrt3/2.0);
    int i = ij.first; int j = ij.second;
    double t1 = ts[0]; double t2 = ts[1]; double t3 = ts[2];

    if (t1 <= t2 && t1 <= t3) {
        return std::make_pair(i + std::copysign(1, dot(n1, php)), j);
    } else if (t2 <= t1 && t2 <= t3) {
        return std::make_pair(i, j + std::copysign(1, dot(n2, php)));
    } else {
        return std::make_pair(i - std::copysign(1, dot(n3, php)), j + std::copysign(1, dot(n3, php)));
    }
}

std::pair<int, int> new_cell_index(const std::pair<int, int> ij, const std::array<double, 3>& ts, const Vector2D& php, LatticeType lattice_type) {
    if (lattice_type == LatticeType::Rectangular) {
        return new_cell_index_rectangle(ij, ts, php);
    } else {
        return new_cell_index_hexagon(ij, ts, php);
    }
}

// Calculate the new photon position in a rectangular grid after stepping to the nearest cell along the photon direction
Vector2D new_photon_position_rectangle(const std::array<double, 3>& ts, const Vector2D& phr, const Vector2D& php) {
    double tx = ts[0]; double ty = ts[1];
    return phr + std::min(tx, ty) * php;
}

// Calculate the new photon position in a hexagonal grid after stepping to the nearest cell along the photon direction
Vector2D new_photon_position_hexagon(const std::array<double, 3>& ts, const Vector2D& phr, const Vector2D& php) {
    double t1 = ts[0]; double t2 = ts[1]; double t3 = ts[2];
    return phr + std::min({t1, t2, t3}) * php;
}

Vector2D new_photon_position(const std::array<double, 3>& ts, const Vector2D& phr, const Vector2D& php, LatticeType lattice_type) {
    if (lattice_type == LatticeType::Rectangular) {
        return new_photon_position_rectangle(ts, phr, php);
    } else {
        return new_photon_position_hexagon(ts, phr, php);
    }
}

// Calculates the rectangular cell index (i, j) in which point `r` lies in a detector with origin `O` with grid spacing (`sx`, `sy`), given a closest fibre to origin `foo`
std::pair<int, int> cell_index_rectangle(const Vector2D& r, const Vector2D& s, const Vector2D& foo) {
    int i = std::floor((r.x - foo.x) / s.x + 0.5);
    int j = std::floor((r.y - foo.y) / s.y + 0.5);
    return std::make_pair(i, j);
}

// Calculates the hexagonal cell index (i, j) in which point `r` lies in a detector with origin `O` with grid spacing `s`, given a closest fibre to origin `foo`
std::pair<int, int> cell_index_hexagon(const Vector2D& r, double s, const Vector2D& foo) {
    double s_inv = 1/s;
    double sqrt3_inv = 1/std::sqrt(3);
    
    // Calculate cube coordinates lying on x + y + z = 0
    double xc = (r.x - foo.x) * s_inv - (r.y - foo.y) * s_inv * sqrt3_inv;
    double yc = 2 * (r.y - foo.y) * s_inv * sqrt3_inv;
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

std::pair<int, int> cell_index(const Vector2D& r, const Spacing& s, const Vector2D& foo, LatticeType lattice_type) {
    if (lattice_type == LatticeType::Rectangular) {
        return cell_index_rectangle(r, std::get<Vector2D>(s), foo);
    } else {
        return cell_index_hexagon(r, std::get<double>(s), foo);
    }
}


// DDA state functions

RectDDAState make_dda_state_rectangle(const Vector2D& php, const Vector2D& s) {
    RectDDAState state;
    // Inverse dot products with each direction eigenvector
    state.inv_dotx = (php.x != 0.0) ? 1.0 / php.x : 0.0;
    state.inv_doty = (php.y != 0.0) ? 1.0 / php.y : 0.0;
    // Cell spans in each direction
    state.dtx = (php.x != 0.0) ? s.x * std::abs(state.inv_dotx) : std::numeric_limits<double>::infinity();
    state.dty = (php.y != 0.0) ? s.y * std::abs(state.inv_doty) : std::numeric_limits<double>::infinity();
    // Cell index changes in each direction
    state.di = static_cast<int>(std::copysign(1.0, php.x));
    state.dj = static_cast<int>(std::copysign(1.0, php.y));

    return state;
}

HexDDAState make_dda_state_hexagon(const Vector2D& php, double s) {
    // Initialise constants and fixed vectors
    double sqrt3 = std::sqrt(3);
    Vector2D n1(1.0, 0.0), n2(0.5, sqrt3/2.0), n3(-0.5, sqrt3/2.0);
    double dot1 = dot(n1, php); double dot2 = dot(n2, php); double dot3 = dot(n3, php);

    HexDDAState state;
    // Inverse dot products with each direction eigenvector
    state.inv_dot1 = (dot1 != 0.0) ? 1.0 / dot1 : 0.0;
    state.inv_dot2 = (dot2 != 0.0) ? 1.0 / dot2 : 0.0;
    state.inv_dot3 = (dot3 != 0.0) ? 1.0 / dot3 : 0.0;
    // Cell spans in each direction
    state.dt1 = (dot1 != 0.0) ? s * std::abs(state.inv_dot1) : std::numeric_limits<double>::infinity();
    state.dt2 = (dot2 != 0.0) ? s * std::abs(state.inv_dot2) : std::numeric_limits<double>::infinity();
    state.dt3 = (dot3 != 0.0) ? s * std::abs(state.inv_dot3) : std::numeric_limits<double>::infinity();
    // Cell index changes in each direction
    state.di1 = static_cast<int>(std::copysign(1.0, dot1)); state.dj1 = 0; // i += 1
    state.di2 = 0; state.dj2 = static_cast<int>(std::copysign(1.0, dot2)); // j += 1
    state.di3 = -static_cast<int>(std::copysign(1.0, dot3)); state.dj3 = static_cast<int>(std::copysign(1.0, dot3)); // i -= 1, j += 1
    // Boundary signs in each direction
    state.bsign1 = static_cast<int>(std::copysign(1.0, dot1));
    state.bsign2 = static_cast<int>(std::copysign(1.0, dot2));
    state.bsign3 = static_cast<int>(std::copysign(1.0, dot3)); // !! `di3` = -`bsign3` !!

    return state;
}

// Initialise the DDA state values: calculating the full-cell spans and the cell index changes in each direction
// Based on photon direction `php` and fibre spacing `s`
DDAState make_dda_state(const Vector2D& php, const Spacing& s, LatticeType lattice_type) {
    if (lattice_type == LatticeType::Rectangular) {
        return make_dda_state_rectangle(php, std::get<Vector2D>(s));
    } else {
        return make_dda_state_hexagon(php, std::get<double>(s));
    }
}

std::array<double, 3> initial_ts_rectangle(const std::pair<int, int>& ij, const Vector2D& phr, const Vector2D& foo, const Vector2D& s, const RectDDAState& state) {
    // Pre-calculate necessary constants
    int i = ij.first; int j = ij.second;
    double fx = foo.x + i * s.x; double fy = foo.y + j * s.y;

    // Distance to each cell boundary
    double tx = (state.inv_dotx != 0.0) ? (fx + (s.x*state.di)/2.0 - phr.x) * state.inv_dotx : std::numeric_limits<double>::infinity();
    double ty = (state.inv_doty != 0.0) ? (fy + (s.y*state.dj)/2.0 - phr.y) * state.inv_doty : std::numeric_limits<double>::infinity();

    return {tx, ty, std::numeric_limits<double>::infinity()};
}

std::array<double, 3> initial_ts_hexagon(const std::pair<int, int>& ij, const Vector2D& phr, const Vector2D& foo, double s, const HexDDAState& state) {
    // Initialise constants and fixed vectors
    int i = ij.first; int j = ij.second;
    double sqrt3 = std::sqrt(3);
    Vector2D n1(1.0, 0.0), n2(0.5, sqrt3/2.0), n3(-0.5, sqrt3/2.0);
    Vector2D f(foo.x + (i + j/2.0) * s, foo.y + (j * sqrt3/2.0) * s);

    // Distance to each cell boundary
    double t1 = (state.inv_dot1 != 0.0) ? (dot(n1, f) + (s*state.bsign1)/2.0 - dot(n1, phr)) * state.inv_dot1 : std::numeric_limits<double>::infinity();
    double t2 = (state.inv_dot2 != 0.0) ? (dot(n2, f) + (s*state.bsign2)/2.0 - dot(n2, phr)) * state.inv_dot2 : std::numeric_limits<double>::infinity();
    double t3 = (state.inv_dot3 != 0.0) ? (dot(n3, f) + (s*state.bsign3)/2.0 - dot(n3, phr)) * state.inv_dot3 : std::numeric_limits<double>::infinity();
    
    return {t1, t2, t3};
}

// Find the distance to each cell boundary for the initial point
// Given that the point is not necessarily on the border of the cell, we can't use the precomputed dts
std::array<double, 3> initial_ts(const std::pair<int, int>& ij, const Vector2D& phr, const Vector2D& foo, const Spacing& s, LatticeType lattice_type, const DDAState& dda_state) {
    if (lattice_type == LatticeType::Rectangular) {
        return initial_ts_rectangle(ij, phr, foo, std::get<Vector2D>(s), std::get<RectDDAState>(dda_state));
    } else {
        return initial_ts_hexagon(ij, phr, foo, std::get<double>(s), std::get<HexDDAState>(dda_state));
    }
}

std::pair<int, int> dda_step_rectangle(std::array<double, 3>& ts, const RectDDAState& state, double& t_exit) {
    double& tx = ts[0]; double& ty = ts[1];
    if (tx <= ty) {
        t_exit = tx; tx += state.dtx;
        return {state.di, 0};
    } else {
        t_exit = ty; ty += state.dty;
        return {0, state.dj};
    }
}

std::pair<int, int> dda_step_hexagon(std::array<double, 3>& ts, const HexDDAState& state, double& t_exit) {
    double& t1 = ts[0]; double& t2 = ts[1]; double& t3 = ts[2];
    if (t1 <= t2 && t1 <= t3) {
        t_exit = t1; t1 += state.dt1;
        return {state.di1, state.dj1};
    } else if (t2 <= t1 && t2 <= t3) {
        t_exit = t2; t2 += state.dt2;
        return {state.di2, state.dj2};
    } else {
        t_exit = t3; t3 += state.dt3;
        return {state.di3, state.dj3};
    }
}

// Return the cell index update values after a single voxel has been traversed
// Modifies `t_exit`
std::pair<int, int> dda_step(std::array<double, 3>& ts, const DDAState& dda_state, double& t_exit) {
    if (auto* rect = std::get_if<RectDDAState>(&dda_state)) {
        return dda_step_rectangle(ts, *rect, t_exit);
    } else {
        return dda_step_hexagon(ts, std::get<HexDDAState>(dda_state), t_exit);
    }
}
