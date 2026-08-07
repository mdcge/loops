#ifndef FIBRES_H_
#define FIBRES_H_

#include <vector>
#include <limits>
#include <variant>
#include <set>

#include "vector.hh"

enum struct LatticeType { Rectangular, Hexagonal };
using Spacing = std::variant<Vector2D, double>;

struct Fibres {
    std::vector<Vector2D> fs;  // fibre positions [mm]
    double fr;  // fibre radius [mm]
    LatticeType lattice_type;  // type of lattice: rectangular or hexagonal
    Vector2D foo;  // 2D position of "central" fibre (closest to detector origin)
    Spacing s;  // fibre spacing
    std::set<std::pair<int, int>> fibre_cells;  // set of cells which contain a fibre

    Fibres(const std::vector<Vector2D>& fibre_positions, double fibre_radius, LatticeType lattice_type)
        : fs(fibre_positions),
          fr(fibre_radius),
          lattice_type(lattice_type),
          s(lattice_type == LatticeType::Rectangular
                ? Spacing(spacing_rectangle())
                : Spacing(spacing_hexagon()))
    {}

    Vector2D f00(const Vector2D&) const;
    Vector2D fij(const Vector2D&, int, int);
    void build_fibre_cells(const Vector2D&);
    
    // Grid-specific functions
    Vector2D fij_rectangle(const Vector2D&, int, int, const Vector2D&) const;
    Vector2D fij_hexagon(const Vector2D&, int, int, double) const;
    Vector2D spacing_rectangle() const;
    double spacing_hexagon() const;
    void build_fibre_cells_rectangle(const Vector2D&, const Vector2D&, const Vector2D&);
    void build_fibre_cells_hexagon(const Vector2D&, const Vector2D&, double);
};

std::pair<int, int> cell_index(const Vector2D&, const Vector2D&, const Spacing&, const Vector2D&, LatticeType);

// Grid-specific functions
std::pair<int, int> cell_index_rectangle(const Vector2D&, const Vector2D&, const Vector2D&, const Vector2D&);
std::pair<int, int> cell_index_hexagon(const Vector2D&, const Vector2D&, double, const Vector2D&);
std::pair<double, double> cell_ts_rectangle(const std::pair<int, int>, const Vector2D&, const Vector2D&, const Vector2D&, const Vector2D&);
std::pair<int, int> new_cell_index_rectangle(const std::pair<int, int>, const std::pair<double, double>, const Vector2D&);
Vector2D new_photon_position_rectangle(const std::pair<double, double>, const Vector2D&, const Vector2D&);

#endif // FIBRES_H_
