#ifndef FIBRES_H_
#define FIBRES_H_

#include <vector>
#include <array>
#include <limits>
#include <variant>
#include <set>

#include "vector.hh"

enum struct LatticeType { Rectangular, Hexagonal };
using Spacing = std::variant<Vector2D, double>;

LatticeType infer_lattice_type(const std::vector<Vector2D>&);

struct Fibres {
    std::vector<Vector2D> fs;  // fibre positions [mm]
    double fr_core;  // fibre core radius [mm]
    double fr_innerclad;  // fibre inner cladding radius [mm]
    double fr_outerclad;  // fibre outer cladding radius [mm]
    LatticeType lattice_type;  // type of lattice: rectangular or hexagonal
    Vector2D foo;  // 2D position of "central" fibre (closest to detector origin)
    Spacing s;  // fibre spacing
    std::set<std::pair<int, int>> fibre_cells;  // set of cells which contain a fibre

    Fibres(const std::vector<Vector2D>& fibre_positions, double fibre_core_radius, double fibre_innerclad_radius, double fibre_outerclad_radius)
        : fs(fibre_positions),
          fr_core(fibre_core_radius),
          fr_innerclad(fibre_innerclad_radius),
          fr_outerclad(fibre_outerclad_radius),
          lattice_type(infer_lattice_type(fibre_positions)),
          s(lattice_type == LatticeType::Rectangular
                ? Spacing(spacing_rectangle())
                : Spacing(spacing_hexagon()))
    {}

    Vector2D f00(const Vector2D&) const;
    Vector2D fij(const Vector2D&, int, int);
    void build_fibre_cells();
    
    // Grid-specific functions
    Vector2D fij_rectangle(const Vector2D&, int, int, const Vector2D&) const;
    Vector2D fij_hexagon(const Vector2D&, int, int, double) const;
    Vector2D spacing_rectangle() const;
    double spacing_hexagon() const;
    void build_fibre_cells_rectangle(const Vector2D&, const Vector2D&);
    void build_fibre_cells_hexagon(const Vector2D&, double);
};

std::pair<int, int> cell_index(const Vector2D&, const Spacing&, const Vector2D&, LatticeType);
std::pair<int, int> new_cell_index(const std::pair<int, int>, const std::array<double, 3>&, const Vector2D&, LatticeType);
std::array<double, 3> cell_ts(const std::pair<int, int>, const Vector2D&, const Vector2D&, const Vector2D&, const Spacing&, LatticeType);
Vector2D new_photon_position(const std::array<double, 3>&, const Vector2D&, const Vector2D&, LatticeType);

// Grid-specific functions
std::pair<int, int> cell_index_rectangle(const Vector2D&, const Vector2D&, const Vector2D&);
std::pair<int, int> cell_index_hexagon(const Vector2D&, double, const Vector2D&);
std::pair<int, int> new_cell_index_rectangle(const std::pair<int, int>, const std::array<double, 3>&, const Vector2D&);
std::pair<int, int> new_cell_index_hexagon(const std::pair<int, int>, const std::array<double, 3>&, const Vector2D&);
std::array<double, 3> cell_ts_rectangle(const std::pair<int, int>, const Vector2D&, const Vector2D&, const Vector2D&, const Vector2D&);
std::array<double, 3> cell_ts_hexagon(const std::pair<int, int>, const Vector2D&, const Vector2D&, const Vector2D&, double);
Vector2D new_photon_position_rectangle(const std::array<double, 3>&, const Vector2D&, const Vector2D&);
Vector2D new_photon_position_hexagon(const std::array<double, 3>&, const Vector2D&, const Vector2D&);


// DDA state objects
struct HexDDAState {
    double inv_dot1, inv_dot2, inv_dot3;  // 1 / (n_k . pip), compute once instead of for every cell
    double dt1, dt2, dt3;  // t increment for crossing a full cell in each direction
    double di1, dj1, di2, dj2, di3, dj3;  // index change when crossing cell boundary in each direction
    int bsign1, bsign2, bsign3;  // sign of n_k . php, i.e. direction of boundary crossing
};

struct RectDDAState {
    double inv_dotx, inv_doty;  // 1 / (n_k . pip), compute once instead of for every cell
    double dtx, dty;  // t increment for crossing a full cell in each direction
    double di, dj;  // index change when crossing cell boundary in each direction
};

using DDAState = std::variant<RectDDAState, HexDDAState>;

DDAState make_dda_state(const Vector2D&, const Spacing&, LatticeType);
std::array<double, 3> initial_ts(const std::pair<int, int>&, const Vector2D&, const Vector2D&, const Spacing&, LatticeType, const DDAState&);
std::pair<int, int> dda_step(std::array<double, 3>&, const DDAState&, double&);

RectDDAState make_dda_state_rectangle(const Vector2D&, const Vector2D& s);
HexDDAState make_dda_state_hexagon(const Vector2D&, double s);
std::array<double, 3> initial_ts_rectangle(const std::pair<int, int>&, const Vector2D&, const Vector2D&, const Vector2D&, const RectDDAState&);
std::array<double, 3> initial_ts_hexagon(const std::pair<int, int>&, const Vector2D&, const Vector2D&, double, const HexDDAState&);
std::pair<int, int> dda_step_rectangle(std::array<double, 3>&, const RectDDAState&, double&);
std::pair<int, int> dda_step_hexagon(std::array<double, 3>&, const HexDDAState&, double&);

#endif // FIBRES_H_
