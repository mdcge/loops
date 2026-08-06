#ifndef FIBRES_H_
#define FIBRES_H_

#include <vector>
#include <limits>
#include <variant>

#include "vector.hh"

struct Fibres {
    std::vector<Vector2D> fs;  // fibre positions [mm]
    double fr;  // fibre radius [mm]

    Fibres(std::vector<Vector2D> fibre_positions, double fibre_radius) : fs(fibre_positions), fr(fibre_radius) {}

    Vector2D f00(Vector2D) const;
    Vector2D spacing_rectangle();
    double spacing_hexagon();
};

std::pair<int, int> cell_index_rectangle(Vector2D, Vector2D, Vector2D, Vector2D);
std::pair<int, int> cell_index_hexagon(Vector2D, Vector2D, double, Vector2D);

#endif // FIBRES_H_
