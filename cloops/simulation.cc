#include "simulation.hh"

#include "geometry.hh"

#include <algorithm>

const double C = 299.8;  // mm/ns
const double C_INV = 1.0 / C;  // ns/mm

Photon Simulation::generate_photon(Vector origin) {
    Photon photon;
    photon.r = origin;
    photon.p = sample_direction(rng);
    photon.wl = params.emission_spectrum.wavelengths[emission_spectrum_distribution(rng)];
    
    int time_component = scintillation_component_distribution(rng);
    std::exponential_distribution<double> rise(1.0 / params.scint_rise_time);
    std::exponential_distribution<double> decay(-1.0 / params.scint_time_cnsts[time_component]);
    photon.t = rise(rng) + decay(rng);
    
    return photon;
}

// Propagate photon by distance `d` [mm] in refractive index `n`
void Simulation::propagate_photon(Photon& photon, double d, double n) {
    photon.r += d * photon.p;
    photon.t += d * n * C_INV;
}

// Track a single photon, from creation to absorption
void Simulation::track_photon(Photon& photon, int max_steps) {
    // === INITIALISATION ===
    double absorption_distance = sample_length(params.absorption_length.at(photon.wl), rng);  // sample absorption distance
    double refractive_index = params.refractive_index.at(photon.wl);  // look up refractive index at this wavelength
    double mean_scattering_length = params.scattering_length.at(photon.wl);  // look up mean scattering length at this wavelength
    Vector2D s = fibres.spacing_rectangle();  // calculate grid spacings
    Vector2D foo = fibres.f00(Vector2D(detector.O().x, detector.O().y));  // calculate (0, 0) fibre position
    
    // === LOOP ===
    for (int i=0; i<max_steps; i++) {
        
        double scattering_distance = sample_length(mean_scattering_length, rng);  // sample scattering length for this step
        double intersection_distance = detector.intersection_distance(photon.r, photon.p);  // calculate minimum distance to detector walls
        double max_interaction_distance = std::min({scattering_distance, absorption_distance, intersection_distance});  // maximum distance needed to be checked for fibre intersections

        // === VOXEL TRAVERSAL ===
        // Photon properties for voxel traversal
        Vector2D photon_r_2d(photon.r.x, photon.r.y);
        Vector2D photon_p_2d(photon.p.x, photon.p.y);
        // Calculate initial cell index (i0, j0)
        std::pair<int, int> ij = cell_index(Vector2D(photon.r.x, photon.r.y), Vector2D(detector.O().x, detector.O().y), fibres.s, fibres.foo, fibres.lattice_type);
        int ci = ij.first; int cj = ij.second;
        // Voxel traversal state variables
        double voxel_traversal_distance = 0.0;  // how far have we travelled in the grid while searching for fibre intersection
        bool fibre_hit = false;  // have we hit a fibre?

        // March through grid of voxels until `max_interaction_distance` is reached, or fibre intersected
        while (voxel_traversal_distance < max_interaction_distance) {
            // Distance to exit current cell `t`
            std::pair<double, double> ts = cell_ts_rectangle(std::make_pair(ci, cj), photon_r_2d, photon_p_2d, foo, s);
            double cell_exit_distance = std::min(ts.first, ts.second);

            // Distance to fibre in this cell (infinity if no intersection)
            Vector2D fibre_centre = fibres.fij(foo, ci, cj);
            bool cell_has_fibre = fibres.fibre_cells.count(std::make_pair(ci, cj)) > 0;  // is this cell in the `fibre_cells` set: does this cell contain a fibre?
            
            // Only continue with checks if cell has a fibre
            if (cell_has_fibre) {
                double fibre_distance = circle_intersection_distance(photon_r_2d, photon_p_2d, fibre_centre, fibres.fr);

                if ((voxel_traversal_distance + fibre_distance) < max_interaction_distance) { // if fibre intersection happens before other processes, set `fibre_hit` flag to true and exit loop
                    fibre_hit = true;
                    voxel_traversal_distance += fibre_distance;
                    break;
                }
            }

            // If no fibre intersection in this cell or no fibre in this cell, increment total distance and step to next cell, updating cell indices and photon position
            voxel_traversal_distance += cell_exit_distance;
            std::pair<int, int> new_ij = new_cell_index_rectangle(std::make_pair(ci, cj), ts, photon_p_2d);
            ci = new_ij.first;
            cj = new_ij.second;
            photon_r_2d = new_photon_position_rectangle(ts, photon_r_2d, photon_p_2d);
        }

        // === INTERACTION TYPE ===
        if (fibre_hit) { // if fibre is hit, break out of loop 
            propagate_photon(photon, voxel_traversal_distance, refractive_index);
            break;
        } else if (absorption_distance < scattering_distance && absorption_distance < intersection_distance) {  // if absorption happens first, break out of the loop
            propagate_photon(photon, absorption_distance, refractive_index);
            break;
        } else if (scattering_distance < absorption_distance && scattering_distance < intersection_distance) {  // if scattering happens first, propagate the photon and resample direction
            absorption_distance -= scattering_distance;
            propagate_photon(photon, scattering_distance, refractive_index);
            photon.p = sample_direction(rng);
        } else {  // if intersection happens first, break out of loop
            propagate_photon(photon, intersection_distance, refractive_index);
            break;
        }
    }
}
