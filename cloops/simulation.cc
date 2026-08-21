#include "simulation.hh"

#include "geometry.hh"

#include <cmath>
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

// Interact the photon with the fibre until it exits the fibre
// Return `true` if photon absorbed (in fibre), `false` if photon still alive and in scintillator
bool Simulation::interact_photon_on_fibre(Photon& photon, Vector2D& fc, double fr, double n_scint, double n_fibre) {
    const double EPSILON = 0.0001; // for nudging photon away from fibre edge

    // Track whether photon is inside or outside fibre
    bool in_fibre = false;

    // Initialise uniform number generator
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // Sample fibre absorption length at this photon wavelength
    double fibre_absorption_length = params.fibre_absorption_length.at(photon.wl);

    while (true) {
        // Normal vector
        Vector2D n0 = (Vector2D(photon.r.x, photon.r.y) - fc).norm();
        Vector2D n2d = -std::copysign(1.0, dot(Vector2D(photon.p.x, photon.p.y), n0)) * n0;
        Vector n(n2d.x, n2d.y, 0.0);

        // Reflection/refraction coefficients
        double n1 = in_fibre ? n_fibre : n_scint;
        double n2 = in_fibre ? n_scint : n_fibre;
        double eta = n1 / n2;
        double costheta1 = -dot(photon.p, n);
        double k = 1 - eta*eta * (1 - costheta1*costheta1);
        double R0 = std::pow((n_scint - n_fibre) / (n_scint + n_fibre), 2);  // reflectance at normal incidence
        double R = R0 + (1 - R0) * std::pow(1 - costheta1, 5);  // reflection probability (Schlick approximation)

        double u = uniform(rng);
        if ((k < 0) || (u < R)) { // reflection (TIR or Fresnel)
            photon.p += 2 * costheta1 * n;
            photon.p = photon.p.norm();  // guard against floating point drift
            photon.r += photon.p * EPSILON;  // nudge photon position away from fibre edge to avoid immediate re-interaction
        } else { // refraction
            photon.p += (eta - 1) * photon.p + (eta * costheta1 - std::sqrt(k)) * n;
            photon.p = photon.p.norm();  // guard against floating point drift
            photon.r += photon.p * EPSILON;  // nudge photon position away from fibre edge to avoid immediate re-interaction
            in_fibre = !in_fibre;  // refraction changes in <-> out of fibre state
        }

        // Exit to main loop if photon outside of fibre
        // Return `false` because photon is still alive
        if (!in_fibre) {
            return false;
        }

        // Propagate inside fibre
        double fibre_absorption_distance = sample_length(fibre_absorption_length, rng);
        double intersection_distance = circle_intersection_distance(Vector2D(photon.r.x, photon.r.y), Vector2D(photon.p.x, photon.p.y), fc, fr);

        if (fibre_absorption_distance < intersection_distance) { // if photon absorbed in fibre, exit this loop and outer loop
            propagate_photon(photon, fibre_absorption_distance, n_fibre);
            return true; // return `true` because photon absorbed
        }

        // If photon reaches fibre edge, we're still in the fibre so allow loop to restart
        propagate_photon(photon, intersection_distance, n_fibre);
    }
}

// Track a single photon, from creation to absorption
void Simulation::track_photon(Photon& photon, int max_steps) {
    // === INITIALISATION ===
    double scint_absorption_distance = sample_length(params.scintillator_absorption_length.at(photon.wl), rng);  // sample scintillator absorption distance
    double scint_refractive_index = params.scintillator_refractive_index.at(photon.wl);  // look up scintillator refractive index at this wavelength
    double fibre_refractive_index = params.fibre_refractive_index.at(photon.wl);  // look up fibre refractive index at this wavelength
    double mean_scattering_length = params.scintillator_scattering_length.at(photon.wl);  // look up mean scintillator scattering length at this wavelength
    
    // === LOOP ===
    for (int i=0; i<max_steps; i++) {
        
        double scattering_distance = sample_length(mean_scattering_length, rng);  // sample scattering length for this step
        double intersection_distance = detector.intersection_distance(photon.r, photon.p);  // calculate minimum distance to detector walls
        double max_interaction_distance = std::min({scattering_distance, scint_absorption_distance, intersection_distance});  // maximum distance needed to be checked for fibre intersections

        // === VOXEL TRAVERSAL ===
        // Photon properties for voxel traversal
        Vector2D photon_r_2d(photon.r.x, photon.r.y);
        Vector2D photon_p_2d(photon.p.x, photon.p.y);
        // Calculate initial cell index (i0, j0)
        std::pair<int, int> ij = cell_index(Vector2D(photon.r.x, photon.r.y), fibres.s, fibres.foo, fibres.lattice_type);
        int ci = ij.first; int cj = ij.second;
        // Voxel traversal state variables
        double voxel_traversal_distance = 0.0;  // how far have we travelled in the grid while searching for fibre intersection
        bool fibre_hit = false;  // have we hit a fibre?

        // March through grid of voxels until `max_interaction_distance` is reached, or fibre intersected
        Vector2D fibre_centre;
        while (voxel_traversal_distance < max_interaction_distance) {
            // Distance to exit current cell `t`
            std::array<double, 3> ts = cell_ts(std::make_pair(ci, cj), photon_r_2d, photon_p_2d, fibres.foo, fibres.s, fibres.lattice_type);
            double cell_exit_distance = *std::min_element(ts.begin(), ts.end());

            // Distance to fibre in this cell (infinity if no intersection)
            fibre_centre = fibres.fij(fibres.foo, ci, cj);
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
            std::pair<int, int> new_ij = new_cell_index(std::make_pair(ci, cj), ts, photon_p_2d, fibres.lattice_type);
            ci = new_ij.first;
            cj = new_ij.second;
            photon_r_2d = new_photon_position(ts, photon_r_2d, photon_p_2d, fibres.lattice_type);
        }

        // === INTERACTION TYPE ===
        if (fibre_hit) { // if fibre is hit, enter fibre interaction simulation
            propagate_photon(photon, voxel_traversal_distance, scint_refractive_index);
            bool absorbed = interact_photon_on_fibre(photon, fibre_centre, fibres.fr, scint_refractive_index, fibre_refractive_index);
            if (absorbed) { // if photon was absorbed in fibre, break out of loop
                break;
            }
        } else if (scint_absorption_distance < scattering_distance && scint_absorption_distance < intersection_distance) {  // if absorption happens first, break out of the loop
            propagate_photon(photon, scint_absorption_distance, scint_refractive_index);
            break;
        } else if (scattering_distance < scint_absorption_distance && scattering_distance < intersection_distance) {  // if scattering happens first, propagate the photon and resample direction
            scint_absorption_distance -= scattering_distance;
            propagate_photon(photon, scattering_distance, scint_refractive_index);
            photon.p = sample_direction(rng);
        } else {  // if intersection happens first, break out of loop
            propagate_photon(photon, intersection_distance, scint_refractive_index);
            break;
        }
    }
}
