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

// Fresnel reflectance from medium 1 -> medium 2
double fresnel_reflectance(double n1, double n2, double costheta1, double costheta2) {
    double Rs = (n1*costheta1 - n2*costheta2) / (n1*costheta1 + n2*costheta2); // s-polarised
    double Rp = (n2*costheta1 - n1*costheta2) / (n2*costheta1 + n1*costheta2); // p-polarised
    return (Rs*Rs + Rp*Rp) / 2.0;
}

// Interact the photon with the fibre until it exits the fibre
// Return `true` if photon absorbed (in fibre), `false` if photon still alive and in scintillator
bool Simulation::interact_photon_on_fibre(Photon& photon, Vector2D& fc, double fr_core, double fr_innerclad, double fr_outerclad, double n_scint, double n_outerclad, double n_innerclad, double n_fibre, double n_g_fibre) {
    const double EPSILON = 1e-6; // for nudging photon away from fibre edge

    // Initialise uniform number generator
    // std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // === OUTSIDE FIBRE ===
    
    // Normal at entry point to fibre
    Vector2D n2d_out = (Vector2D(photon.r.x, photon.r.y) - fc).norm();  // 2D normal (always pointing outwards from the fibre)
    Vector n_out(n2d_out.x, n2d_out.y, 0.0);
    // Incident angle and TIR condition
    double costhetas = -dot(photon.p, n_out);
    double eta_s_oc = n_scint / n_outerclad;
    double k = 1 - eta_s_oc*eta_s_oc * (1 - costhetas*costhetas);
    if (k <= 0) { // TIR back into scintillator
        photon.p += 2 * costhetas * n_out;
        photon.p = photon.p.norm();
        photon.r += photon.p * EPSILON; // nudge photon to avoid immediate re-intersection
        return false;  // return `false` because photon is still alive -> back to main simulation
    }
    // Remaining angles
    double eta_s_ic = n_scint / n_innerclad;
    double eta_s_f = n_scint / n_fibre;
    double costhetaoc = std::sqrt(k);
    double costhetaic = std::sqrt(1 - eta_s_ic*eta_s_ic * (1 - costhetas*costhetas));
    double costhetaf = std::sqrt(1 - eta_s_f*eta_s_f * (1 - costhetas*costhetas));
    // Reflectances    
    double R_s_oc = fresnel_reflectance(n_scint, n_outerclad, costhetas, costhetaoc);
    double R_oc_ic = fresnel_reflectance(n_outerclad, n_innerclad, costhetaoc, costhetaic);
    double R_ic_f = fresnel_reflectance(n_innerclad, n_fibre, costhetaic, costhetaf);
    double Reff_s_ic = R_s_oc + (1 - R_s_oc) * (1 - R_s_oc) * (R_oc_ic / (1 - R_s_oc * R_oc_ic));
    double Reff_s_f = Reff_s_ic + (1 - Reff_s_ic) * (1 - Reff_s_ic) * (R_ic_f / (1 - Reff_s_ic * R_ic_f));
    // Non-TIR reflection back into scintillator or refraction into fibre?
    double u = uniform01(rng);
    if (u < Reff_s_f) { // non-TIR reflection back into scintillator
        photon.p += 2 * costhetas * n_out;
        photon.p = photon.p.norm();
        photon.r += photon.p * EPSILON; // nudge photon to avoid immediate re-intersection
        return false;  // return `false` because photon is still alive -> back to main simulation
    }

    
    // === INSIDE FIBRE ===

    // Teleport photon to fibre core edge
    double deltar = fr_outerclad - fr_core;
    photon.r -= deltar * n_out;
    // Update direction to refracted direction
    photon.p += (eta_s_f - 1) * photon.p + (eta_s_f * costhetas - costhetaf) * n_out;
    photon.p = photon.p.norm();
    photon.r += photon.p * EPSILON; // nudge photon to avoid immediate re-intersection
    // Distance to exit based on number of reflections
    double dc = circle_intersection_distance(Vector2D(photon.r.x, photon.r.y), Vector2D(photon.p.x, photon.p.y), fc, fr_core);  // chord length
    double NR = std::floor(std::log(uniform01(rng)) / std::log(Reff_s_f)); // number of reflections before exiting
    double fibre_exit_distance = (NR + 1) * dc;
    // Distance to absorption in fibre
    double fibre_absorption_length = params.fibre_absorption_length.at(photon.wl);
    double fibre_absorption_distance = sample_length(fibre_absorption_length, rng);
    // Absorption in fibre or refraction into scintillator?
    if (fibre_absorption_distance < fibre_exit_distance) { // absorption in fibre
        double distance_along_chord = std::fmod(fibre_absorption_distance, dc);
        photon.r.x += distance_along_chord * photon.p.x;
        photon.r.y += distance_along_chord * photon.p.y;
        photon.r.z += fibre_absorption_distance * photon.p.z;
        photon.t += n_g_fibre * C_INV * fibre_absorption_distance;
        return true; // return `true` because photon is dead -> exit main simulation
    }

    // Refraction out to scintillator
    // CONSIDERING ANGLE
    double thetaf = std::acos(costhetaf);
    Vector2D rho = Vector2D(photon.r.x, photon.r.y) - fc;  // radial vector to first entry point
    double sigma = std::copysign(1.0, rho.x*photon.p.y - rho.y*photon.p.x);  // rotation handedness
    double Phi = sigma * NR * (M_PI - 2*thetaf);
    double cosPhi = std::cos(Phi);
    double sinPhi = std::sin(Phi);
    // Rotated vector components
    double rot_rx = cosPhi * rho.x - sinPhi * rho.y;
    double rot_ry = sinPhi * rho.x + cosPhi * rho.y;
    double rot_px = cosPhi * photon.p.x - sinPhi * photon.p.y;
    double rot_py = sinPhi * photon.p.x + cosPhi * photon.p.y;
    // Update position and direction to exit point
    photon.r.x = fc.x + rot_rx + dc * rot_px;
    photon.r.y = fc.y + rot_ry + dc * rot_py;
    photon.r.z += fibre_exit_distance * photon.p.z;
    photon.p.x = rot_px;
    photon.p.y = rot_py;

    // Normal at exit point of fibre
    Vector2D n2d_in = (Vector2D(photon.r.x, photon.r.y) - fc).norm();  // 2D normal (always pointing outwards from the fibre)
    Vector n_in(n2d_in.x, n2d_in.y, 0.0);
    // Update direction and time
    double eta_f_s = n_fibre / n_scint;
    photon.p += (eta_f_s - 1) * photon.p + (costhetas - eta_f_s * costhetaf) * n_in;
    photon.p = photon.p.norm();
    photon.t += n_g_fibre * C_INV * fibre_exit_distance;
    // Teleport photon to fibre edge
    photon.r += deltar * n_in;
    photon.r += photon.p * EPSILON; // nudge photon to avoid immediate re-intersection
    return false;  // return `false` because photon is still alive -> back to main simulation
}

// Track a single photon, from creation to absorption
void Simulation::track_photon(Photon& photon, int max_steps) {
    // === INITIALISATION ===
    double scint_absorption_distance = sample_length(params.scintillator_absorption_length.at(photon.wl), rng);  // sample scintillator absorption distance
    double scint_refractive_index = params.scintillator_refractive_index.at(photon.wl);  // look up scintillator refractive index at this wavelength
    double scint_group_refractive_index = params.scintillator_group_refractive_index.at(photon.wl);  // look up scintillator group velocity refractive index at this wavelength
    double fibre_core_refractive_index = params.fibre_core_refractive_index.at(photon.wl);  // look up fibre core refractive index at this wavelength
    double fibre_core_group_refractive_index = params.fibre_core_group_refractive_index.at(photon.wl);  // look up fibre core group velocity refractive index at this wavelength
    double fibre_inner_cladding_refractive_index = params.fibre_inner_cladding_refractive_index.at(photon.wl);  // look up fibre inner cladding refractive index at this wavelength
    double fibre_outer_cladding_refractive_index = params.fibre_outer_cladding_refractive_index.at(photon.wl);  // look up fibre outer cladding refractive index at this wavelength
    double mean_scattering_length = params.scintillator_scattering_length.at(photon.wl);  // look up mean scintillator scattering length at this wavelength

    // Setup variables to avoid re-computing cell index when exiting fibre
    bool photon_fibre_to_scintillator = false; 
    int last_ci = 0; int last_cj = 0;
    
    // === LOOP ===
    for (int i=0; i<max_steps; i++) {
        
        double scattering_distance = sample_length(mean_scattering_length, rng);  // sample scattering length for this step
        double intersection_distance = detector.intersection_distance(photon.r, photon.p);  // calculate minimum distance to detector walls
        double max_interaction_distance = std::min({scattering_distance, scint_absorption_distance, intersection_distance});  // maximum distance needed to be checked for fibre intersections

        // === VOXEL TRAVERSAL ===
        // Photon properties for voxel traversal
        Vector2D X0(photon.r.x, photon.r.y);
        Vector2D photon_p_2d(photon.p.x, photon.p.y);
        // Calculate initial cell index (i0, j0)
        // If exiting fibre, just reuse fibre entry cell index
        int ci, cj;
        if (photon_fibre_to_scintillator) {
            ci = last_ci;
            cj = last_cj;
        } else {
            std::pair<int, int> ij = cell_index(Vector2D(photon.r.x, photon.r.y), fibres.s, fibres.foo, fibres.lattice_type);
            ci = ij.first; cj = ij.second;
        }
        // Voxel traversal state
        DDAState dda_state = make_dda_state(photon_p_2d, fibres.s, fibres.lattice_type);
        std::array<double, 3> ts = initial_ts(std::make_pair(ci, cj), X0, fibres.foo, fibres.s, fibres.lattice_type, dda_state);

        double t_entry = 0.0;  // distance from photon position at start of ray cast to entry point of current cell
        double voxel_traversal_distance = 0.0;  // how far have we travelled in the grid while searching for fibre intersection
        bool fibre_hit = false;  // have we hit a fibre?
        // March through grid of voxels until `max_interaction_distance` is reached, or fibre intersected
        Vector2D fibre_centre;
        bool first_cell_in_step = true; // is this the first cell in this step?

        // Fibre intersection loop
        while (t_entry < max_interaction_distance) {
            // Does this cell contain a fibre?
            bool cell_has_fibre = fibres.fibre_cells.count(std::make_pair(ci, cj)) > 0;
            
            // Skip this fibre check if
            //     - the photon has been reflected/refracted back into scintillator from the fibre -> `photon_fibre_to_scintillator` == true
            //     AND
            //     - this is the first cell in voxel traversal -> `first_cell_in_step` == true
            bool skip_this_fibre = photon_fibre_to_scintillator && first_cell_in_step;
            
            // Only continue with checks if cell has a fibre AND we're not still same cell (after fibre exit)
            if (cell_has_fibre && !skip_this_fibre) {
                fibre_centre = fibres.fij(fibres.foo, ci, cj);  // fibre position in this cell
                Vector2D r_entry(X0.x + t_entry * photon_p_2d.x, X0.y + t_entry * photon_p_2d.y);  // position of photon at entrance of current cell
                bool already_inside = (r_entry - fibre_centre).mag() <= fibres.fr_outerclad;  // is photon somehow still on inside of fibre?

                if (!already_inside) {
                    double fibre_distance = circle_intersection_distance(r_entry, photon_p_2d, fibre_centre, fibres.fr_outerclad);

                    if ((t_entry + fibre_distance) < max_interaction_distance) { // if fibre intersection happens before other processes, set `fibre_hit` flag to true and exit loop
                        fibre_hit = true;
                        voxel_traversal_distance = t_entry + fibre_distance;
                        break;
                    }
                }
            }

            // If no fibre intersection in this cell or no fibre in this cell, increment total distance and step to next cell, updating cell indices
            double t_exit;
            std::pair<int, int> cell_index_steps = dda_step(ts, dda_state, t_exit);
            ci += cell_index_steps.first;
            cj += cell_index_steps.second;
            t_entry = t_exit;
            first_cell_in_step = false;
        }
        // Reset to `false` by default: only photons exiting fibres can set to `true`
        photon_fibre_to_scintillator = false;

        // === INTERACTION TYPE ===
        if (fibre_hit) { // if fibre is hit, enter fibre interaction simulation
            propagate_photon(photon, voxel_traversal_distance, scint_group_refractive_index);
            bool absorbed = interact_photon_on_fibre(photon, fibre_centre, fibres.fr_core, fibres.fr_innerclad, fibres.fr_outerclad, scint_refractive_index, fibre_outer_cladding_refractive_index, fibre_inner_cladding_refractive_index, fibre_core_refractive_index, fibre_core_group_refractive_index);
            if (absorbed) { // if photon was absorbed in fibre, break out of loop
                break;
            }
            // Record fibre entry/exit cell and set flag to `true`
            photon_fibre_to_scintillator = true;
            last_ci = ci;
            last_cj = cj;
        } else if (scint_absorption_distance < scattering_distance && scint_absorption_distance < intersection_distance) {  // if absorption happens first, break out of the loop
            propagate_photon(photon, scint_absorption_distance, scint_group_refractive_index);
            break;
        } else if (scattering_distance < scint_absorption_distance && scattering_distance < intersection_distance) {  // if scattering happens first, propagate the photon and resample direction
            scint_absorption_distance -= scattering_distance;
            propagate_photon(photon, scattering_distance, scint_group_refractive_index);
            photon.p = sample_direction(rng);
        } else {  // if intersection happens first, break out of loop
            propagate_photon(photon, intersection_distance, scint_group_refractive_index);
            break;
        }
    }
}
