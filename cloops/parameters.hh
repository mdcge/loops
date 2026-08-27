#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <random>
#include <vector>

#include "vector.hh"
#include "xoshiro.hh"

// Pair of wavelength and value arrays
struct OpticalProperty {
    std::vector<double> wavelengths;  // wavelengths at which values are evaluated [nm]
    std::vector<double> values;  // values at reference wavelengths

    OpticalProperty(std::vector<double> wavelengths, std::vector<double> values) : wavelengths(wavelengths), values(values) {}

    // Return value of `values` at `wavelength`
    double at(double) const;
};

struct Parameters {
    OpticalProperty emission_spectrum;  // emission spectrum of the scintillator [nm] -> [1]
    OpticalProperty scintillator_scattering_length;  // scattering length of the scintillator [nm] -> [mm]
    OpticalProperty scintillator_absorption_length;  // absorption length of the scintillator [nm] -> [mm]
    OpticalProperty fibre_absorption_length;  // absorption length of the fibres [nm] -> [mm]
    OpticalProperty scintillator_refractive_index;  // refractive index of the scintillator [nm] -> [1]
    OpticalProperty fibre_core_refractive_index;  // refractive index of the fibre core [nm] -> [1]
    OpticalProperty fibre_inner_cladding_refractive_index;  // refractive index of the fibre inner cladding [nm] -> [1]
    OpticalProperty fibre_outer_cladding_refractive_index;  // refractive index of the fibre inner cladding [nm] -> [1]
    double scint_rise_time;  // scintillation rise time [ns]
    std::vector<double> scint_time_cnsts;  // scintillation time constants (used in exponential decay) [ns]
    std::vector<double> scint_time_amps;  // relative amplitude of each scintillation time constant

    Parameters(
        OpticalProperty emission_spectrum,
        OpticalProperty scintillator_scattering_length,
        OpticalProperty scintillator_absorption_length,
        OpticalProperty fibre_absorption_length,
        OpticalProperty scintillator_refractive_index,
        OpticalProperty fibre_core_refractive_index,
        OpticalProperty fibre_inner_cladding_refractive_index,
        OpticalProperty fibre_outer_cladding_refractive_index,
        double scintillation_rise_time,
        std::vector<double> scintillation_time_constants,
        std::vector<double> scintillation_time_amplitudes
    ) : emission_spectrum(emission_spectrum),
        scintillator_scattering_length(scintillator_scattering_length),
        scintillator_absorption_length(scintillator_absorption_length),
        fibre_absorption_length(fibre_absorption_length),
        scintillator_refractive_index(scintillator_refractive_index),
        fibre_core_refractive_index(fibre_core_refractive_index),
        fibre_inner_cladding_refractive_index(fibre_inner_cladding_refractive_index),
        fibre_outer_cladding_refractive_index(fibre_outer_cladding_refractive_index),
        scint_rise_time(scintillation_rise_time),
        scint_time_cnsts(scintillation_time_constants),
        scint_time_amps(scintillation_time_amplitudes) {}
};

double sample_length(double, xoshiro256pp&);
Vector sample_direction(xoshiro256pp&);

#endif // PARAMETERS_H_
