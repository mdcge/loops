#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <random>
#include <vector>

#include "vector.hh"

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
    OpticalProperty scattering_length;  // scattering length of the medium [nm] -> [mm]
    OpticalProperty absorption_length;  // absorption length of the medium [nm] -> [mm]
    OpticalProperty refractive_index;  // refractive index of the medium [nm] -> [1]
    double scint_rise_time;  // scintillation rise time [ns]
    std::vector<double> scint_time_cnsts;  // scintillation time constants (used in exponential decay) [ns]
    std::vector<double> scint_time_amps;  // relative amplitude of each scintillation time constant

    Parameters(
        OpticalProperty emission_spectrum,
        OpticalProperty scattering_length,
        OpticalProperty absorption_length,
        OpticalProperty refractive_index,
        double scintillation_rise_time,
        std::vector<double> scintillation_time_constants,
        std::vector<double> scintillation_time_amplitudes
    ) : emission_spectrum(emission_spectrum),
        scattering_length(scattering_length),
        absorption_length(absorption_length),
        refractive_index(refractive_index),
        scint_rise_time(scintillation_rise_time),
        scint_time_cnsts(scintillation_time_constants),
        scint_time_amps(scintillation_time_amplitudes) {}

    double sample_scattering_length(double, std::mt19937&);
    double sample_absorption_length(double, std::mt19937&);
    double lookup_refractive_index(double);
};

double sample_length(double, std::mt19937&);
Vector sample_direction(std::mt19937&);

#endif // PARAMETERS_H_
