#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <random>
#include <vector>

#include "vector.hh"

struct Parameters {
    std::vector<double> emission_wls;  // wavelengths from which photon wavelength is sampled [nm]
    std::vector<double> emission_prob;  // relative probabilities of photon wavelength (i.e. emission spectrum)
    std::vector<double> ls_wls;  // wavelengths at which scattering length is evaluated [nm]
    std::vector<double> ls;  // corresponding absorption lengths [mm]
    std::vector<double> la_wls;  // wavelengths at which absorption length is evaluated [nm]
    std::vector<double> la;  // corresponding scattering lengths [mm]
    std::vector<double> n_wls;  // wavelengths at which refractive index is evaluated [nm]
    std::vector<double> n;  // corresponding refractive indices
    double scint_rise_time;  // scintillation rise time [ns]
    std::vector<double> scint_time_cnsts;  // scintillation time constants (used in exponential decay) [ns]
    std::vector<double> scint_time_amps;  // relative amplitude of each scintillation time constant

    Parameters(
        std::vector<double> wavelengths_emission, std::vector<double> probability_emission,
        std::vector<double> wavelengths_scattering, std::vector<double> scattering_lengths,
        std::vector<double> wavelengths_absorption, std::vector<double> absorption_lengths,
        std::vector<double> wavelengths_refractive, std::vector<double> refractive_indices,
        double scintillation_rise_time, std::vector<double> scintillation_time_constants, std::vector<double> scintillation_time_amplitudes
    ) : emission_wls(wavelengths_emission),
        emission_prob(probability_emission),
        ls_wls(wavelengths_scattering),
        ls(scattering_lengths),
        la_wls(wavelengths_absorption),
        la(absorption_lengths),
        n_wls(wavelengths_refractive),
        n(refractive_indices),
        scint_rise_time(scintillation_rise_time),
        scint_time_cnsts(scintillation_time_constants),
        scint_time_amps(scintillation_time_amplitudes) {}

    double sample_scattering_length(double, std::mt19937&);
    double sample_absorption_length(double, std::mt19937&);
    double lookup_refractive_index(double);
};

double lookup_value(const std::vector<double>&, const std::vector<double>&, double);
double sample_length(double, std::mt19937&);
Vector sample_direction(std::mt19937&);

#endif // PARAMETERS_H_
