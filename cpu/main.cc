#include <iostream>
#include <chrono>
#include <random>

#include "vector.hh"
#include "photon.hh"
#include "parameters.hh"

int main() {
    // Simulation parameters
    const int nb_events = 1;
    const int nb_photons = 1000;
    const int max_steps = 10000;
    
    Vector origin(0.0, 0.0, 0.0);
    Parameters parameters(1.0, 9.0);
    std::mt19937 rng(42);

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    for (int e=0; e<nb_events; e++) {
        for (int n=0; n<nb_photons; n++) {
            Photon photon(origin, 0.0);
            photon.track(parameters, rng, max_steps);
        }
    }

    // Stop timer
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " ms\n";
    
    return 0;
}
