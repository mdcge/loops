#include <iostream>
#include <chrono>
#include <random>
#include <fstream>

#include "vector.hh"
#include "photon.hh"
#include "parameters.hh"

int main() {
    std::ofstream file("filename.txt");
    
    // Simulation parameters
    const int nb_events = 1;
    const int nb_photons = 1000;
    const int max_steps = 100000;
    
    Vector origin(0.0, 0.0, 0.0);
    Parameters parameters(1.649, 22100.0);
    std::mt19937 rng(42);

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    for (int e=0; e<nb_events; e++) {
        for (int n=0; n<nb_photons; n++) {
            Photon photon(origin, 0.0);
            photon.track(parameters, 1.496, rng, max_steps);

            file << photon.r.x << "," << photon.r.y << "," << photon.r.z << "," << photon.t << "\n";
        }
    }

    file.close();

    // Stop timer
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " ms\n";
    
    return 0;
}
