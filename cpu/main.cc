#include <iostream>
#include <chrono>
#include <random>

#include "vector.hh"
#include "photon.hh"
#include "parameters.hh"

int main() {
    Vector pos(0.0, 0.0, 0.0);
    Vector dir(1.0, 0.0, 0.0);

    Photon photon(pos, dir, 0.0);
    Parameters parameters(1.0, 9.0);

    std::mt19937 rng(42);

    for (int i=0; i<5; i++) {
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    double absorption_length = parameters.sample_absorption_length(rng);
        double scattering_length = parameters.sample_scattering_length(rng);
        std::cout << "Position: (" << photon.r.x << ", " << photon.r.y << ", " << photon.r.z << ") at time " << photon.t << std::endl;
        std::cout << "Scat. length: " << scattering_length << std::endl;
        photon.propagate(scattering_length);
    }

    // Stop timer
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " ms\n";
    
    return 0;
}
