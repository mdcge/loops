#include <iostream>
#include <chrono>
#include <random>

#include "vector.hh"
#include "photon.hh"
#include "parameters.hh"

int main() {
    Vector origin(0.0, 0.0, 0.0);
    
    Parameters parameters(1.0, 9.0);

    std::mt19937 rng(42);

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    for (int n=0; n<10; n++) {
        Photon photon(origin, 0.0);
        photon.track(parameters, rng, 1000);
        std::cout << "Photon final position: (" << photon.r.x << ", " << photon.r.y << ", " << photon.r.z << ")" << std::endl;
    }

    // Stop timer
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " ms\n";
    
    return 0;
}
