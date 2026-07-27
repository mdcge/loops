#include <iostream>
#include <chrono>

#include "vector.hh"
#include "photon.hh"

int main() {
    Vector pos(0.0, 0.0, 0.0);
    Vector dir(1.0, 0.0, 0.0);

    Photon photon(pos, dir, 0.0);

    for (int i=0; i<5; i++) {
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

        std::cout << "Position: (" << photon.r.x << ", " << photon.r.y << ", " << photon.r.z << ") at time " << photon.t << std::endl;
        photon.propagate(i);
    }

    // Stop timer
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Elapsed time: " << duration.count() << " ms\n";
    
    return 0;
}
