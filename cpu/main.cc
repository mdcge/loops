#include <iostream>

#include "vector.hh"
#include "photon.hh"

int main() {
    Vector pos(0.0, 0.0, 0.0);
    Vector dir(1.0, 0.0, 0.0);

    Photon photon(pos, dir, 0.0);

    for (int i=0; i<5; i++) {
        std::cout << "Position: (" << photon.r.x << ", " << photon.r.y << ", " << photon.r.z << ") at time " << photon.t << std::endl;
        photon.propagate(i);
    }
    
    return 0;
}
