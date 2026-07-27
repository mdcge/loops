#include <iostream>

#include "vector.hh"

int main() {
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(3.0, -4.0, 0.0);

    std::cout << "Magnitude: " << v2.mag() << std::endl;

    std::cout << "Normalised: (" << v2.norm().x << ", " << v2.norm().y << ", " << v2.norm().z << ")" << std::endl;
    
    return 0;
}
