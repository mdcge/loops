#include <iostream>

#include "vector.hh"

int main() {
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(3.0, -4.0, 0.0);

    std::cout << "Multiply: (" << (v1*3.0).x << ", " << (3.0*v1).y << ", " << (-v1).z << ")" << std::endl;
    
    return 0;
}
