#include <iostream>

#include "vector.hh"

int main() {
    Vector v1(1.0, 2.0, 3.0);
    Vector v2(2.0, -3.0, 0.0);

    std::cout << "Dot product: " << dot(v1, v2) << std::endl;
    
    return 0;
}
