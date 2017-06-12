#include <iostream>
#include "Constants.h"
#include "IntegralMomentum.h"

int main() {
    std::cout << "Hello, World!" << std::endl;
    Constants constants(3,3,4);
    IntegralMomentum integralMomentum(constants);

    std::cout << " countour " << integralMomentum.contourIntegral(1.) << std::endl;
    return 0;
}