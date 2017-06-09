#include <iostream>
#include "program/Constants.h"

#include "AdlerFunction.h"

int main() {
    std::cout << "Hello, World!" << std::endl;

    std::cout << "Adler Function :" << AdlerFunction::D0(-2., AdlerFunction::mz) << " " << AdlerFunction::D0(-10., 25) << " " << AdlerFunction::D0(-10., 10) << std::endl;
    return 0;
}