#include <iostream>
#include "program/Constants.h"

#include "RunAlpha.h"

int main() {
    std::cout << "Hello, World!" << std::endl;
    RunAlpha runAlpha(10.);
    std::cout << runAlpha.runAlpha() << std::endl;
    return 0;
}