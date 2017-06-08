//
// Created by Dirk Hornung on 30/5/17.
//

#include "NumericalMethods.h"
#include <math.h>
#include <cmath>
#include <iostream>

double NumericalMethods::NewtonRaphson(double xi, double (*func)(double), double (*Dfunc)(double)) {
    double y = 1;
    int i = 0;
    while ( std::abs(y) > 1e-15 ) {
        xi = xi - func(xi)/Dfunc(xi);
        y = func(xi);
        std::cout << "f(x) " << y << " - " << xi << std::endl;
        i++;
    }


    return xi;
}

