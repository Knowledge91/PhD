//
// Created by Dirk Hornung on 12/6/17.
//

#ifndef PHD_INTEGRALMOMENTUM_H
#define PHD_INTEGRALMOMENTUM_H

#include <complex>
#include <AdlerFunction.h>
#include <typeinfo>

using namespace std::complex_literals;

class IntegralMomentum {
public:
    IntegralMomentum(Constants constants) : constants(constants), adlerFunction(AdlerFunction(constants)) {};

    std::complex<double> contourIntegral(double s0) {
        std::complex<double> lowerLimit, upperLimit, ii;
        lowerLimit = 0.;
        upperLimit = 2*Constants::pi();
        ii = 1.i;


        auto preparedContourFunction = [&](std::complex<double> x) -> std::complex<double> {
            return adlerFunction.D0(s0*exp(ii * x), Constants::mz());
        };

       return NumericalMethods::integrate(preparedContourFunction, lowerLimit, upperLimit);
    }

private:
    AdlerFunction adlerFunction;
    Constants constants;
};


#endif //PHD_INTEGRALMOMENTUM_H
