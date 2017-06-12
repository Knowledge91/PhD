//
// Created by Dirk Hornung on 12/6/17.
//

#ifndef PHD_INTEGRALMOMENTUM_H
#define PHD_INTEGRALMOMENTUM_H

#import <complex>
#import <AdlerFunction.h>
#include <typeinfo>

using namespace std::complex_literals;

class IntegralMomentum {
public:
    IntegralMomentum(Constants constants) : constants(constants), adlerFunction(AdlerFunction(constants)) {};

    std::complex<double> contourIntegral(double s0) {
        std::complex<double> scale, lowerLimit, upperLimit;
        scale = 1.;
        lowerLimit = 0.;
        upperLimit = Constants::pi();


        auto preparedContourFunction = [&]() -> std::complex<double> {
            std::complex<double> x;
            return adlerFunction.D0(s0*exp(x), Constants::mz());
        }();
       return preparedContourFunction;
//        return NumericalMethods::integrate(adlerFunction.D0, scale, lowerLimit, upperLimit);
    }

private:
    AdlerFunction adlerFunction;
    Constants constants;
};


#endif //PHD_INTEGRALMOMENTUM_H
