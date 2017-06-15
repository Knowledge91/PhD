//
// Created by Dirk Hornung on 30/5/17.
//

#include <complex>
#include <cmath>
#include "gtest/gtest.h"
#include "NumericalMethods.h"
#include "RunAlpha.h"
#include "Constants.h"

using namespace std::complex_literals;


class NumericalMethodsFixture : public ::testing::Test {
public:
    static double function(double x) {
        return 2.*std::pow(x, 3) + 7.*std::pow(x, 2) - 13.;
    }

    static std::complex<double> complexFunction(std::complex<double> x) {
        std::complex<double> a = 1.i;
        return exp(- a * x );
    }

protected:
};



TEST_F(NumericalMethodsFixture, GaussianQuadrature) {
    EXPECT_NEAR(NumericalMethods::integrate(function, -5., 13.), 19152., 1e-9);

    std::complex<double> b, c;
    b = 3. + 4.i;
    c = -12. + 13.i;
    EXPECT_NEAR(NumericalMethods::integrate(complexFunction, b, c).real(), 237379.339821324 , 1 );
    EXPECT_NEAR(NumericalMethods::integrate(complexFunction, b, c).imag(), 373386.344001859 , 1e-1 );
}



