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

TEST_F(NumericalMethodsFixture, GaussJordan) {
//                        3  0  2             0.2  0.2 0
//     the inverse of     2  0 -2 should be  -0.2  0.3 1
//                        2  1  1             0.2 -0.3 0
    std::vector<std::vector<double> > a(3, std::vector<double>(3));
    a[0][0] = 3.; a[0][1] = 0.; a[0][2] = 2.;
    a[1][0] = 2.; a[1][1] = 0.; a[1][2] = -2.;
    a[2][0] = 0.; a[2][1] = 1.; a[2][2] = 1.;

    NumericalMethods::gaussj(a);
    EXPECT_NEAR(a[0][0], 0.2, 1e-15);
    EXPECT_NEAR(a[0][1], 0.2, 1e-15);
    EXPECT_NEAR(a[0][2], 0.0, 1e-15);
    EXPECT_NEAR(a[1][0], -0.2, 1e-15);
    EXPECT_NEAR(a[1][1], 0.3, 1e-15);
    EXPECT_NEAR(a[1][2], 1, 1e-15);
    EXPECT_NEAR(a[2][0], 0.2, 1e-15);
    EXPECT_NEAR(a[2][1], -0.3, 1e-15);
    EXPECT_NEAR(a[2][2], 0, 1e-15);
}




