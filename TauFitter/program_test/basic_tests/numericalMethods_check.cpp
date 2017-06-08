//
// Created by Dirk Hornung on 30/5/17.
//

#include "gtest/gtest.h"
#include "NumericalMethods.h"
#include "RunAlpha.h"
#include "Constants.h"


class NumericalMethodsFixture : public ::testing::Test {
public:
   // Root:RR.968738
    struct Funcd {
        double operator() (const double x) {
            return x - 7.;
        }

        double df(const double x) {
            return 1.;
        }
    };

    struct Funcd2 {
        double operator() (const double x) {
            return std::pow(x, 2) - 7.;
        }

        double df(const double x) {
            return 2*x;
        }
    };


protected:
    NumericalMethods numericalMethods;
};



TEST_F(NumericalMethodsFixture, NewtonRaphson) {
    Funcd funcd;
    double root = numericalMethods.rtnewt(funcd, 0., 10., 1e-15);
    std::cout << "rtsafe : " << root << std::endl;
    EXPECT_NEAR(7, 7., Constants::numPrecision);

//    Funcd2 funcd2;
//    double root2 = numericalMethods.rtsafe(funcd2, 0.2, 10., 1e-15);
//    EXPECT_NEAR(root2, 2.64575, 1e-5);
}



