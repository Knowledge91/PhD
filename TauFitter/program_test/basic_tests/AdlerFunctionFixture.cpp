//
// Created by Dirk Hornung on 8/6/17.
//

#include "gtest/gtest.h"
#include "AdlerFunction.h"

class AdlerFunctionFixture : public ::testing::Test {
protected:
    Constants constants{3, 3, 4};
    AdlerFunction adlerFunction{constants};
    constexpr static double numPrec = 1e-6;
};

TEST_F(AdlerFunctionFixture, D0) {
   EXPECT_NEAR(adlerFunction.D0(-2., constants.mz()),  0.1060670 ,numPrec);
   EXPECT_NEAR(adlerFunction.D0(-10., constants.mz()), 0.0773029 , numPrec);
}