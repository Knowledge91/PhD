//
// Created by Dirk Hornung on 1/6/17.
//

#include "gtest/gtest.h"
#include "RunAlpha.h"
#include "Constants.h"


class RunAlphaFixture : public ::testing::Test {
protected:
    RunAlpha runAlpha;
};


TEST_F(RunAlphaFixture, getApproximateAlpha) {
    // Compare with Mathematica Package: RunDec, Function: AlphasLam[...]
    EXPECT_NEAR(runAlpha.getApproximatedAlpha(10.), 0.165139, 1e-6);
    EXPECT_NEAR(runAlpha.getApproximatedAlpha(25.), 0.133393, 1e-6);
    EXPECT_NEAR(runAlpha.getApproximatedAlpha(50.), 0.116702, 1e-6);
}

TEST_F(RunAlphaFixture, runAlpha) {
    EXPECT_NEAR(runAlpha.runAlpha(10.), 0.20018369730618291319, 1e-4);
    EXPECT_NEAR(runAlpha.runAlpha(25.), 0.15476396933502079243, 1e-3);
    EXPECT_NEAR(runAlpha.runAlpha(50.), 0.13252265820843999954, 1e-3);
}


