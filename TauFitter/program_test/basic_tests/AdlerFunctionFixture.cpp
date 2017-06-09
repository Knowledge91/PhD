//
// Created by Dirk Hornung on 8/6/17.
//

#include "gtest/gtest.h"
#include "AdlerFunction.h"

class AdlerFunctionFixture : public ::testing::Test {
protected:
};

TEST_F(AdlerFunctionFixture, D0) {
    EXPECT_NEAR(AdlerFunction::D0(-2., AdlerFunction::mz), , 1e-2);
}