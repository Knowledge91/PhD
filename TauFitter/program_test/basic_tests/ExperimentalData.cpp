//
// Created by Dirk Hornung on 30/5/17.
//

#include "experimentalData/ExperimentalData.h"
#include "gtest/gtest.h"

class ExperimentalDataFixture : public ::testing::Test {
protected:
    ExperimentalData experimentalData;
};



TEST_F(ExperimentalDataFixture, integralMomentum) {
   EXPECT_NEAR(experimentalData.integralMomentum(1.), 2.62334, 1e-1);
   EXPECT_NEAR(experimentalData.integralMomentum(1.5), 3.81741, 1e-1);
   EXPECT_NEAR(experimentalData.integralMomentum(10.), 7.40945, 1e-5);
}



