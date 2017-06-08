//
// Created by Dirk Hornung on 30/5/17.
//

#include "gtest/gtest.h"
#include "Constants.h"


class ConstantsFixture : public ::testing::Test {
protected:
    Constants constants = Constants(3);
};

TEST_F(ConstantsFixture, zetaFunction) {
    EXPECT_EQ(Constants::zetaFunction(3), 1.202056903159594);
    EXPECT_EQ(Constants::zetaFunction(4), 1.082323233711138);
    EXPECT_EQ(Constants::zetaFunction(5), 1.036927755143370);
}

TEST_F(ConstantsFixture, betaFunction) {
    EXPECT_EQ(Constants::beta(1), 4.5);
    EXPECT_EQ(Constants::beta(2), 8.);
    EXPECT_EQ(Constants::beta(3), 3863./192.);
    EXPECT_EQ(Constants::beta(4), 94.45607914690399);
    EXPECT_NEAR(Constants::beta(5), 254.6443492569230, Constants::numPrecision);
}

TEST_F(ConstantsFixture, adlerCoefficients) {
    EXPECT_NEAR(Constants::c(2,1), 1.63982120489698, 1e-14);
    EXPECT_NEAR(Constants::c(2,2), -1.125, 1e-15);

    EXPECT_NEAR(Constants::c(3,1), 6.37101448310094, 1e-13);
    EXPECT_NEAR(Constants::c(3,2), -5.68959771101822, 1e-14);
    EXPECT_NEAR(Constants::c(3,3), 1.6875, 1e-15);

    EXPECT_NEAR(Constants::c(4,1), 49.07570000294799, 1e-12);
    EXPECT_NEAR(Constants::c(4,2), -33.0914066167203, 1e-12);
    EXPECT_NEAR(Constants::c(4,3), 15.8015948497910, 1e-14);
    EXPECT_NEAR(Constants::c(4,4), -2.84765625, 1e-14);
    EXPECT_NEAR(Constants::c(5,2), -299.177187205152, 1e-11);
    EXPECT_NEAR(Constants::c(5,3), 129.577532569234, 1e-12);
    EXPECT_NEAR(Constants::c(5,4), -40.6160884120297, 1e-13);
    EXPECT_NEAR(Constants::c(5,5), 5.12578125, 1e-15);


}
