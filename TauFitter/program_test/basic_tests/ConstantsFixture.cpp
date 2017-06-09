//
// Created by Dirk Hornung on 30/5/17.
//

#include "gtest/gtest.h"
#include "Constants.h"


class ConstantsFixture : public ::testing::Test {
protected:
    Constants constants{3, 3, 4};
    constexpr static double numPrecision = 1e-12;
};

TEST_F(ConstantsFixture, zetaFunction) {
    std::cout << constants.pi << std::endl;
    EXPECT_EQ(constants.zetaFunction(3), 1.202056903159594);
    EXPECT_EQ(constants.zetaFunction(4), 1.082323233711138);
    EXPECT_EQ(constants.zetaFunction(5), 1.036927755143370);
}

TEST_F(ConstantsFixture, betaFunction) {
    EXPECT_EQ(constants.beta(1), 4.5);
    EXPECT_EQ(constants.beta(2), 8.);
    EXPECT_EQ(constants.beta(3), 3863./192.);
    EXPECT_EQ(constants.beta(4), 94.45607914690399);
    EXPECT_NEAR(constants.beta(5), 254.6443492569230, numPrecision);
}

TEST_F(ConstantsFixture, adlerCoefficients) {
    EXPECT_NEAR(constants.c(0,0), -5./3., numPrecision);
    EXPECT_NEAR(constants.c(0,1), 1., 1e-15);

    EXPECT_NEAR(constants.c(1,1), 1., 1e-15);
    EXPECT_NEAR(constants.c(1,2), 0., 1e-15);

    EXPECT_NEAR(constants.c(2,1), 1.63982120489698, 1e-14);
    EXPECT_NEAR(constants.c(2,2), -1.125, 1e-15);
    EXPECT_NEAR(constants.c(2,3), 0, 1e-15);

    EXPECT_NEAR(constants.c(3,1), 6.37101448310094, 1e-13);
    EXPECT_NEAR(constants.c(3,2), -5.68959771101822, 1e-14);
    EXPECT_NEAR(constants.c(3,3), 1.6875, 1e-15);
    EXPECT_NEAR(constants.c(3,4), 0., 1e-15);

    EXPECT_NEAR(constants.c(4,1), 49.07570000294799, 1e-12);
    EXPECT_NEAR(constants.c(4,2), -33.0914066167203, 1e-12);
    EXPECT_NEAR(constants.c(4,3), 15.8015948497910, 1e-14);
    EXPECT_NEAR(constants.c(4,4), -2.84765625, 1e-14);
    EXPECT_NEAR(constants.c(4,5), 0, 1e-15);


    EXPECT_NEAR(constants.c(5,2), -299.177187205152, 1e-11);
    EXPECT_NEAR(constants.c(5,3), 129.577532569234, 1e-12);
    EXPECT_NEAR(constants.c(5,4), -40.6160884120297, 1e-13);
    EXPECT_NEAR(constants.c(5,5), 5.12578125, 1e-15);
    EXPECT_NEAR(constants.c(5,6), 0, 1e-15);
}
