//
// Created by Dirk Hornung on 25/7/17.
//

#ifndef PHD_CHISQUARED_H
#define PHD_CHISQUARED_H

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/RosenBrockFCN.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

class Chisquared {
    void start() {

    }

    double RosenBrock(const double *xx ) {
        const double x = xx[0];
        const double y = xx[1];
        const double tmp1 = y-x*x;
        const double tmp2 = 1-x;
        return 100*tmp1*tmp1+tmp2*tmp2;
    }
};


#endif //PHD_CHISQUARED_H
