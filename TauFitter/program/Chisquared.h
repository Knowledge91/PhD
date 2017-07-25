//
// Created by Dirk Hornung on 25/7/17.
//

#ifndef PHD_CHISQUARED_H
#define PHD_CHISQUARED_H

#include "experimentalData/ExperimentalData.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/RosenBrockFCN.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

class Chisquared {
 public:
    void minimize() {
        // MINUIT2 TEST
        ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
        min.SetMaxFunctionCalls(1000000);
        min.SetMaxIterations(100000);
        min.SetTolerance(1e-15);

        ROOT::Math::Functor f(&chisquaredFunction, 2);
        double step[2] = {0.01,0.01};
        double variable[2] = { -1.,1.2};

        min.SetFunction(f);

        // Set the free variables to be minimized!
        min.SetVariable(0,"x",variable[0], step[0]);
        min.SetVariable(1,"y",variable[1], step[1]);

        min.Minimize();

        const double *xs = min.X();
        std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
                  << chisquaredFunction(xs) << std::endl;
    }

    /*
     * Method : chisquaredFunction
     * Usage: Chisquared::chisquaredFunction()
     * -------------------------
     * Prepares the chisquare function to minimize.
     * (I_exp(s) - I_th)^T Cor^-1 (I_exp(s) - I_th)
     * PDG 2016 p. 524 eq. 39.20
     */
    static double chisquaredFunction(const double *xx ) {
        ExperimentalData experimentalData;

        double i_exp = 0, i_th = 0, sbin;
        for(int i = 0; i < 80; i++) {
            sbin = experimentalData.alephData.sbin[i];
            i_exp += experimentalData.integralMomentum(sbin);
        }

        const double x = xx[0];
        const double y = xx[1];

        const double tmp1 = y-x*x;
        const double tmp2 = 1-x;
        return 100*tmp1*tmp1+tmp2*tmp2;
    }

 private:
};


#endif //PHD_CHISQUARED_H
