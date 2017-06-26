#include <iostream>
#include <vector>
#include "Constants.h"
#include "IntegralMomentum.h"
#include "experimentalData/ExperimentalData.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/RosenBrockFCN.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

double RosenBrock(const double *xx )
{
    const double x = xx[0];
    const double y = xx[1];
    const double tmp1 = y-x*x;
    const double tmp2 = 1-x;
    return 100*tmp1*tmp1+tmp2*tmp2;
}

extern"C" {
double vphlmntv2_(double *energy, double *vprehadsp, double *vprehadtm, double *vpimhad, double *vprelepsp, double *vpreleptm, double *vpimlep, double *vpretopsp, double *vpretoptm, int *nrflag);
}
extern"C" {
void aleph_vplusa_(double *sbin, double *dsbin, double *sfm2, double *derr, double (*corerr)[80]);
}

int main() {
    std::cout << "Tau fitter start..." << std::endl;

    Constants constants(3, 3, 4);
    IntegralMomentum integralMomentum(constants);

    std::cout << "countour " << integralMomentum.contourIntegral(2.5) << std::endl;

    ExperimentalData experimentalData;
    std::cout << "Exp. Integral Momentum " << experimentalData.integralMomentum(1.9) << std::endl;
    experimentalData.exportDataPoints();
//    experimentalData.plotDataPoints();


    // MINUIT2 TEST
    ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
    min.SetMaxFunctionCalls(1000000);
    min.SetMaxIterations(100000);
    min.SetTolerance(1e-15);

    ROOT::Math::Functor f(&RosenBrock, 2);
    double step[2] = {0.01,0.01};
    double variable[2] = { -1.,1.2};

    min.SetFunction(f);

    // Set the free variables to be minimized!
    min.SetVariable(0,"x",variable[0], step[0]);
    min.SetVariable(1,"y",variable[1], step[1]);

    min.Minimize();

    const double *xs = min.X();
    std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): "
         << RosenBrock(xs) << std::endl;
    return 0;
}