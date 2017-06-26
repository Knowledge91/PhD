#include <iostream>
#include <vector>
#include "Constants.h"
#include "IntegralMomentum.h"
#include "experimentalData/ExperimentalData.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/RosenBrockFCN.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"


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
    ROOT::Minuit2::MnUserParameters userParameters;
    //MnUserParameters.Add(name, initial Value, STEP)
    userParameters.Add("x", -1., 0.01);
    userParameters.Add("y", 1.2, 0.01);

    RosenBrockFCN rosenBrockFCN;
    ROOT::Minuit2::MnMigrad migrad(rosenBrockFCN, userParameters);
    migrad.Minimizer();

    return 0;
}