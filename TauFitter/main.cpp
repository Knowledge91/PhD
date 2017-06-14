#include <iostream>
#include "Constants.h"
#include "IntegralMomentum.h"
#include "experimentalData/ExperimentalData.h"

extern"C" {
double vphlmntv2_(double *energy, double *vprehadsp, double *vprehadtm, double *vpimhad, double *vprelepsp, double *vpreleptm, double *vpimlep, double *vpretopsp, double *vpretoptm, int *nrflag);
}


extern"C" {
void aleph_vplusa_(double *sbin, double *dsbin, double *sfm2, double *derr, double (*corerr)[80]);
}

int main() {
    std::cout << "Tau fitter start..." << std::endl;

    Constants constants(3,3,4);
    IntegralMomentum integralMomentum(constants);

    std::cout << "countour " << integralMomentum.contourIntegral(2.5) << std::endl;

    ExperimentalData experimentalData;
    std::cout << "Exp. Integral Momentum " << experimentalData.vectorPlusAxialvectorIntegralMomentum(1.9) << std::endl;
    experimentalData.exportDataPoints();
//    experimentalData.plotDataPoints();

    return 0;
}