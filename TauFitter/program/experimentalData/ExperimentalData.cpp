//
// Created by Dirk Hornung on 13/6/17.
//

#include <cmath>
#include <fstream>
#include "ExperimentalData.h"


extern"C" {
double vphlmntv2_(double *energy, double *vprehadsp, double *vprehadtm, double *vpimhad, double *vprelepsp, double *vpreleptm, double *vpimlep, double *vpretopsp, double *vpretoptm, int *nrflag);
}

extern"C" {
void aleph_vplusa_(double *sbin, double *dsbin, double *sfm2, double *derr, double (*corerr)[80]);
}


ExperimentalData::ExperimentalData() {
    double sbin[80], dsbin[80], sfm2[80], derr[80], corerr[80][80];
    aleph_vplusa_(sbin, dsbin, sfm2, derr, corerr);

    for(int i = 0; i < 80; i++) {
        alephData.s.push_back(sbin[i]);
        alephData.binWidth.push_back(dsbin[i]);
        alephData.sfm2.push_back(sfm2[i]);
    }
}

double ExperimentalData::vectorPlusAxialvectorIntegralMomentum(double s0) {
    int Nmax = 0;

    // Set Nmax to next bigger sbin
    for (int i = 80; i > 0; i--) {
        //std::cout << sbin[i] << std::endl;
        if(fabs(s0-alephData.s[i]) < alephData.binWidth[i]/2.) {
            Nmax = i;
            // std::cout << "Nmax: " << Nmax << std::endl;
            break;
        }
    }

    // Integrate momenta up to s0
    double mom = 0.;
    for (int i=1; i < Nmax; i++) {
        mom += alephData.sfm2[i];
    }
    return mom;
}

void ExperimentalData::exportDataPoints() {
    std::ofstream dataFile;
    dataFile.open(Constants::generatedPath()+"data.dat");
    for (int i = 0; i < 80; i++) {
        dataFile << alephData.s[i] << "\t" << alephData.sfm2[i] << "\n";
    }
    dataFile.close();
}
