//
// Created by Dirk Hornung on 13/6/17.
//

#include <cmath>
#include <fstream>
#include <vector>
#include "ExperimentalData.h"
#include "NumericalMethods.h"

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
        alephData.sbin.push_back(sbin[i]);
        alephData.binWidth.push_back(dsbin[i]);
        alephData.sfm2.push_back(sfm2[i]);

        alephData.cor.resize(80);
        for(int j = 0; j < 80; j++) {
            alephData.cor[i].push_back(corerr[i][j]);
//            std::cout << alephData.cor[i][j] << " ";
        }
//        std::cout << std::endl;
    }
}

double ExperimentalData::integralMomentum(double s0) {
    int nMax = 80;

    // Set Nmax to next bigger sbin
    for (int i = 80; i > 0; i--) {
        //std::cout << sbin[i] << std::endl;
        if(fabs(s0-alephData.sbin[i]) < alephData.binWidth[i]/2.) {
            nMax = i;
            // std::cout << "Nmax: " << Nmax << std::endl;
            break;
        }
    }

    // Integrate momenta up to s0
    double mom = 0.;

    double prefactor = Constants::sTau() / 12. / std::pow(Constants::pi(), 2) / 100. / Constants::Be();

    // prepare 1/wTau(s/sTau) scale
    auto weight = [&](double s) -> double {
        return 1./Constants::wTau(s/Constants::sTau());
    };

    auto fkt = [](double x) {
        return std::pow((1-x),2)*(1.+2.*x);
    };

    for (int i=0; i <= nMax; i++) {
        double lowerLimit = alephData.sbin[i] - alephData.binWidth[i]/2.;
        double upperLimit = alephData.sbin[i] + alephData.binWidth[i]/2.;
        double weightInt = NumericalMethods::integrate(weight, lowerLimit, upperLimit);
        mom += weightInt*alephData.sfm2[i];
    }
    return prefactor*mom;
}

void ExperimentalData::exportDataPoints() {
    std::ofstream dataFile;
    dataFile.open(Constants::generatedPath()+"data.tsv");
    for (int i = 0; i < 80; i++) {
        dataFile << alephData.sbin[i] << "\t" << alephData.sfm2[i] << "\t" << alephData.binWidth[i]  << "\n";
    }
    dataFile.close();
}
