//
// Created by Dirk Hornung on 13/6/17.
//

#ifndef PHD_DATA_H
#define PHD_DATA_H

#include <iostream>
#include <map>
#include <vector>

class ExperimentalData {
public:
    /*
     * Initializer:
     * Usage: ExperimentalData experimentalData()
     * ------------------------------------------
     * Fills the AlephData struct with data from experimentalData/aleph14_vpa.f90
     */
    ExperimentalData();


    double vectorPlusAxialvectorIntegralMomentum(double s0);

    /*
     * Methods :
     * Usage: exportDataPoints()
     * -------------------------
     * Exports table of energy and smf2 value into data.dat file ( s | sfm2 ).
     */
    void exportDataPoints();

private:
    struct AlephData {
    public:
        std::vector<double> s;          // mass squared bin center in GeV^2
        std::vector<double> binWidth;   // bin size in Gev^2
        std::vector<double> sfm2;       // normalized invariant mass squared distribution
    };

    void initData(double maxEnergy);

    AlephData alephData;
};


#endif //PHD_DATA_H
