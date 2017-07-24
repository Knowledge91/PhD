//
// Created by Dirk Hornung on 13/6/17.
//

#ifndef PHD_DATA_H
#define PHD_DATA_H

#include <iostream>
#include <map>
#include <vector>

//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/lu.hpp>
//#include <boost/numeric/ublas/io.hpp>

#include "Constants.h"

class ExperimentalData {
public:
    /*
     * Initializer:
     * Usage: ExperimentalData experimentalData()
     * ------------------------------------------
     * Fills the AlephData struct with data from experimentalData/aleph14_vpa.f90
     */
    ExperimentalData();

    /*
     * Method: integralMomentum
     * Usage: integralMomentum(double s0)
     * ------------------------------------------
     * Outputs the experimental integral Momentum for the invariant mass squared s0
     */
    double integralMomentum(double s0);

    /*
     * Method : exportDataPoints()
     * Usage: this->exportDataPoints()
     * -------------------------
     * Exports table of energy and smf2 value into data.dat file ( s | sfm2 ).
     */
    void exportDataPoints();

    /*
     * Method : plotDataPoints
     * Usage: this->plotDataPoints()
     * -------------------------
     * Plots previously exported data points.
     */
    void plotDataPoints() {
        std::string str = "gnuplot> load '"+Constants::programPath()+"/program/experimentalData/plot.gp'";
        const char * c = str.c_str();
        std::system(c);
    }

private:
    struct AlephData {
    public:
        std::vector<double> sbin;       // mass squared bin center in GeV^2
        std::vector<double> binWidth;   // bin size in Gev^2
        std::vector<double> sfm2;       // normalized invariant mass squared distribution
    };

    /*
     *
     */


    void initData(double maxEnergy);

    AlephData alephData;
};


#endif //PHD_DATA_H
