//
// Created by Dirk Hornung on 22/6/17.
//

#ifndef PHD_ROSENBROCKFCN_H
#define PHD_ROSENBROCKFCN_H

#include <vector>
#include <cmath>
#include "Minuit2/FCNBase.h"

class RosenBrockFCN : public ROOT::Minuit2::FCNBase {
public:
    RosenBrockFCN() {}
    ~RosenBrockFCN() {}

    virtual double Up() const { return 1.; }
    virtual double operator()(const std::vector<double>& parameter) const {
        double x = parameter[0];
        double y = parameter[1];
        // RosenBrockFCN Function
        return 100*std::pow((y-x*x),2) + std::pow((1-x), 2);
    }
};

#endif //PHD_ROSENBROCKFCN_H
