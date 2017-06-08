//
// Created by Dirk Hornung on 7/6/17.
//

#ifndef PHD_ADLERFUNCTION_H
#define PHD_ADLERFUNCTION_H


#include <cmath>
#include "Constants.h"
#include "RunAlpha.h"

class AdlerFunction: public RunAlpha {
public:
    AdlerFunction() : nc(3), order(1) { };

    double D0(double s, double mu) {
        double result = nc / 12. / pi;
        for(int n=0; n<=order; n++) {
            for(int k=1; k<=n+1; k++) {
                result += std::pow(runAlpha(mu), n) * k * c(n,k) * std::pow(l(s, mu), k-1);
            }
        }
        return result;
    }

    double l(double s, double mu) {
        return log(-s / std::pow(mu, 2));
    }



    int nc;
    int order;
};


#endif //PHD_ADLERFUNCTION_H
