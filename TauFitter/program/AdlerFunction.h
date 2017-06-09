//
// Created by Dirk Hornung on 7/6/17.
//

#ifndef PHD_ADLERFUNCTION_H
#define PHD_ADLERFUNCTION_H


#include <cmath>
#include "RunAlpha.h"

class AdlerFunction: public RunAlpha {

public:
    AdlerFunction() { };

    static double D0(double s, double mu) {
        double result = 0.;
        double factor = nc / 12. / std::pow(pi, 2);
        for(int n=0; n<=order; n++) {
            for(int k=1; k<=n+1; k++) {
                std::cout << "n:k=>" << n << ":" << k << " c(n,k)=" << c(n,k) << " l(s,mu)^{k-1}=" << std::pow(l(s,mu), k-1) << std::endl;
                result += std::pow(runAlpha(mu), n) * k * c(n,k) * std::pow(l(s, mu), k-1);
            }
        }
        return factor*result;
    }

    static double l(double s, double mu) {
        return log(-s / std::pow(mu, 2));
    }
};


#endif //PHD_ADLERFUNCTION_H
