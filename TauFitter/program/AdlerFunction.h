//
// Created by Dirk Hornung on 7/6/17.
//

#ifndef PHD_ADLERFUNCTION_H
#define PHD_ADLERFUNCTION_H


#include <cmath>
#include <complex>
#include "RunAlpha.h"
#include "NumericalMethods.h"

class AdlerFunction: public RunAlpha {

public:
    AdlerFunction(Constants constants) : constants(constants) { };

    template <typename T>
    T D0(T s, double mu) {
        T result = 0.;
        T factor = constants.nc / 12. / std::pow(constants.pi(), 2);
        for(int n=0; n<=constants.order; n++) {
            for(int k=1; k<=n+1; k++) {
//              std::cout << "n:k=>" << n << ":" << k << " c(n,k)=" << c(n,k) << " l(s,mu)^{k-1}=" << std::pow(l(s,mu), k-1) << std::endl;
                result += std::pow(runAlpha(mu), n) * double(k) * constants.c(n,k) * std::pow(l(s, mu), k-1);
            }
        }
        return factor*result;
    }

    template <typename T>
    static T l(T s, double mu) {
        return log(-s / std::pow(mu, 2));
    }

private:
    Constants constants;
};


#endif //PHD_ADLERFUNCTION_H
