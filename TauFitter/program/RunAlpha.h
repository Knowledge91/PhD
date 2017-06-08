//
// Created by Dirk Hornung on 30/5/17.
//

#ifndef PHD_RUNALPHA_H
#define PHD_RUNALPHA_H

#include <iostream>
#include <cmath>
#include "Constants.h"
#include "NumericalMethods.h"
#include <boost/math/tools/roots.hpp>

class RunAlpha: public Constants {
public:
    /*
     * Constructor: RunAlpha
     * Usage: RunAlpha()
     * -------------------------
     * Set /alpha_{mz} and /mu_{mz} from PDG 2016 p. 29, Gauge & Higgs-Boson summary table
     */
    //RunAlpha()  {}

    /*
     * Method: runAlpha
     * Usage: RunAlpha::runAlpha(double mu)
     * -------------------------
     * Solves and return alpha_s iteratively
     */
     double runAlpha(double mu) {
        double approximatedAs = getApproximatedAlpha(mu)/Constants::pi;
        Beta_functor_deriv beta_functor_deriv = Beta_functor_deriv(mu);
        return Constants::pi * boost::math::tools::newton_raphson_iterate(beta_functor_deriv, approximatedAs, approximatedAs-10, approximatedAs+10,20);
    }


    /*
     * Method: getApproximatedAlpha
     * Usage:  RunAlpha::getApproximatedAlpha(double mu)
     * ---------------------------------------
     * Get an approximated value for \alpha_s for nf=3 (\Lambda = 0.332),
     * forumla (and \Lambda) taken from PDG p.132 (p.146)
     */
    double getApproximatedAlpha(double mu) {
        return (16.71314886448282 + log(std::pow(mu,2))*(20.731909992671614
                + log(std::pow(mu,2)) *(9.23729030880304 + 1.3962634015954634*log(std::pow(mu,2))))
                + (-8.83285552360219 + (-5.7374135345463095 - 1.103220465458144*log(std::pow(mu,2)))*log(std::pow(mu,2)))
                  *log(2.205240620131297 + log(std::pow(mu,2))) 
                + (3.6441027182711707 + 0.8716803677693978*log(std::pow(mu,2)))*std::pow(log(2.205240620131297 + log(std::pow(mu,2))),2) -
                0.6887351053980427*std::pow(log(2.205240620131297 + log(std::pow(mu,2))),3))/std::pow(2.205240620131297 + log(std::pow(mu,2)),4);
    }

private:
    struct Beta_functor_deriv: Constants {
    public:
       Beta_functor_deriv(double mu) : mu(mu) {}
        std::pair<double, double> operator()(double const& x) {
            double fx = betaIntegral(alphaMz/pi, x, mz, mu);
            double dx = betaDerivative(x);
//            std::cout << "fx : " << fx << " dx : " << dx << std::endl;
            return std::make_pair(fx, dx);
        }

        // \int_{as1}^{as2}  1/betaFunction(alpha_s) - to the fith order, evaluated with Mathematica (runningAlpha.nb)
        double betaIntegral(double as1, double as2, double mu1, double mu2) {
            return (7493*atan((768 + 3863.*as1)/(12.*sqrt(19082))))/(162.*sqrt(19082)) - (7493*atan((768 + 3863*as2)/(12.*sqrt(19082))))/(162.*sqrt(19082)) +
                   (2*(9 + 16*as1*log(as1) - 8*as1*log(864 + 1536*as1 + 3863*std::pow(as1,2))))/(81.*as1) - (2*(9 + 16*as2*log(as2) - 8*as2*log(864 + 1536*as2 + 3863*std::pow(as2,2))))/(81.*as2)
                   - log(mu1/mu2);
        }
        // Derivative needed for Raphson-Netwon Root finder - evaluated with Mathematica (runningAlpha.nb)
        // D[\int_{as1}^{as2}-Log(mu1, mu2), as2]
        double betaDerivative(double as2) {
            return  192./(std::pow(as2,2)*(864 + as2*(1536 + 3863*as2)));
        }

    private:
        double mu;
    };
};


#endif //PHD_RUNALPHA_H
