//
// Created by Dirk Hornung on 30/5/17.
//

#ifndef PHD_CONSTANTS_H
#define PHD_CONSTANTS_H

#include <string>
#include <cmath>

class Constants {
public:
    Constants(int nc, int nf, int order) : nc(nc), nf(nf), order(order) {};

    static double zetaFunction(int i) {
        switch(i) {
            case 3: return 1.202056903159594;
            case 4: return 1.082323233711138;
            case 5: return 1.036927755143370;
            case 7: return 1.008349277381923;
        }
        return 0;
    }

    double beta(int i) {
        switch(i) {
            case 1: return 11./2. - nf/3.;
            case 2: return 51./4. - 19./12.*nf;
            case 3: return 2857. / 64. - 5033. / 576. * nf + 325. / 1728. * std::pow(nf, 2);
            case 4: return 149753./768. + 891./32.*zetaFunction(3)
                           - (1078361./20736. + 1627./864.*zetaFunction(3))
                             * nf + (50065./20736. + 809./1296.*zetaFunction(3))*std::pow(nf, 2)
                           + 1093./93312.*std::pow(nf, 3);
            case 5: return 2./std::pow(4.,5)*(8157455./16. + 621885./2.*zetaFunction(3) - 88209./2.*zetaFunction(4)
                           - 288090.*zetaFunction(5) + nf*(-336460813./1944. - 4811164./81.*zetaFunction(3)
                           + 33935./6.*zetaFunction(4) + 1358995./27.*zetaFunction(5))
                           + std::pow(nf, 2)*(25960913./1944. + 698531./81.*zetaFunction(3)
                           - 10526./9.*zetaFunction(4) - 381760./81.*zetaFunction(5))
                           + std::pow(nf, 3)*(-630559./5832. - 48722./243.*zetaFunction(3)
                           + 1618./27.*zetaFunction(4) + 460./9.*zetaFunction(5))
                           + std::pow(nf, 4)*(1205./2916. - 152./81.*zetaFunction(3)));
        }
        return 0;
    }

    double betaFunktion(double a) {
        double betaFunction;
        for(int i; i<=5; i++) {
            betaFunction += std::pow(a, i+1)*beta(i);
        }
        return betaFunction;
    }

    double c(int n, int k) {
        switch(n) {
            case 0:
                switch(k) {
                    case 0: return - 5./3.;
                    case 1: return 1.;
                    default:
                        return 0;
                }
            case 1:
                switch(k) {
                    case 1: return 1.;
                    case 2: return 0.;
                }
            case 2:
                switch(k) {
                    case 1: return 365./24. - 11*zetaFunction(3) - (11./12. - 2./3. * zetaFunction(3))*nf;
                    case 2: return - 1./4.*beta(1)*c(1,1);
                    case 3: return 0.;
                }
            case 3:
                switch(k) {
                    case 1: return 87029./288. - 1103./4. * zetaFunction(3) + 275./6.*zetaFunction(5)
                                   - ( 7847./216. - 262./9.*zetaFunction(3) + 25./9.*zetaFunction(5) )*nf
                                       + (151./162. - 19./27.*zetaFunction(3))*std::pow(nf, 2);
                    case 2: return 1./4.*(-beta(2)*c(1,1) - 2*beta(1)*c(2,1));
                    case 3: return 1./12.*std::pow(beta(1), 2)*c(1,1);
                    default:
                        return 0.;
                }
            case 4: switch(k) {
                    case 1: return 78631453./20736 - 1704247./432.*zetaFunction(3) + 4185./8.*std::pow(zetaFunction(3), 2)
                                   + 34165./96. * zetaFunction(5) - 1995./16.*zetaFunction(7);
                    case 2: return 1./4.*(-beta(3)*c(1,1)-2*beta(2)*c(2,1)-3*beta(1)*c(3,1));
                    case 3: return 1./24.*(6.*c(2,1)*std::pow(beta(1), 2) + 5*beta(2)*beta(1)*c(1,1));
                    case 4: return - 1./32.*std::pow(beta(1), 3)*c(1,1);
                    default:
                        return 0.;
                }
            case 5: switch(k) {
                    case 2: return 1./4.*(-beta(4)*c(1,1) - 2*beta(3)*c(2,1) - 3*beta(2)*c(3,1) - 4*beta(1)*c(4,1));
                    case 3: return 1./24.*(12.*c(3,1)*std::pow(beta(1), 2) + 6*beta(1)*beta(3)*c(1,1)
                                           + 14*beta(2)*beta(1)*c(2,1) + 3*std::pow(beta(2), 2)*c(1,1));
                    case 4: return 1./96.*(-12.*std::pow(beta(1), 3)*c(2,1) - 13.*beta(2)*std::pow(beta(1),2)*c(1,1));
                    case 5: return 1./80.*std::pow(beta(1), 4)*c(1,1);
                    default:
                        return 0.;
                }
            default:
                return 0;
        }
    }


    int nc;
    int nf;
    int order;
    static double pi() {
        return 3.141592653589793;
    }
    constexpr static double lamba = 0.332;
    static double mz() { return 91.1876; }
    static double alphaMz() { return 0.1181; }
    static std::string programPath() { return "/Users/Knowledge/Developer/PhD/TauFitter/"; }
    static std::string generatedPath() { return programPath()+"/generatedFiles"; }
};


#endif //PHD_CONSTANTS_H
