//
// Created by Dirk Hornung on 30/5/17.
//

#ifndef PHD_NUMERICALMETHODS_H
#define PHD_NUMERICALMETHODS_H

#include <cmath>


class NumericalMethods {
public:
    // Finds Root to solve non-linear equation
    // based on "Netwon-Raphson" method
    static double NewtonRaphson(double xi, double (*func)(double), double (*Dfunc)(double));

    /*
     * Numerical Recipes: p. 460, Chapter 9. Root Finding and Nonlinear Sets of Equations
     * Using the Newton-Raphson method, return the root of a function known to lie in the interval Œx1;x2 .
     * The root will be refined until its accuracy is known within  ̇xacc.
     * funcd is a user- supplied struct that returns the function value as a functor and the first derivative
     * of the function at the point x as the function df (see text).
    */
    template<class T>
    double rtnewt(T &funcd, const double x1, const double x2, const double xacc) {
        const int JMAX = 20;
        double rtn = 0.5 * (x1 + x2);
        for (int j = 0; j < JMAX; j++) {
            double f = funcd(rtn);
            double df = funcd.df(rtn);
            double dx = f / df;
            rtn -= dx;
            if ((x1 - rtn) * (rtn - x2) < 0.0)
                throw ("Jumped out of brackets in rtnewt");
            if (std::abs(dx) < xacc) return rtn;
        }
        throw ("Maximum number of iterations exceeded in rtnewt");
    }

    template<class T>
    double rtsafe(T &funcd, const double x1, const double x2, const double xacc) {
        const int MAXIT = 10000;
        double xh, xl;
        double fl = funcd(x1);
        double fh = funcd(x2);
        if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) throw ("Root must be bracketed in rtsafe");
        if (fl == 0.0) return x1;
        if (fh == 0.0) return x2;
        if (fl < 0.0) {
            xl = x1;
            xh = x2;
        } else {
            xh = x1;
            xl = x2;
        }
        double rts = 0.5 * (x1 + x2);
        double dxold = std::abs(x2 - x1);
        double dx = dxold;
        double f = funcd(rts);
        double df = funcd.df(rts);
        for (int j = 0; j < MAXIT; j++) {
            if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) || (std::abs(2.0 * f) > std::abs(dxold * df))) {
                dxold = dx;
                dx = 0.5 * (xh - xl);
                rts = xl + dx;
                if (xl == rts) return rts;
            } else {
                dxold = dx;
                dx = f / df;
                double temp = rts;
                rts -= dx;
                if (temp == rts) return rts;
            }
            if (std::abs(dx) < xacc) return rts;
            double f = funcd(rts);
            double df = funcd.df(rts);
            if (f < 0.0)
                xl = rts;
            else
                xh = rts;
            throw ("Maximum number of iterations exceeded in rtsafe");
        }
    }

};

#endif //PHD_NUMERICALMETHODS_H
