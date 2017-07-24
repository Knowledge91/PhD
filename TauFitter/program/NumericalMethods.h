//
// Created by Dirk Hornung on 30/5/17.
//

#ifndef PHD_NUMERICALMETHODS_H
#define PHD_NUMERICALMETHODS_H

#include <cmath>
#include <vector>
#include <iostream>

template<class T>
inline void SWAP(T &a, T &b)
{T dum=a; a=b; b=dum;}

class NumericalMethods {
public:
    // 3rd Numerical Recpies p.184 gauleg
    // Gauss Quadratur
    template <typename T, typename Func>
    static T integrate(Func function, T a, T b) {
        std::vector<T> x(10);
        std::vector<T> w(10);

        T x1 = a;
        T x2 = b;

        // GAULEG weights and abciss nodes
        // initialise integration
        // find abscissas + weights
        const double EPS = 1.0e-14;
        T z1, z, xm, xl, pp, p3, p2, p1;
        int n = x.size();
        int m = (n+1.)/2.;
        xm = 0.5*(x2+x1);
        xl = 0.5*(x2-x1);
        for (int i = 0; i < m; i++) {
            z = cos(3.141592654*(i+0.75)/(n+0.5));
            do {
                p1 = 1.0;
                p2 = 0.0;
                for (int j = 0; j < n; j++) {
                    p3 = p2;
                    p2 = p1;
                    T jComplx(j);
                    p1 = ((2.0*j+1.0)*z*p2-jComplx*p3) / (jComplx+1.);
                }
                T nComplx(n);
                pp = nComplx*(z*p1-p2)/(z*z-1.0);
                z1 = z;
                z = z1 - p1/pp;
            } while (std::abs(z-z1) > EPS);
            x[i] = xm - xl*z;
            x[n-1-i] = xm+xl*z;
            w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
            w[n-1-i] = w[i];
        }


        // qgauss Integrate
        T s = 0;
        for (int i = 0; i < x.size(); i++) {
            s += w[i]*function(x[i]);
        }

        return s;
    }

    /*
     * gaussj:
     * Usage: ExperimentalData experimentalData()
     * ------------------------------------------
     * Fills the AlephData struct with data from experimentalData/aleph14_vpa.f90
     * taken from: Numerical Recipies, Third Edition, p. 44-45
     */
    static void gaussj(std::vector<std::vector<double> > &a, std::vector<std::vector<double> > &b) {
        int i, icol, irow, j, k, l, ll, n = a.size(), m = b[0].size();
        double big, dum, pivinv;

        std::vector<int> indxc(n), indxr(n), ipiv(n);
        for (j=0; j<n; j++) ipiv[j]=0;
        for (i=0; i<n; i++) {
            big = 0.0;
            for (j = 0; j < n; j++)
                if (ipiv[j] != 1)
                    for (k = 0; k < n; k++) {
                        if (ipiv[k] == 0) {
                            if (abs(a[j][k]) >= big) {
                                big = abs(a[j][k]);
                                irow = j;
                                icol = k;
                            }
                        }
                    }
            ++(ipiv[icol]);

            if (irow != icol) {
                for (l = 0; l < n; l++) SWAP(a[irow][l], a[icol][l]);
                for (l = 0; l < m; l++) SWAP(b[irow][l], b[icol][l]);
            }
            indxr[i] = irow;
            indxc[i] = icol;
            if (a[icol][icol] == 0.0) throw ("gaussj: Singular Matrix");
            pivinv = 1.0 / a[icol][icol];
            a[icol][icol] = 1.0;
            for (l = 0; l < n; l++) a[icol][l] *= pivinv;
            for (l = 0; l < m; l++) b[icol][l] *= pivinv;
            for (ll = 0; ll < n; ll++)
                if (ll != icol) {
                    dum = a[ll][icol];
                    a[ll][icol] = 0.0;
                    for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
                    for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
                }
        }

        for (l=n-1; l>=0; l--) {
            if (indxr[l] != indxc[l])
                for (k=0; k<n; k++)
                    SWAP(a[k][indxr[l]], a[k][indxc[l]]);
        }
    }

    static void gaussj(std::vector<std::vector<double> > &a) {
        std::vector<std::vector<double> > b(a.size(), std::vector<double>(a.size()));
        gaussj(a, b);
    }

};

#endif //PHD_NUMERICALMETHODS_H
