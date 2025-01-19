#include <iostream>

#include "ellpack.h"

int main() {
    size_t npoints = 50;
    double length = 50.0;
    double h = length / (npoints - 1);
    double d1 = 1.2;
    double xs_a1 = 0.03;
    double nu1 = 0.0;
    double xs_f1 = 0.0;
    double d2 = 1.2;
    double xs_a2 = 0.03;
    double nu2 = 1.0;
    double xs_f2 = 0.3;
    double xs_s12 = 0.08;
    size_t ngroups = 2;

    Ellpack A(npoints * ngroups, npoints * ngroups, 4);
    for (int i = 0; i < npoints; i++) {
        double a = -d1 / (h * h);
        double b;
        if (i * h > 15.0 && i * h < 35) {  // fuel
            b = 2 * d1 / (h * h) + xs_a1;
        } else {  // water (moderator)
            b = 2 * d1 / (h * h) + xs_a1 + xs_s12;
        }
        A.insert(i, i - 1, a);
        A.insert(i, i, b);
        A.insert(i, i + 1, a);
    }
    for (int i = npoints; i < 2 * npoints; i++) {
        double a = -d2 / (h * h);
        double b = 2 * d2 / (h * h) + xs_a2;
        double c;
        if ((i - npoints) * h > 15.0 && (i - npoints) * h < 35) {  // fuel
            c = 0.0;
        } else {  // water (moderator)
            c = -xs_s12;
        }
        A.insert(i, i - 1, a);
        A.insert(i, i, b);
        A.insert(i, i + 1, a);
        A.insert(i, i - npoints, c);
    }

    A.deleteRow(0);
    A.insert(0, 0, 1.0);
    A.deleteRow(npoints - 1);
    A.insert(npoints - 1, npoints - 1, 1.0);
    A.deleteRow(npoints);
    A.insert(npoints, npoints, 1.0);
    A.deleteRow(2 * npoints - 1);
    A.insert(2 * npoints - 1, 2 * npoints - 1, 1.0);

    Ellpack F(npoints * ngroups, npoints * ngroups, 2);
    for (int i = 0; i < npoints; i++) {
        double a = nu2 * xs_f2;
        if (i * h > 15.0 && i * h < 35) {  // fuel
            a = nu2 * xs_f2;
        } else {  // water (moderator)
            a = 0.0;
        }
        F.insert(i, i + npoints, a);
    }
    F.deleteRow(0);
    F.insert(0, 0, 1.0);
    F.deleteRow(npoints - 1);
    F.insert(npoints - 1, npoints - 1, 1.0);
    F.deleteRow(npoints);
    F.insert(npoints, npoints, 1.0);
    F.deleteRow(2 * npoints - 1);
    F.insert(2 * npoints - 1, 2 * npoints - 1, 1.0);

    std::vector<double> source(npoints * ngroups);
    std::vector<double> source_new(npoints * ngroups);
    std::vector<double> phi(npoints * ngroups, 1.0);
    double keff = 1.0;
    int iters = 0;

    while (true) {
        phi[0] = 0;
        phi[npoints - 1] = 0;
        phi[npoints] = 0;
        phi[2 * npoints - 1] = 0;

        // s = F * phi
        F.mvp(source, phi);
        for (auto &val : source) {
            val /= keff;
        }

        // A phi = 1 / keff * s
        A.solve_jacobi(phi, source);

        // s = F * phi
        F.mvp(source_new, phi);

        // keff = integral (source_new) / integral (source)
        double power_new = 0, power = 0;
        for (int i = 0; i < npoints; i++) {
            power_new += source_new[i];
            power += source[i];
        }
        double keff_new = power_new / power;

        if (std::abs(keff_new - keff) / keff < 0.0001) break;

        keff = keff_new;
        iters++;
    }

    std::cout << "keff: " << keff << " iters: " << iters << std::endl;
    for (int i = 0; i < npoints; i++) {
        std::cout << i * h << " " << phi[i] << " " << phi[npoints + i] << std::endl;
    }

    return 0;
}