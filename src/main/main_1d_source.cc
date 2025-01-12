#include <iostream>

#include "ellpack.h"

int main() {
    size_t npoints = 100;
    double length = 10.0;
    double h = length / (npoints - 1);
    double d = 1.0;
    double xs_a = 1.0;

    Ellpack A(npoints, npoints, 3);
    for (int i = 0; i < npoints; i++) {
        double a = -d / (h * h);
        double b = 2 * d / (h * h) + xs_a;
        A.insert(i, i - 1, a);
        A.insert(i, i, b);
        A.insert(i, i + 1, a);
    }
    A.deleteRow(0);
    A.insert(0, 0, 1.0);
    A.deleteRow(npoints - 1);
    A.insert(npoints - 1, npoints - 1, 1.0);

    std::vector<double> source(npoints, 1.0);
    source[0] = 0.0;
    source[npoints - 1] = 0.0;

    std::vector<double> phi(npoints, 0.0);

    A.solve_cg(phi, source);

    for (int i = 0; i < npoints; i++) {
        std::cout << i * h << " " << phi[i] << std::endl;
    }

    return 0;
}