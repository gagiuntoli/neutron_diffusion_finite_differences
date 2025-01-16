#include <iostream>

#include "ellpack.h"

int main() {
    size_t nx = 50;
    size_t ny = 50;
    size_t nz = 50;
    double lx = 50.0;
    double ly = 50.0;
    double lz = 50.0;
    double hx = lx / (nx - 1);
    double hy = ly / (ny - 1);
    double hz = lz / (nz - 1);
    double d = 1.2;
    double xs_a = 0.03;
    double xs_f = 0.04;
    double nu = 1.0;
    size_t npoints = nx * ny * nz;

    Ellpack A(npoints, npoints, 7);
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                double ax = -d / (hx * hx);
                double ay = -d / (hy * hy);
                double az = -d / (hz * hz);
                double b = 2 * d / (hx * hx) + 2 * d / (hy * hy) + 2 * d / (hz * hz) + xs_a;

                size_t row = k * nx * ny + j * nx + i;
                A.insert(row, row - 1, ax);
                A.insert(row, row, b);
                A.insert(row, row + 1, ax);

                A.insert(row, row - nx, ay);
                A.insert(row, row + nx, ay);

                A.insert(row, row - nx * ny, az);
                A.insert(row, row + nx * ny, az);
            }
        }
    }
    Ellpack F(npoints, npoints, 1);
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                size_t row = k * nx * ny + j * nx + i;
                double a = nu * xs_f;
                F.insert(row, row, a);
            }
        }
    }

    std::vector<size_t> dirichlet;
    // k=0
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            size_t row = 0 * nx * ny + j * nx + i;
            A.deleteRow(row);
            A.insert(row, row, 1.0);
            F.deleteRow(row);
            F.insert(row, row, 1.0);
            dirichlet.push_back(row);
        }
    }
    // k=nz-1
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            size_t row = (nz - 1) * nx * ny + j * nx + i;
            A.deleteRow(row);
            A.insert(row, row, 1.0);
            F.deleteRow(row);
            F.insert(row, row, 1.0);
            dirichlet.push_back(row);
        }
    }
    // j=0
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            size_t row = k * nx * ny + 0 * nx + i;
            A.deleteRow(row);
            A.insert(row, row, 1.0);
            F.deleteRow(row);
            F.insert(row, row, 1.0);
            dirichlet.push_back(row);
        }
    }
    // j=ny-1
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            size_t row = k * nx * ny + (ny - 1) * nx + i;
            A.deleteRow(row);
            A.insert(row, row, 1.0);
            F.deleteRow(row);
            F.insert(row, row, 1.0);
            dirichlet.push_back(row);
        }
    }
    // i=0
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            size_t row = k * nx * ny + j * nx + 0;
            A.deleteRow(row);
            A.insert(row, row, 1.0);
            F.deleteRow(row);
            F.insert(row, row, 1.0);
            dirichlet.push_back(row);
        }
    }
    // i=nx-1
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            size_t row = k * nx * ny + j * nx + nx - 1;
            A.deleteRow(row);
            A.insert(row, row, 1.0);
            F.deleteRow(row);
            F.insert(row, row, 1.0);
            dirichlet.push_back(row);
        }
    }

    std::vector<double> source(npoints);
    std::vector<double> source_new(npoints);
    std::vector<double> phi(npoints, 1.0);
    double keff = 1.0;
    int iters = 0;

    while (true) {
        for (const auto &node : dirichlet) {
            phi[node] = 0.0;
        }

        // s = F * phi
        F.mvp(source, phi);
        for (auto &val : source) {
            val /= keff;
        }

        // A phi = 1 / keff * s
        A.solve_cg(phi, source);

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
    // for (int k = 0; k < nz; k++) {
    //     for (int j = 0; j < ny; j++) {
    //         for (int i = 0; i < nx; i++) {
    //             size_t row = k * nx * ny + j * nx + i;
    //             std::cout << i * hx << " " << j * hy << " " << k * hz << " " << phi[row] << std::endl;
    //         }
    //     }
    // }

    return 0;
}