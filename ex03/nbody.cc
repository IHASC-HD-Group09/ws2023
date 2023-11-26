#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "nbody_generate.hh"
#include "nbody_io.hh"
#include "time_experiment.hh"
#include "vcl/vectorclass.h"

#ifndef __GNUC__
#define __restrict__
#endif

// basic data type for position, velocity, acceleration
const int M = 4;
using double3 = double[M];  // pad up for later use with SIMD
const int B = 64;           // block size for tiling

/*const double gamma = 6.674E-11;*/
const double G = 1.0;
const double epsilon2 = 1E-10;

void acceleration(int n, double3* __restrict__ x, double* __restrict__ m, double3* __restrict__ a) {
    Vec4d pos_i1, pos_i2, pos_i3, pos_i4;
    Vec4d pos_j1, pos_j2, pos_j3, pos_j4;
    Vec4d acc_i1, acc_i2, acc_i3, acc_i4;
    Vec4d acc_j1, acc_j2, acc_j3, acc_j4;
    Vec4d diff_pos11, diff_pos12, diff_pos13, diff_pos14, diff_pos21, diff_pos22, diff_pos23,
        diff_pos24, diff_pos31, diff_pos32, diff_pos33, diff_pos34, diff_pos41, diff_pos42,
        diff_pos43, diff_pos44;
    Vec4d diff_pos_i12, diff_pos_i13, diff_pos_i14, diff_pos_i23, diff_pos_i24, diff_pos_i34;

    for (int i = 0; i < n; i += 4) {
        pos_i1.load(&x[i][0]);
        acc_i1.load(&a[i][0]);
        pos_i2.load(&x[i + 1][0]);
        acc_i2.load(&a[i + 1][0]);
        pos_i3.load(&x[i + 2][0]);
        acc_i3.load(&a[i + 2][0]);
        pos_i4.load(&x[i + 3][0]);
        acc_i4.load(&a[i + 3][0]);

        diff_pos_i12 = pos_i2 - pos_i1;
        diff_pos_i13 = pos_i3 - pos_i1;
        diff_pos_i14 = pos_i4 - pos_i1;
        diff_pos_i23 = pos_i3 - pos_i2;
        diff_pos_i24 = pos_i4 - pos_i2;
        diff_pos_i34 = pos_i4 - pos_i3;

        double r2_i12 = horizontal_add(diff_pos_i12 * diff_pos_i12) + epsilon2;
        double r2_i13 = horizontal_add(diff_pos_i13 * diff_pos_i13) + epsilon2;
        double r2_i14 = horizontal_add(diff_pos_i14 * diff_pos_i14) + epsilon2;
        double r2_i23 = horizontal_add(diff_pos_i23 * diff_pos_i23) + epsilon2;
        double r2_i24 = horizontal_add(diff_pos_i24 * diff_pos_i24) + epsilon2;
        double r2_i34 = horizontal_add(diff_pos_i34 * diff_pos_i34) + epsilon2;
        double r_i12 = sqrt(r2_i12);
        double r_i13 = sqrt(r2_i13);
        double r_i14 = sqrt(r2_i14);
        double r_i23 = sqrt(r2_i23);
        double r_i24 = sqrt(r2_i24);
        double r_i34 = sqrt(r2_i34);

        double invfact_i12 = G / (r_i12 * r2_i12);
        double invfact_i13 = G / (r_i13 * r2_i13);
        double invfact_i14 = G / (r_i14 * r2_i14);
        double invfact_i23 = G / (r_i23 * r2_i23);
        double invfact_i24 = G / (r_i24 * r2_i24);
        double invfact_i34 = G / (r_i34 * r2_i34);

        double m_i1 = m[i];
        double m_i2 = m[i + 1];
        double m_i3 = m[i + 2];
        double m_i4 = m[i + 3];

        // i1 -> i2, i3, i4
        acc_i1 += m_i2 * invfact_i12 * diff_pos_i12;
        acc_i2 -= m_i1 * invfact_i12 * diff_pos_i12;
        acc_i1 += m_i3 * invfact_i13 * diff_pos_i13;
        acc_i3 -= m_i1 * invfact_i13 * diff_pos_i13;
        acc_i1 += m_i4 * invfact_i14 * diff_pos_i14;
        acc_i4 -= m_i1 * invfact_i14 * diff_pos_i14;
        // i2 -> i3, i4
        acc_i2 += m_i3 * invfact_i23 * diff_pos_i23;
        acc_i3 -= m_i2 * invfact_i23 * diff_pos_i23;
        acc_i2 += m_i4 * invfact_i24 * diff_pos_i24;
        acc_i4 -= m_i2 * invfact_i24 * diff_pos_i24;
        // i3 -> i4
        acc_i3 += m_i4 * invfact_i34 * diff_pos_i34;
        acc_i4 -= m_i3 * invfact_i34 * diff_pos_i34;

        for (int j = i + 4; j < n; j += 2) {
            pos_j1.load(&x[j][0]);
            acc_j1.load(&a[j][0]);
            pos_j2.load(&x[j + 1][0]);
            acc_j2.load(&a[j + 1][0]);
            diff_pos11 = pos_j1 - pos_i1;
            diff_pos12 = pos_j2 - pos_i1;
            diff_pos21 = pos_j1 - pos_i2;
            diff_pos22 = pos_j2 - pos_i2;
            diff_pos31 = pos_j1 - pos_i3;
            diff_pos32 = pos_j2 - pos_i3;
            diff_pos41 = pos_j1 - pos_i4;
            diff_pos42 = pos_j2 - pos_i4;

            double r2_11 = horizontal_add(diff_pos11 * diff_pos11) + epsilon2;
            double r2_12 = horizontal_add(diff_pos12 * diff_pos12) + epsilon2;
            double r2_21 = horizontal_add(diff_pos21 * diff_pos21) + epsilon2;
            double r2_22 = horizontal_add(diff_pos22 * diff_pos22) + epsilon2;
            double r2_31 = horizontal_add(diff_pos31 * diff_pos31) + epsilon2;
            double r2_32 = horizontal_add(diff_pos32 * diff_pos32) + epsilon2;
            double r2_41 = horizontal_add(diff_pos41 * diff_pos41) + epsilon2;
            double r2_42 = horizontal_add(diff_pos42 * diff_pos42) + epsilon2;

            double r_11 = sqrt(r2_11);
            double r_12 = sqrt(r2_12);
            double r_21 = sqrt(r2_21);
            double r_22 = sqrt(r2_22);
            double r_31 = sqrt(r2_31);
            double r_32 = sqrt(r2_32);
            double r_41 = sqrt(r2_41);
            double r_42 = sqrt(r2_42);

            double invfact_11 = G / (r_11 * r2_11);
            double invfact_12 = G / (r_12 * r2_12);
            double invfact_21 = G / (r_21 * r2_21);
            double invfact_22 = G / (r_22 * r2_22);
            double invfact_31 = G / (r_31 * r2_31);
            double invfact_32 = G / (r_32 * r2_32);
            double invfact_41 = G / (r_41 * r2_41);
            double invfact_42 = G / (r_42 * r2_42);

            double m_j1 = m[j];
            double m_j2 = m[j + 1];
            double factorj1 = m_j1 * invfact_11;
            double factorj2 = m_j1 * invfact_12;
            double factorj3 = m_j2 * invfact_21;
            double factorj4 = m_j2 * invfact_22;
            double factorj5 = m_j1 * invfact_31;
            double factorj6 = m_j1 * invfact_32;
            double factorj7 = m_j2 * invfact_41;
            double factorj8 = m_j2 * invfact_42;

            acc_i1 += factorj1 * diff_pos11;
            acc_i1 += factorj2 * diff_pos12;
            acc_i2 += factorj3 * diff_pos21;
            acc_i2 += factorj4 * diff_pos22;
            acc_i3 += factorj5 * diff_pos31;
            acc_i3 += factorj6 * diff_pos32;
            acc_i4 += factorj7 * diff_pos41;
            acc_i4 += factorj8 * diff_pos42;

            acc_j1 -= factorj1 * diff_pos11;
            acc_j1 -= factorj3 * diff_pos21;
            acc_j2 -= factorj2 * diff_pos12;
            acc_j2 -= factorj4 * diff_pos22;
            acc_j1 -= factorj5 * diff_pos31;
            acc_j1 -= factorj7 * diff_pos41;
            acc_j2 -= factorj6 * diff_pos32;
            acc_j2 -= factorj8 * diff_pos42;

            acc_j1.store(&a[j][0]);
            acc_j2.store(&a[j + 1][0]);
        }
        acc_i1.store(&a[i][0]);
        acc_i2.store(&a[i + 1][0]);
        acc_i3.store(&a[i + 2][0]);
        acc_i4.store(&a[i + 3][0]);
    }
}

/** \brief do one time step with leapfrog
 *
 * does n*(n-1)*13 + 12n flops
 */
void leapfrog(int n,
              double dt,
              double3* __restrict__ x,
              double3* __restrict__ v,
              double* __restrict__ m,
              double3* __restrict__ a) {
    // update position: 6n flops
    for (int i = 0; i < n; i++) {
        x[i][0] += dt * v[i][0];
        x[i][1] += dt * v[i][1];
        x[i][2] += dt * v[i][2];
    }

    // save and clear acceleration
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < M; j++) {
            a[i][j] = 0.0;
        }
    }

    // compute new acceleration: n*(n-1)*13 flops
    acceleration(n, x, m, a);

    // update velocity: 6n flops
    for (int i = 0; i < n; i++) {
        v[i][0] += dt * a[i][0];
        v[i][1] += dt * a[i][1];
        v[i][2] += dt * a[i][2];
    }
}

template <typename T>
size_t alignment(const T* p) {
    for (size_t m = 64; m > 1; m /= 2) {
        if (((size_t)p) % m == 0) {
            return m;
        }
    }
    return 1;
}

int main(int argc, char** argv) {
    int n;               // number of bodies in the system
    double* m;           // array for maasses
    double3* x;          // array for positions
    double3* v;          // array for velocites
    double3* a;          // array for accelerations
    int timesteps;       // final time step number
    int k;               // time step number
    int mod;             // files are written when k is a multiple of mod
    char basename[256];  // common part of file name
    char name[256];      // filename with number
    FILE* file;          // C style file hande
    double t;            // current time
    double dt;           // time step

    // command line for restarting
    if (argc == 5) {
        sscanf(argv[1], "%s", &basename);
        sscanf(argv[2], "%d", &k);
        sscanf(argv[3], "%d", &timesteps);
        sscanf(argv[4], "%d", &mod);
    } else if (argc == 6)  // command line for starting with initial condition
    {
        sscanf(argv[1], "%s", &basename);
        sscanf(argv[2], "%d", &n);
        sscanf(argv[3], "%d", &timesteps);
        sscanf(argv[4], "%lg", &dt);
        sscanf(argv[5], "%d", &mod);
    } else  // invalid command line, print usage
    {
        std::cout << "usage: " << std::endl;
        std::cout << "nbody_vanilla <basename> <load step> <final step> <every>" << std::endl;
        std::cout << "nbody_vanilla <basename> <nbodies> <timesteps> <timestep> <every>"
                  << std::endl;
        return 1;
    }

    // set up computation from file
    if (argc == 5) {
        sprintf(name, "%s_%06d.vtk", basename, k);
        file = fopen(name, "r");
        if (file == NULL) {
            std::cout << "could not open file " << std::string(basename) << " aborting"
                      << std::endl;
            return 1;
        }
        n = get_vtk_numbodies(file);
        rewind(file);
        x = new (std::align_val_t(64)) double3[n];
        v = new (std::align_val_t(64)) double3[n];
        m = new (std::align_val_t(64)) double[n];
        read_vtk_file_double(file, n, x, v, m, &t, &dt);
        fclose(file);
        k *= mod;  // adjust step number
        std::cout << "loaded " << n << "bodies from file " << std::string(basename) << std::endl;
    }
    // set up computation from initial condition
    if (argc == 6) {
        x = new (std::align_val_t(64)) double3[n];
        v = new (std::align_val_t(64)) double3[n];
        m = new (std::align_val_t(64)) double[n];
        // plummer(n,17,x,v,m);
        two_plummer(n, 17, x, v, m);
        // cube(n,17,1.0,100.0,0.1,x,v,m);
        std::cout << "initialized " << n << " bodies" << std::endl;
        k = 0;
        t = 0.0;
        printf("writing %s_%06d.vtk \n", basename, k);
        sprintf(name, "%s_%06d.vtk", basename, k);
        file = fopen(name, "w");
        write_vtk_file_double(file, n, x, v, m, t, dt);
        fclose(file);
    }
    if (n % B != 0) {
        std::cout << n << " is not a multiple of the block size " << B << std::endl;
        exit(1);
    }
    if (B % 4 != 0) {
        std::cout << B << "=B is not a multiple of 4 " << std::endl;
        exit(1);
    }
    if (M != 4) {
        std::cout << M << "=M is not 4 " << std::endl;
        exit(1);
    }

    // allocate acceleration vector
    a = new (std::align_val_t(64)) double3[n];
    // explicitly fill/clear padded values
    for (int i = 0; i < n; i++) {
        for (int j = 3; j < M; j++) {
            x[i][j] = v[i][j] = a[i][j] = 0.0;
        }
    }

    // report alignment
    std::cout << "x aligned at " << alignment(x) << std::endl;
    std::cout << "v aligned at " << alignment(v) << std::endl;
    std::cout << "a aligned at " << alignment(a) << std::endl;
    std::cout << "m aligned at " << alignment(m) << std::endl;

    // initialize timestep and write first file
    std::cout << "step=" << k << " finalstep=" << timesteps << " time=" << t << " dt=" << dt
              << std::endl;
    auto start = get_time_stamp();

    // do time steps
    k += 1;
    for (; k <= timesteps; k++) {
        leapfrog(n, dt, x, v, m, a);
        t += dt;
        if (k % mod == 0) {
            auto stop = get_time_stamp();
            double elapsed = get_duration_seconds(start, stop);
            double flop = mod * (13.0 * n * (n - 1.0) + 12.0 * n);
            printf("%g seconds for %g ops = %g GFLOPS\n", elapsed, flop, flop / elapsed / 1E9);

            printf("writing %s_%06d.vtk \n", basename, k / mod);
            sprintf(name, "%s_%06d.vtk", basename, k / mod);
            file = fopen(name, "w");
            write_vtk_file_double(file, n, x, v, m, t, dt);
            fclose(file);

            start = get_time_stamp();
        }
    }

    return 0;
}
