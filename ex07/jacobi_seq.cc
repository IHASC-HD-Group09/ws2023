#include <stdio.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "time_experiment.hh"

#define HAVE_NEON 1

#ifdef HAVE_VCL
#include "vcl/vectorclass.h"
#endif

#ifndef __GNUC__
#define __restrict__
#endif

const int B = 20;

struct GlobalContext {
    // input data
    int n;           // nxn lattice of points
    int iterations;  // number of iterations to do
    double* u;       // the initial guess
    // output data

    GlobalContext(int n_) : n(n_) {}
};

void jacobi_vanilla_kernel(int n, int iterations, double* __restrict__ u) {
    // do iterations
    for (int i = 0; i < iterations; i++) {
        for (int i1 = 1; i1 < n - 1; i1++) {
            for (int i0 = 1; i0 < n - 1; i0++) {
                u[i1 * n + i0] = 0.25
                                 * (u[i1 * n + i0 - n] + u[i1 * n + i0 - 1] + u[i1 * n + i0 + 1]
                                    + u[i1 * n + i0 + n]);
            }
        }
    }
}

void jacobi_vanilla(const std::shared_ptr<GlobalContext>& context) {
    jacobi_vanilla_kernel(context->n, context->iterations, context->u);
}

// main function runs the experiments and outputs results as csv
int main(int argc, char** argv) {
    // read parameters
    int n = 1024;
    int iterations = 1000;
    if (argc != 3) {
        std::cout << "usage: ./jacobi_vanilla <size> <iterations>" << std::endl;
        exit(1);
    }
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &iterations);
    // std::cout << "jacobi_vanilla: n=" << n
    //        << " iterations=" << iterations
    //        << " memory (mbytes)=" << (n*n)*8.0*2.0/1024.0/1024.0
    //        << std::endl;

    // get global context shared by aall threads
    auto context = std::make_shared<GlobalContext>(n);
    context->iterations = iterations;

    // allocate aligned arrays
    context->u = new (std::align_val_t(64)) double[n * n];

    // fill boundary values and initial values
    auto g = [&](int i0, int i1) {
        return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1) ? 0.0 : ((double)(i0 + i1)) / n;
    };

    std::cout << "N,";
    std::cout << "vanilla";


    std::cout << std::endl;
    std::cout << n * n;

    // warmup
    for (int i1 = 0; i1 < n; i1++) {
        for (int i0 = 0; i0 < n; i0++) {
            context->u[i1 * n + i0] = g(i0, i1);
        }
    }
    auto start = get_time_stamp();
    jacobi_vanilla(context);
    auto stop = get_time_stamp();
    double elapsed = get_duration_seconds(start, stop);
    double updates = context->iterations;
    updates *= (n - 2) * (n - 2);

    // vanilla
    for (int i1 = 0; i1 < n; i1++) {
        for (int i0 = 0; i0 < n; i0++) {
            context->u[i1 * n + i0] = g(i0, i1);
        }
    }
    start = get_time_stamp();
    jacobi_vanilla(context);
    stop = get_time_stamp();
    elapsed = get_duration_seconds(start, stop);
    std::cout << "," << updates / elapsed / 1e9;

    std::cout << std::endl;

    // deallocate arrays
    delete[] context->u;

    return 0;
}
