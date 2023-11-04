#include <cmath>
#include <iostream>

#include "time_experiment.hh"

#include "vcl/vectorclass.h"
#include "vcl/vectormath_exp.h"

int n = 64;  // NOLINT(readability-magic-numbers)

float func1(float x) {
    return powf(x, 3) - (2 * powf(x, 2)) + (3 * x) - 1;
}

float func2(float x) {
    float sum = 0;
    for (int i = 0; i <= 15; ++i) {  // NOLINT(readability-magic-numbers)
        sum += powf(x, i);           // NOLINT(cppcoreguidelines-narrowing-conversions)
    }
    return sum;
}

float func_1_vectorized(float x) {
    auto const bases = Vec4f(x, x, x, x);
    auto const powers = Vec4f(3, 2, 1, 0);  // NOLINT(readability-magic-numbers)
    auto const results = pow(bases, powers);
    auto const scalars = Vec4f(1, -2, 3, -1);  // NOLINT(readability-magic-numbers)
    return horizontal_add(results * scalars);
}

float func_2_vectorized(float x) {
    auto const base1 = Vec8f(x);
    auto const base2 = Vec8f(x);
    auto const powers1 = Vec8f(0, 1, 2, 3, 4, 5, 6, 7);        // NOLINT(readability-magic-numbers)
    auto const powers2 = Vec8f(8, 9, 10, 11, 12, 13, 14, 15);  // NOLINT(readability-magic-numbers)
    auto results1 = pow(base1, powers1);
    auto results2 = pow(base2, powers2);
    return horizontal_add(results1 + results2);
}

class Experiment {
    int _n;
    float (*_func)(float);

   public:
    Experiment(int n, float (*func)(float)) : _n(n), _func(func) {}

    [[nodiscard]] float midpoint(int a, int b) const {
        float const step_size =
            (float((b - a)) / _n);  // NOLINT(cppcoreguidelines-narrowing-conversions)
        float sum = 0;
        auto x_i = float(a);
        for (int i = 0; i < _n; ++i) {
            float const midpoint = x_i + (step_size / 2);
            sum += _func(midpoint);
            x_i += step_size;
        }
        return step_size * sum;
    }

    void run() const { midpoint(0, 1); }
};

int main() {
    for (int i = 8; i <= 64; i *= 2) {  // NOLINT(readability-magic-numbers)
        Experiment const e1(i, &func1);
        Experiment const e2(i, &func2);
        Experiment const e1_vec(i, &func_1_vectorized);
        Experiment const e2_vec(i, &func_2_vectorized);
        auto d1 = time_experiment(e1);
        auto d2 = time_experiment(e2);
        auto d1_vec = time_experiment(e1_vec);
        auto d2_vec = time_experiment(e2_vec);
        double const time_per_experiment_1 = d1.second * 1.0e-6 / d1.first * 1e9;
        double const time_per_experiment_2 = d2.second * 1.0e-6 / d2.first * 1e9;
        double const time_per_experiment_1_vec = d1_vec.second * 1.0e-6 / d1_vec.first * 1e9;
        double const time_per_experiment_2_vec = d2_vec.second * 1.0e-6 / d2_vec.first * 1e9;

        std::cout << "N = " << i << ":\n";
        std::cout << "func1: " << time_per_experiment_1 << " ns\n";
        std::cout << "func2: " << time_per_experiment_2 << " ns\n";
        std::cout << "func1_vec: " << time_per_experiment_1_vec << " ns\n";
        std::cout << "func2_vec: " << time_per_experiment_2_vec << " ns\n";
    }
}