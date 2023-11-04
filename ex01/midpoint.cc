#include <math.h>
#include <iostream>
#include <immintrin.h>

int n = 64;


float func1(float x) {
    return pow(x, 3) - (2 * pow(x, 2)) + (3 * x) - 1;
}

float func2(float x) {
    float sum = 0;
    for(int i=0; i < 15; ++i) {
        sum += pow (x, i);
    }
    return sum;
}

__m256 func1_vectorized(__m256 x) {
    float result;
    __m256 x_squared = _mm256_mul_ps(x, x);
    __m256 x_cubed = _mm256_mul_ps(x_squared, x);
    __m256 two_x_squared = _mm256_add_ps(x_squared, x_squared);
    __m256 three_x = _mm256_mul_ps(_mm256_set1_ps(3.0f), x);
    return _mm256_sub_ps(x_cubed, _mm256_add_ps(two_x_squared, _mm256_sub_ps(three_x, _mm256_set1_ps(1.0f))));
}

__m256 func2_vectorized(__m256 x) {
    __m256 sum = _mm256_set1_ps(0.0f);
    __m256 x_powers = _mm256_set1_ps(1.0f);

    for (int i = 0; i < 15; ++i) {
        sum = _mm256_add_ps(sum, x_powers);
        x_powers = _mm256_mul_ps(x_powers, x);
    }

    return sum;
}



class Experiment {
    int _n;
    float (*_func)(float);

    public:
     Experiment(int n, float (*func)(float)) : _n(n), _func(func) {

     }

     float midpoint(int a, int b) const {
        float step_size = (float((b - a)) / _n);
        float sum = 0;
        float x_i = float(a);
        for(int i=0; i < _n; ++i) {
            float midpoint = x_i + (step_size / 2);
            sum += _func(midpoint);
            x_i += step_size;
        }
        return step_size * sum;
      }
    
     void run() const {
        std::cout << midpoint(0, 1) << std::endl;
     }
};

class VectorizedExperiment {
    int _n;
    __m256 (*_func)(__m256);

    public:
     VectorizedExperiment(int n, __m256 (*func)(__m256)) : _n(n), _func(func) {

     }

     __m256 midpoint(int a, int b) const {
        __m256 step_size = _mm256_set1_ps(float((b - a)) / _n);
        __m256 sum = _mm256_set1_ps(0.0f);
        __m256 x_i = _mm256_set1_ps(float(a));
        for(int i=0; i < _n; ++i) {
            __m256 midpoint =_mm256_add_ps(x_i, _mm256_div_ps(step_size, _mm256_set1_ps(2.0f)));
            sum = _mm256_add_ps(sum, _func(midpoint));
            x_i = _mm256_add_ps(x_i, step_size);
        }
        return _mm256_mul_ps(step_size, sum);
      }
    
     void run() const {
        float result;
        _mm256_storeu_ps(&result, midpoint(0, 1));
        std::cout << result << std::endl;
     }
};

int main() {
    Experiment const e1(n, &func1);
    Experiment const e2(n, &func2);
    VectorizedExperiment const e3(n, &func1_vectorized);
    e1.run();
    e2.run();
    e3.run();

}