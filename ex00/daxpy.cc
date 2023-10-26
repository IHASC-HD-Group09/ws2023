#include <iostream>
#include <thread>
#include <vector>



// NOLINTBEGIN(readability-identifier-length, readability-magic-numbers)
const unsigned int P = std::thread::hardware_concurrency();  // number of threads
const int N = 12;                                            // problem size
std::vector<double> x(N, 2.0);                               // first vector
std::vector<double> y(N, 2.0);                               // second vector
// NOLINTEND(readability-identifier-length, readability-magic-numbers)
const double ALPHA = 5.0;                                   // scalar alpha

void daxpy(int rank) {
    for (unsigned int i = (N * rank) / P; i < (N * (rank + 1)) / P; ++i) {
        y[i] = x[i] * ALPHA + y[i];
    }
}

int main() {
    std::vector<std::thread> threads;
    threads.reserve(P);
    for (int rank = 0; rank < P; ++rank) {
        threads[rank] = std::thread{daxpy, rank};
    }
    for (int rank = 0; rank < P; ++rank) {
        threads[rank].join();
    }
    std::cout << "y after daxpy: ";
    for (const double i : y) {  // NOLINT(readability-identifier-length)
        std::cout << i << " ";
    }
}
