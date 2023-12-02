#include <atomic>
#include <cmath>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#include "time_experiment.hh"

const int P = 4;                              // number of threads and number of leafs
std::vector<std::atomic_int> sync_levels(P);  // counter shared between 2 threads for one index
std::atomic_int counter = 0;                  // shared counter

// Combining Tree Barrier
void f(int rank) {
    if (P < 3) {
        counter++;
        return;
    }

    int level = 0;
    int index = rank;
    int partner = 0;
    int sync_index = 0;  // sync to this index

    if (index % 2 == 0) {
        partner = index + 1;
        sync_index = index / 2;
    } else {
        partner = index - 1;
        sync_index = (index - 1) / 2;
    }

    sync_levels[index]++;  // each thread initializes its own index
    // wait for partner to initialize its index
    while (sync_levels[partner] == 0) {
    }

    



}

void fill_atomic_vector(std::vector<std::atomic_int>& vec, int value) {
    for (auto& v : vec) {
        v = value;
    }
}

int main() {
    std::vector<std::thread> threads;
    fill_atomic_vector(sync_levels, 0);

    for (int rank = 0; rank < P; ++rank) {
        threads.push_back(std::thread{f, rank});
    }
    for (int rank = 0; rank < P; ++rank) {
        threads[rank].join();
    }
    std::cout << "Counter = " << counter << std::endl;
}
