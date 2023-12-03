#include <atomic>
#include <cmath>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#include "time_experiment.hh"

const int P = 512*2*2*2;  // number of threads and number of leafs
std::vector<std::atomic_int> sync_counter(P);
std::vector<std::atomic_bool> sleep_flag(P);
std::atomic_int counter = 0;  // shared counter

std::atomic<double> time_barrier = 0.0;

void wake_up(int rank) {
    int partner{};
    for (int k = std::log2(P) - 1; k >= 0; k--) {
        partner = std::pow(2, k) + rank;

        if (partner < P) {
            sleep_flag[partner] = false;
        }
    }
}

void lock(int rank) {
    int partner{};
    for (int k = 0; k < std::log2(P); k++) {
        partner = std::pow(2, k) + rank;

        if ((rank & (1 << k)) >> k) {
            sleep_flag[rank] = true;
            break;
        }

        if (partner < P) {
            while (!sleep_flag[partner]) {}  // wait for sleep of partner
            sync_counter[rank] += sync_counter[partner];
        }
    }
    if (!sleep_flag[rank]) {
        counter += sync_counter[rank];
    }
}

// Combining Tree Barrier
void f(int rank) {
    if (P == 1) {
        counter++;
        return;
    }

    sync_counter[rank] = 1;  // initialize own value
    auto start = std::chrono::high_resolution_clock::now();
    lock(rank);  // start lock from leafs
    while (sleep_flag[rank]) {}
    wake_up(rank);  // start wake up from root
    auto end = std::chrono::high_resolution_clock::now();
    time_barrier = std::max(time_barrier.load(), static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>((end - start)).count()));
    
}

void fill_atomic_vector(std::vector<std::atomic_int>& vec, int value) {
    for (auto& v : vec) {
        v = value;
    }
}

int main() {
    std::vector<std::thread> threads;
    fill_atomic_vector(sync_counter, 0);

    for (int rank = 0; rank < P; ++rank) {
        threads.push_back(std::thread{f, rank});
    }
    for (int rank = 0; rank < P; ++rank) {
        threads[rank].join();
    }
    std::cout << "Counter = " << counter << std::endl;
    std::cout << "Time = " << time_barrier << std::endl;
}
