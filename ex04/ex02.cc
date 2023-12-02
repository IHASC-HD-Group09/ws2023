#include <future>
#include <iostream>
#include <numeric>
#include <vector>

const int N = 512;

int STEP_SIZE = std::max(1, static_cast<int>(N / (std::thread::hardware_concurrency())));

double parallelScalarProduct(std::vector<double> vec1, std::vector<double> vec2) {
    std::vector<std::future<double>> futures;

    for (int i = 0; i < vec1.size(); i += STEP_SIZE) {
        auto start_vec_1 = std::begin(vec1) + i;
        auto end_vec_1 = start_vec_1 + STEP_SIZE;
        end_vec_1 = std::min(end_vec_1, std::end(vec1));
        auto start_vec_2 = std::begin(vec2) + i;
        auto end_vec_2 = start_vec_2 + STEP_SIZE;
        end_vec_2 = std::min(end_vec_2, std::end(vec2));
        std::vector<double> sliced_vec1(start_vec_1, end_vec_1);
        std::vector<double> sliced_vec2(start_vec_2, end_vec_2);

        auto future = std::async(
            std::launch::async,
            [](std::vector<double> a, std::vector<double> b) {
                return std::inner_product(std::begin(a), std::end(a), std::begin(b), 0.0);
            },
            sliced_vec1, sliced_vec2);

        futures.push_back(std::move(future));
    }

    double result = 0;
    for (auto& future : futures) {
        result += future.get();
    }
    return result;
}

void linear_fill(std::vector<double>& vec) {
    for (int i = 0; i < vec.size(); ++i) {
        vec[i] = i;
    }
}

int main() {
    std::vector<double> test(N, 2.0);
    std::vector<double> test2(N, 2.0);
    linear_fill(test);

    auto result = parallelScalarProduct(test, test2);
    std::cout << result << std::endl;
}
