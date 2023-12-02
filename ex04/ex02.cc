#include <iostream>
#include <vector>
#include <future>
#include <numeric>

int TASK_SIZE = 4;

double parallelScalarProduct(std::vector<double> vec1, std::vector<double> vec2) {
    std::vector<std::future<double>> futures;

    for (int i = 0; i < vec1.size(); i+= TASK_SIZE)
    {
        std::cout << i << std::endl;
        auto start_vec_1 = std::begin(vec1) + i;
        auto end_vec_1 = start_vec_1 + TASK_SIZE;
        auto start_vec_2 = std::begin(vec2) + i;
        auto end_vec_2 = start_vec_2 + TASK_SIZE;
        std::vector<double> sliced_vec1(end_vec_1 - start_vec_1);
        std::vector<double> sliced_vec2(end_vec_2 - start_vec_2);
        std::copy(start_vec_1, end_vec_1, std::begin(sliced_vec1));
        std::copy(start_vec_2, end_vec_2, std::begin(sliced_vec2));

        auto future = std::async(std::launch::async, [](std::vector<double> a, std::vector<double> b) {
            return std::inner_product(std::begin(a), std::end(b), std::begin(b), 0.0);
        }, sliced_vec1, sliced_vec2);

        futures.push_back(std::move(future));
    }

    double result = 0;

    for (auto& future: futures) {
        result += future.get();
    }

    return result;
}

int main() {
    std::vector<double> test(16, 2.0);
    std::vector<double> test2(16, 2.0);

    auto result = parallelScalarProduct(test, test2);
    std::cout << result << std::endl;
}