#include <iostream>
#include <vector>
#include <thread>

const int P = std::thread::hardware_concurrency(); // number of threads
const int N=25600;            // problem size
std::vector<double> x(N,2.0); // first vector
std::vector<double> y(N,2.0); // second vector
double alpha = 5.0; // scalar alpha

void daxpy (int rank)
{
  for (int i=(N*rank)/P; i<(N*(rank+1))/P; ++i) 
    y[i] = x[i]*alpha + y[i];
}

int main ()
{
  std::vector<std::thread> threads;
  for (int rank=0; rank<P; ++rank)
    threads.push_back(std::thread{daxpy,rank});
  for (int rank=0; rank<P; ++rank)
    threads[rank].join();
  std::cout << "y after daxpy: ";
  for (double i: y)
    std::cout << i << " ";
}
