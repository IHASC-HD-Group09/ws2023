cmake_minimum_required(VERSION 3.16.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(ex01 VERSION 1.0)
add_compile_options(-Wall -O1 -march=native)

include_directories(vcl)


add_executable(midpoint midpoint.cc)


