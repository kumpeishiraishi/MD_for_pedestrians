cmake_minimum_required(VERSION 3.10)

project(test CXX)

add_executable(a.out main.cpp)
target_compile_features(a.out PRIVATE cxx_std_17)