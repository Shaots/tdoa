cmake_minimum_required(VERSION 3.24)
project(tdoa)

set(CMAKE_CXX_STANDARD 17)

add_executable(tdoa main.cpp point.cpp Solver.cpp uniTest.cpp SupportFunc.cpp)
