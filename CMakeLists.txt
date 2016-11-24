cmake_minimum_required(VERSION 3.6)
project(KPZ_Compact)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O0")

set(SOURCE_FILES main.cpp)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "C:/Users/nisha/ClionProjects/KPZ_Compact")
add_executable(KPZ_Compact ${SOURCE_FILES})