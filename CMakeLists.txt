cmake_minimum_required(VERSION 3.6)
project(KPZ_Compact)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(KPZ_Compact ${SOURCE_FILES})