cmake_minimum_required(VERSION 3.9)

SET(CMAKE_CXX_FLAGS_RELEASE "-std=c++11")
SET(CMAKE_CXX_FLAGS_DEBUG "-g")
SET(CMAKE_BUILD_TYPE RELEASE)

ADD_EXECUTABLE(laser_line_cluster main.cpp)
