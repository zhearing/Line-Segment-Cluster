cmake_minimum_required(VERSION 3.9)
project(laser_line_cluster)

set(CMAKE_CXX_STANDARD 11)

add_executable(laser_line_cluster
        lines.txt
        main.cpp clustering.h)
