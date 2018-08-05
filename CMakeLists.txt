cmake_minimum_required(VERSION 3.9)
project(line_segment_cluster)

set(CMAKE_CXX_STANDARD 11)

add_executable(line_segment_cluster
        lines.txt
        main.cpp clustering.h)
