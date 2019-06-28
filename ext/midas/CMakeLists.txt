cmake_minimum_required(VERSION 3.1)

project(MIDAs)

add_library(midas "${CMAKE_CURRENT_SOURCE_DIR}/src/MIDAs.cpp")
target_include_directories(midas PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
