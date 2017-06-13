cmake_minimum_required(VERSION 2.6)
project(g2o_trajopt)

add_executable(g2o_trajopt main.cpp)

install(TARGETS g2o_trajopt RUNTIME DESTINATION bin)
