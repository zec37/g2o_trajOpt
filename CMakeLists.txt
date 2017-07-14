cmake_minimum_required(VERSION 2.6)
project(g2o_trajopt)
set ( CMAKE_BUILD_TYPE "Debug" )
set( CMAKE_CXX_FLAGS "-O3" )

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

# g2o
find_package(G2O REQUIRED)
include_directories(${G2O_INCLUDE_DIR})

# cholmod
find_package(Cholmod REQUIRED)
include_directories(${CHOLMOD_INCLUDE_DIR})

# eigen
include_directories("/usr/local/include/eigen3")
include_directories("/usr/local/include/eigen3/Eigen")

#find_package(OpenGL)
#find_package(GLUT)

add_executable( g2o_trajopt traj_match.cpp)
target_link_libraries( ${PROJECT_NAME} 
    g2o_core g2o_stuff g2o_types_slam2d
    ${CHOLMOD_LIBRARIES}
    ${EIGEN_LIBRARIES}
    #${GLUT_LIBRARIES}
)

install(TARGETS g2o_trajopt RUNTIME DESTINATION bin)
