cmake_minimum_required(VERSION 3.4)
project(vkr_sim)

## Compile as C++11
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(submodules/kdl_parser/kdl_parser)

# Debug option
ADD_DEFINITIONS(-DDEBUG)

find_package(Eigen3 REQUIRED)
find_package(orocos_kdl REQUIRED)
find_package(Boost COMPONENTS system thread REQUIRED)

message("message ${orocos_kdl_INCLUDE_DIRS}")

include_directories(
    ${Eigen3_INCLUDE_DIRS}
    ${orocos_kdl_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    submodules/kdl_parser/kdl_parser/include
)

add_executable(vkr_test src/test.cpp)
target_link_libraries(vkr_test ${orocos_kdl_LIBRARIES} kdl_parser Eigen3::Eigen)