cmake_minimum_required(VERSION 3.0)
project(SC4PS-01)
set(CMAKE_CXX_STANDARD 11)

# include boost program options
set(BOOSTROOT /usr/lib64)
find_package(Boost 1.70 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Path: 01-hello-world/CMakeLists.txt
add_executable(hello-world main.cpp)
target_link_libraries(hello-world ${Boost_LIBRARIES})
