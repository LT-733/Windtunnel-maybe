cmake_minimum_required(VERSION 3.31)
project(windtunnel_maybe)

set(CMAKE_CXX_STANDARD 20)

find_package(SDL2 REQUIRED)
find_package(OpenMP REQUIRED)

include_directories(${SDL2_INCLUDE_DIRS})

add_executable(windtunnel_maybe src/main3.cpp)

add_compile_definitions(rsize_t=size_t)

target_link_libraries(windtunnel_maybe PRIVATE
        ${SDL2_LIBRARIES}
        OpenMP::OpenMP_CXX
)

target_include_directories(windtunnel_maybe PRIVATE include)