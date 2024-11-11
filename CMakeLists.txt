cmake_minimum_required(VERSION 3.10)

project(GasFlow VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(src)

file(GLOB SOURCES "src/*.cpp")

add_executable(GasFlow ${SOURCES})

find_package(Threads REQUIRED)
target_link_libraries(GasFlow Threads::Threads)
