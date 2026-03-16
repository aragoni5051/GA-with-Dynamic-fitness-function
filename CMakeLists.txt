cmake_minimum_required(VERSION 3.16)
project(GA_WITH_DYNAMIC_FITNESS_FUNCTION LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SRC_FILES CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/src/*.cpp")

add_executable(ga ${SRC_FILES})
target_include_directories(ga PRIVATE "${CMAKE_SOURCE_DIR}/include")