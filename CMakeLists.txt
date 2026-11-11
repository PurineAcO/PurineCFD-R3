cmake_minimum_required(VERSION 3.20)
project(PurineCFD-R3 VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PURINE_WERROR "Treat warnings as errors" ON)

# 算例输出根目录, 程序在此下按 3D-<n>-<tag> 自动编号
set(PURINE_DATA_DIR "C:/Users/Administrator/Desktop/PurineCFD-Plugin/data"
    CACHE PATH "Root directory for solver output cases")

set(PURINE_WARNINGS
    -Wall -Wextra -Wpedantic
    -Wno-unused-parameter -Wno-unused-variable
    -fno-math-errno)
if(PURINE_WERROR)
  list(APPEND PURINE_WARNINGS -Werror)
endif()

file(GLOB_RECURSE PURINE_SRC CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp")

add_executable(PurineCFD-R3 src/main.cpp ${PURINE_SRC})
target_compile_options(PurineCFD-R3 PRIVATE ${PURINE_WARNINGS} -fopenmp)
target_include_directories(PurineCFD-R3 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_compile_definitions(PurineCFD-R3 PRIVATE PURINE_DATA_DIR="${PURINE_DATA_DIR}")
target_link_options(PurineCFD-R3 PRIVATE -fopenmp)
