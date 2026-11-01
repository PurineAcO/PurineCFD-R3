cmake_minimum_required(VERSION 3.20)
project(PurineCFD-R3 VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(PURINE_WERROR "Treat warnings as errors" ON)
option(PURINE_NATIVE "Optimize for the build machine" OFF)
option(PURINE_SANITIZE "Enable AddressSanitizer and UBSan" OFF)

if(MSVC)
  set(PURINE_WARNINGS /W4 /utf-8 /permissive- /wd4100)
  set(PURINE_OPENMP /openmp)
  if(PURINE_WERROR)
    list(APPEND PURINE_WARNINGS /WX)
  endif()
else()
  set(PURINE_WARNINGS
      -Wall -Wextra -Wpedantic
      -Wno-unused-parameter -Wno-unused-variable
      -fno-math-errno)
  set(PURINE_OPENMP -fopenmp)
  if(PURINE_WERROR)
    list(APPEND PURINE_WARNINGS -Werror)
  endif()
endif()

function(purine_add_executable name)
  add_executable(${name} ${ARGN})
  target_compile_features(${name} PRIVATE cxx_std_17)
  target_compile_options(${name} PRIVATE ${PURINE_WARNINGS} ${PURINE_OPENMP})
  if(NOT MSVC)
    target_link_options(${name} PRIVATE -fopenmp)
  endif()
  target_include_directories(${name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
    target_include_directories(${name} SYSTEM PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/third_party")
  endif()
  if(PURINE_NATIVE AND NOT MSVC)
    target_compile_options(${name} PRIVATE -march=native)
  endif()
  if(PURINE_SANITIZE AND NOT MSVC)
    target_compile_options(${name} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(${name} PRIVATE -fsanitize=address,undefined)
  endif()
endfunction()

purine_add_executable(PurineCFD-R3 src/main.cpp
  src/physic.hpp
  src/grad.hpp
  src/readmesh.hpp)

include(CTest)

if(BUILD_TESTING)
  purine_add_executable(mesh_check tests/mesh_check.cpp)
  add_test(NAME mesh_check COMMAND mesh_check WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()
