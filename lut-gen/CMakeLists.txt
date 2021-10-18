cmake_minimum_required(VERSION 3.16)
project(lut-gen)

set(CMAKE_CXX_STANDARD 14)

include_directories(ext)

add_executable(lut-Emu-MC Emu_MC.cpp)
add_executable(lut-Eavg-MC Eavg_MC.cpp)
add_executable(lut-Emu-IS Emu_IS.cpp)
add_executable(lut-Eavg-IS Eavg_IS.cpp)

if(MSVC)
	target_compile_definitions(lut-Emu-MC PUBLIC _USE_MATH_DEFINES)
	target_compile_definitions(lut-Eavg-MC PUBLIC _USE_MATH_DEFINES)
	target_compile_definitions(lut-Emu-IS PUBLIC _USE_MATH_DEFINES)
	target_compile_definitions(lut-Eavg-IS PUBLIC _USE_MATH_DEFINES)
endif()
