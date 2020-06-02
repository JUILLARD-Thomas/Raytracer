PROJECT(RayTrace)
cmake_minimum_required(VERSION 3.1)



set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11
set(CMAKE_CXX_FLAGS "-ljsoncpp")

MESSAGE(STATUS "C++11 enabled. ")



SET(SRCS
  Vector.cpp
  Raytracer.cpp
  )
  
add_executable(lray ${SRCS})
target_link_libraries(lray "-ljsoncpp -lglut -lGLU -lGL" )