CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
PROJECT(ball_atomtyping_validation CXX)
SET(CMAKE_CXX_STANDARD 14)

# BALL-specific setup
SET(BALL_DIR "/mnt/c/Users/Sam/source/repos/ball/build/cmake")
FIND_PACKAGE(BALL REQUIRED)

# Add required BALL include directories
INCLUDE_DIRECTORIES(${BALL_INCLUDE_DIRS})

# Add required BALL compiler flags and definitions
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BALL_CXX_FLAGS}")
ADD_DEFINITIONS("${BALL_COMPILE_DEFINITIONS}")

# Tell cmake to generate an executable called example that depends
# on the file example.cpp
ADD_EXECUTABLE(ball_atomtyping_validation main.cpp)

# Link example against BALL
TARGET_LINK_LIBRARIES(ball_atomtyping_validation BALL)