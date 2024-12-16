cmake_minimum_required(VERSION 3.10)

# Set the project name
project(BFS)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(BFS src/main.cpp)

# Add ParlayLib
include(FetchContent)
FetchContent_Declare(
  parlaylib
  GIT_REPOSITORY https://github.com/cmuparlay/parlaylib.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(parlaylib)
target_link_libraries(BFS PUBLIC parlay)

# Datasets
add_custom_target(test
    COMMAND BFS datasets/test.txt 11 26
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(powergrid
    COMMAND BFS datasets/USpowergrid_n4941.txt 4942 13190
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(random
    COMMAND BFS datasets/random_1k_5k.txt 1000 9940
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)