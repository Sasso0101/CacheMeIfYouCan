cmake_minimum_required(VERSION 3.10)

# Set the project name
project(BFS)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -O3 -march=native -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-fopt-info-omp-vec-optimized-missed=optinfo.txt")

include_directories(include)

# Add sources
add_executable(BFS src/bfs.cpp src/benchmark.cpp src/main.cpp)

include(FetchContent)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(BFS PRIVATE OpenMP::OpenMP_CXX)
endif()

# Add ParlayLib
FetchContent_Declare(
  parlaylib
  GIT_REPOSITORY https://github.com/cmuparlay/parlaylib.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(parlaylib)
target_link_libraries(BFS PRIVATE parlay)

# Make a list of datasets and download them if they don't exist
set(DATASETS
    test.txt
    powergrid.txt
    random_1k_5k.txt
    twitch.txt
    epinions.txt
    unconnected.txt
    pokec.txt
    roadnet-ca.txt
)

foreach(dataset ${DATASETS})
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/datasets/${dataset})
        file(DOWNLOAD http://salvatoreandaloro.altervista.org/datasets/${dataset} ${PROJECT_SOURCE_DIR}/datasets/${dataset})
    endif()
endforeach()

# Datasets
add_custom_target(test
    COMMAND BFS datasets/test.txt 16 32 $(SOURCE)
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(powergrid
    COMMAND BFS datasets/powergrid.txt 4942 13190 $(SOURCE)
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(random
    COMMAND BFS datasets/random_1k_5k.txt 1000 9940 $(SOURCE)
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(twitch
    COMMAND BFS datasets/twitch.txt 168114 13595114 $(SOURCE)
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(epinions
    COMMAND BFS datasets/epinions.txt 75888 811480 $(SOURCE)
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(unconnected
    COMMAND BFS datasets/unconnected.txt 4 2 $(SOURCE)
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(pokec
    COMMAND BFS datasets/pokec.txt 1632804 44603928 $(SOURCE)
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(roadnet-ca
    COMMAND BFS datasets/roadnet-ca.txt 1971281 5533214 $(SOURCE)
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)