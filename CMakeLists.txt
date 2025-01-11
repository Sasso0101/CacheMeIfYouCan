cmake_minimum_required(VERSION 3.10)

# Set the project name
project(BFS)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

include_directories(include)

# Add sources
add_executable(BFS src/bfs.cpp src/benchmark.cpp src/bfs_cached.cpp src/bfs_hybrid_bitmap.cpp src/test_openmp.cpp src/main.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(BFS PRIVATE OpenMP::OpenMP_CXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
endif()

find_package(likwid)
if(LIKWID_FOUND)
    message(STATUS "Found LIKWID: ${LIKWID_LIBRARIES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -DLIKWID_PERFMON")
    target_link_libraries(BFS PRIVATE likwid::likwid)
    target_link_options(BFS PRIVATE -lm)
endif()

# find_package(absl REQUIRED)
# target_link_libraries(BFS PRIVATE absl::flat_hash_map)

# Add ParlayLib
# include(FetchContent)
# FetchContent_Declare(
#   parlaylib
#   GIT_REPOSITORY https://github.com/cmuparlay/parlaylib.git
#   GIT_TAG        master
# )
# FetchContent_GetProperties(parlaylib)
# if(NOT parlaylib_POPULATED)
#   FetchContent_Populate(parlaylib)  
#   add_subdirectory(${parlaylib_SOURCE_DIR} EXCLUDE_FROM_ALL)
# endif()
# target_link_libraries(BFS PRIVATE parlay)

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