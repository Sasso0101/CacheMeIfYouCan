cmake_minimum_required(VERSION 3.10)

# Set the project name
project(BFS)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_COMPILER "/apps/gcc-9.1.0/local/bin/g++-9.1.0")
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(LIKWID_DIR "/apps/likwid-4.3.4/lib")

include_directories(include)

# Add sources
add_executable(BFS src/complete.cpp src/reference.cpp src/main.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(BFS PRIVATE OpenMP::OpenMP_CXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp") #=libomp
endif()

find_package(LIKWID)
if(LIKWID_FOUND)
    message(STATUS "Found LIKWID: ${LIKWID_LIBRARIES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -DLIKWID_PERFMON")
    target_include_directories(BFS PRIVATE ${LIKWID_INCLUDES})
    target_link_libraries(BFS PRIVATE likwid)
    target_link_options(BFS PRIVATE -L${LIKWID_DIR})
endif()

include(FetchContent)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)

target_link_libraries(BFS PRIVATE nlohmann_json::nlohmann_json)

# Datasets
add_custom_target(speedcode
    COMMAND BFS
    DEPENDS BFS
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)