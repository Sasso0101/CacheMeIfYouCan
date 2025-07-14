cmake_minimum_required(VERSION 3.20)

# Set the project name
project(BFS)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -march=native")

include_directories(include)

# Add sources
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/implementations/*.cpp")
add_executable(BFS ${SOURCES})

add_subdirectory(distributed_mmio)
target_link_libraries(BFS PRIVATE distributed_mmio)

# Add OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp") #=libomp
    target_link_libraries(BFS PRIVATE OpenMP::OpenMP_CXX)
endif()

include(FetchContent)

# Add nlohmann/json
include(FetchContent)
FetchContent_Declare(json 
    URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP FALSE)
FetchContent_MakeAvailable(json)
target_link_libraries(BFS PRIVATE nlohmann_json::nlohmann_json)

# Add Debug Target
add_custom_target(debug
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_BUILD_TOOL} -j
    COMMAND gdb --args $<TARGET_FILE:BFS>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Building and running BFS in debug mode with gdb"
)