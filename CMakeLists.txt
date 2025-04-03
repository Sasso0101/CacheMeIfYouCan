cmake_minimum_required(VERSION 3.15)

# Set the project name
project(BFS)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -march=native")

include_directories(include)
include_directories(/usr/lib/gcc/x86_64-redhat-linux/4.8.2/include/)

# Add sources
file(GLOB_RECURSE SOURCES "src/*.cpp" "src/implementations/*.cpp")
add_executable(BFS ${SOURCES})

# Add OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(BFS PRIVATE OpenMP::OpenMP_CXX)

include(FetchContent)

# Add nlohmann/json
include(FetchContent)
FetchContent_Declare(json 
    URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
)
FetchContent_MakeAvailable(json)
target_link_libraries(BFS PRIVATE nlohmann_json::nlohmann_json)

FetchContent_Declare(
    cli11_proj
    QUIET
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(cli11_proj)
target_link_libraries(BFS PRIVATE CLI11::CLI11)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(tests)
endif()