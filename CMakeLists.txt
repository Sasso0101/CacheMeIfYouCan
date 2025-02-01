cmake_minimum_required(VERSION 3.10)

# Set the project name
project(BFS)

set(DBG_FRONTIER_SIZE OFF CACHE BOOL "Value for DBG_FRONTIER_SIZE")
if(DBG_FRONTIER_SIZE)
    message("DBG_FRONTIER_SIZE is ON")
else()
    message("DBG_FRONTIER_SIZE is OFF")
endif()

set(DBG_THREAD_BALANCE OFF CACHE BOOL "Value for DBG_THREAD_BALANCE")
if(DBG_THREAD_BALANCE)
    message("DBG_THREAD_BALANCE is ON")
else()
    message("DBG_THREAD_BALANCE is OFF")
endif()

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_COMPILER "/apps/gcc-9.1.0/local/bin/g++-9.1.0")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -march=native")
# set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(LIKWID_DIR "/apps/likwid-4.3.4/lib")

include_directories(include)
include_directories(deps/json/include)

# Add sources
add_executable(BFS src/complete.cpp src/reference.cpp src/main.cpp)
add_executable(BFS_CACHE src/complete.cpp src/reference.cpp src/main.cpp)
add_executable(BFS_REF src/complete.cpp src/reference.cpp src/main_reference.cpp)
add_executable(DEG_DISTR src/utils/test_degree.cpp)

if(DBG_THREAD_BALANCE)
    message("Setting DBG_THREAD_BALANCE")
    target_compile_definitions(BFS PRIVATE DBG_THREAD_BALANCE=1)
endif()

if(DBG_FRONTIER_SIZE)
    message("Setting DBG_FRONTIER_SIZE")
    target_compile_definitions(BFS_REF PRIVATE DBG_FRONTIER_SIZE=1)
endif()

target_compile_definitions(BFS_CACHE PRIVATE DBG_CACHE=1)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp") #=libomp
    target_link_libraries(BFS PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(BFS_CACHE PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(BFS_REF PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(DEG_DISTR PRIVATE OpenMP::OpenMP_CXX)
endif()

# find_package(LIKWID)
# if(LIKWID_FOUND)
#     message(STATUS "Found LIKWID: ${LIKWID_LIBRARIES}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -DLIKWID_PERFMON")
#     target_include_directories(BFS PRIVATE ${LIKWID_INCLUDES})
#     target_include_directories(BFS_REF PRIVATE ${LIKWID_INCLUDES})
#     target_link_libraries(BFS PRIVATE likwid)
#     target_link_options(BFS PRIVATE -L${LIKWID_DIR})
#     target_link_options(BFS_REF PRIVATE -L${LIKWID_DIR})
# endif()

# include(FetchContent)
# set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
# FetchContent_Declare(
#     json 
#     URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
#     SOURCE_SUBDIR single_include
# )
# set(FETCHCONTENT_FULLY_DISCONNECTED ON)
# set(FETCHCONTENT_SOURCE_DIR ${CMAKE_BINARY_DIR}/_deps)
# set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR}/_deps)

# # FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
# FetchContent_MakeAvailable(json)
# target_include_directories(BFS PRIVATE ${json_SOURCE_DIR}/single_include)

# target_link_libraries(BFS PRIVATE json)

# target_link_libraries(BFS PRIVATE nlohmann_json::nlohmann_json)
# target_link_libraries(BFS_REF PRIVATE nlohmann_json::nlohmann_json)
# target_link_libraries(DEG_DISTR PRIVATE nlohmann_json::nlohmann_json)