include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Set the project source path as an environment variable
set(PROJECT_SOURCE_PATH "${CMAKE_SOURCE_DIR}" CACHE INTERNAL "Project source path")
set(ENV{PROJECT_SOURCE_PATH} "${PROJECT_SOURCE_PATH}")

include_directories(include)

# Enable testing
enable_testing()

# Add the test executable
file(GLOB_RECURSE SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/implementations/*.cpp")
list(FILTER SOURCES EXCLUDE REGEX ".*main\\.cpp$")
add_executable(tests ${SOURCES} run_all.cpp)

target_link_libraries(tests PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(tests PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(tests PRIVATE GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(tests)