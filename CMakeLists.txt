project(grid)

cmake_minimum_required(VERSION 3.1)

# Set project properties.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set up compiler options based on platform.
if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()
if (MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

# Include library target.
add_library(grid
    "${CMAKE_CURRENT_SOURCE_DIR}/src/grid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/grid_file.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/xml_reader.cpp"
    )
target_include_directories(grid PUBLIC src)

# Include tests if testing library is present.
unset(DOCTEST_LIBRARY_PATH CACHE)
find_path(
    DOCTEST_LIBRARY_PATH
    "doctest.h"
    PATHS "${CMAKE_CURRENT_SOURCE_DIR}"
    PATH_SUFFIX "ext/doctest/doctest"
    )
if (DOCTEST_LIBRARY_PATH)
    add_library(doctest INTERFACE)
    target_include_directories(doctest INTERFACE ext/doctest/doctest)
    enable_testing()
    add_executable(
        grid_test
        tests/main.cpp
        tests/grid_test.cpp
        tests/grid_file_test.cpp
        tests/xml_reader_test.cpp
        )
    add_test(NAME grid_test COMMAND grid_test)
    target_link_libraries(grid_test doctest grid)

    add_executable(
        test-utils_test
        tests/utils/utils_test.cpp
        )
    add_test(NAME test-utils_test COMMAND test-utils_test)
    target_link_libraries(test-utils_test doctest)
else()
    message("-- Testing library not found. Ignoring tests...")
endif()
