project(tapp-grid)

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
add_library(tapp-grid
    "${CMAKE_CURRENT_SOURCE_DIR}/src/tapp-grid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/tapp-datfile.cpp"
    )
target_include_directories(tapp-grid PUBLIC src)

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
        tapp-grid_test
        tests/main.cpp
        tests/tapp-grid_test.cpp
        tests/tapp-datfile_test.cpp
        )
    add_test(NAME tapp-grid_test COMMAND tapp-grid_test)
    target_link_libraries(tapp-grid_test doctest tapp-grid)

    add_executable(
        test-utils_test
        tests/utils/utils_test.cpp
        )
    add_test(NAME test-utils_test COMMAND test-utils_test)
    target_link_libraries(test-utils_test doctest)
else()
    message("-- Testing library not found. Ignoring tests...")
endif()
