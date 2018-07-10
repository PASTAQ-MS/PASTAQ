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

find_package(Threads REQUIRED)

# Include library target.
add_library(gridlib
    "${CMAKE_CURRENT_SOURCE_DIR}/src/grid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/grid_files.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/grid_runners.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/xml_reader.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/base64.cpp"
    )
target_include_directories(gridlib PUBLIC src )
target_link_libraries(gridlib ${CMAKE_THREAD_LIBS_INIT})

# Include executable target.
add_executable(grid
    "${CMAKE_CURRENT_SOURCE_DIR}/src/grid_command.cpp"
    )
target_link_libraries(grid stdc++ stdc++fs gridlib ${CMAKE_THREAD_LIBS_INIT})

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
        tests/grid_runners_test.cpp
        tests/grid_files_test.cpp
        tests/xml_reader_test.cpp
        )
    add_test(NAME grid_test COMMAND grid_test)
    target_link_libraries(grid_test stdc++ doctest gridlib)

    add_executable(
        test-utils_test
        tests/utils/utils_test.cpp
        )
    add_test(NAME test-utils_test COMMAND test-utils_test)
    target_link_libraries(test-utils_test stdc++ doctest)
else()
    message("-- Testing library not found. Ignoring tests...")
endif()
