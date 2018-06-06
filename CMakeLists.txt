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
add_library(tapp-mesh
    "${CMAKE_CURRENT_SOURCE_DIR}/src/tapp-mesh.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/tapp-weighted-mesh.cpp"
    )
target_include_directories(tapp-mesh PUBLIC src)

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
        tapp-mesh_test
        tests/main.cpp
        tests/tapp-mesh_test.cpp
        tests/tapp-weighted-mesh_test.cpp
        )
    add_test(NAME tapp-mesh_test COMMAND tapp-mesh_test)
    target_link_libraries(tapp-mesh_test doctest tapp-mesh)
else()
    message("-- Testing library not found. Ignoring tests...")
endif()
