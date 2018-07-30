project(tapp)

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

# Configure the tests.
option(CENTROID_ADD_TESTS "run centroid tests" ON)
option(GRID_ADD_TESTS "run grid tests" ON)
option(UTILS_ADD_TESTS "run utils tests" ON)

# Include library target.
add_subdirectory(${CMAKE_SOURCE_DIR}/grid)
add_subdirectory(${CMAKE_SOURCE_DIR}/centroid)
add_subdirectory(${CMAKE_SOURCE_DIR}/warp2d)

# TODO(alex): Ensure library can be linked against using:
#     target_link_libraries(my_lib ${TAPP_LIBS})
set(TAPP_LIBS "gridlib centroidlib")

unset(DOCTEST_LIBRARY_PATH CACHE)
find_path(
    DOCTEST_LIBRARY_PATH
    "doctest.h"
    PATHS "${CMAKE_CURRENT_SOURCE_DIR}"
    PATH_SUFFIX "${CMAKE_CURRENT_SOURCE_DIR}/ext/doctest/doctest"
    )
if (DOCTEST_LIBRARY_PATH)
    enable_testing()
else()
    message("-- [${PROJECT_NAME}] Testing library not found. Ignoring tests...")
endif()
