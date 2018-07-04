project(centroid)

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
add_library(centroidlib
    "${CMAKE_CURRENT_SOURCE_DIR}/src/centroid.cpp"
    )
target_include_directories(centroidlib PUBLIC src )
target_link_libraries(centroidlib ${CMAKE_THREAD_LIBS_INIT})

# Include executable target.
#add_executable(centroid
    #"${CMAKE_CURRENT_SOURCE_DIR}/src/centroid_command.cpp"
    #)
#target_link_libraries(centroid stdc++ stdc++fs centroidlib ${CMAKE_THREAD_LIBS_INIT})

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
        centroid_test
        tests/main.cpp
        tests/centroid_test.cpp
        )
    add_test(NAME centroid_test COMMAND centroid_test)
    target_link_libraries(centroid_test stdc++ doctest centroidlib)
else()
    message("-- Testing library not found. Ignoring tests...")
endif()
