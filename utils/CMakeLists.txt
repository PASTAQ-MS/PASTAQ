project(utils)

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
if (NOT TARGET testutilslib)
    add_library(testutilslib INTERFACE)
    target_include_directories(testutilslib INTERFACE src/testutils)
endif()

# Include tests if testing library is present.
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} OR (${UTILS_ADD_TESTS}))
    unset(DOCTEST_LIBRARY_PATH CACHE)
    find_path(
        DOCTEST_LIBRARY_PATH
        "doctest.h"
        PATHS "${CMAKE_CURRENT_SOURCE_DIR}"
        PATH_SUFFIX "${CMAKE_CURRENT_SOURCE_DIR}/../ext/doctest/doctest"
        )
    if (DOCTEST_LIBRARY_PATH)
        if(NOT TARGET doctest)
            add_library(doctest INTERFACE)
            target_include_directories(doctest
                INTERFACE
                "${CMAKE_CURRENT_SOURCE_DIR}/../ext/doctest/doctest")
        endif()
        enable_testing()
        if (NOT TARGET utils_test)
            add_executable(
                utils_test
                tests/main.cpp
                tests/mock_stream_test.cpp
                )
            add_test(NAME utils_test COMMAND utils_test)
            target_link_libraries(utils_test stdc++ doctest testutilslib)
        endif()
    else()
        message("-- [${PROJECT_NAME}] Testing library not found. Ignoring tests...")
    endif()
endif()
