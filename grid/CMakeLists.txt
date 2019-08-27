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

# Add utils dependency.
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../utils" ${CMAKE_CURRENT_BINARY_DIR}/utils)

# Include library target.
if (NOT TARGET gridlib)
    add_library(gridlib
        "${CMAKE_CURRENT_SOURCE_DIR}/src/grid/grid.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/grid/grid_serialize.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/grid/raw_data_serialize.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/grid/grid_files.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/grid/grid_runners.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/grid/xml_reader.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/grid/raw_data.cpp"
        )
    target_include_directories(gridlib PUBLIC src)
    target_link_libraries(gridlib ${CMAKE_THREAD_LIBS_INIT} utilslib)
endif()

# Include executable target.
if (NOT TARGET grid)
    add_executable(grid
        "${CMAKE_CURRENT_SOURCE_DIR}/src/grid/grid_command.cpp"
        )
    target_link_libraries(grid stdc++ stdc++fs gridlib ${CMAKE_THREAD_LIBS_INIT})
endif()

# Include tests if testing library is present.
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} OR (${GRID_ADD_TESTS}))
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
        if (NOT TARGET grid_test)
            add_executable(
                grid_test
                tests/main.cpp
                tests/grid_test.cpp
                tests/grid_runners_test.cpp
                tests/grid_files_test.cpp
                tests/xml_reader_test.cpp
                )
            add_test(NAME grid_test COMMAND grid_test)
            target_link_libraries(grid_test stdc++ doctest gridlib testutilslib)
        endif()
    else()
        message("-- [${PROJECT_NAME}] Testing library not found. Ignoring tests...")
    endif()
endif()
