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

# Add utils dependency.
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../utils" ${CMAKE_CURRENT_BINARY_DIR}/utils)

# Add grid dependency
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../grid" ${CMAKE_CURRENT_BINARY_DIR}/grid)

# Include library target.
if (NOT TARGET centroidlib)
    add_library(centroidlib
        "${CMAKE_CURRENT_SOURCE_DIR}/src/centroid/centroid.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/centroid/centroid_files.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/centroid/centroid_serialize.cpp"
        )
    target_include_directories(centroidlib PUBLIC src)
    target_link_libraries(centroidlib gridlib ${CMAKE_THREAD_LIBS_INIT})
endif()

# Include executable target.
if (NOT TARGET centroid)
    add_executable(centroid
        "${CMAKE_CURRENT_SOURCE_DIR}/src/centroid/centroid_command.cpp"
        )
    target_link_libraries(centroid stdc++ stdc++fs centroidlib ${CMAKE_THREAD_LIBS_INIT})
endif()

# Include tests if testing library is present.
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} OR (${CENTROID_ADD_TESTS}))
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
        if (NOT TARGET centroid_test)
            add_executable(
                centroid_test
                tests/main.cpp
                tests/centroid_test.cpp
                tests/centroid_files_test.cpp
                )
            add_test(NAME centroid_test COMMAND centroid_test)
            target_link_libraries(centroid_test stdc++ doctest testutilslib centroidlib)
        endif()
    else()
        message("-- [${PROJECT_NAME}] Testing library not found. Ignoring tests...")
    endif()
endif()
