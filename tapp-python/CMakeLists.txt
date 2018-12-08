project(tapp-python)

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
#add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../utils" ${CMAKE_CURRENT_BINARY_DIR}/utils)

## Add grid dependency
#add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../grid" ${CMAKE_CURRENT_BINARY_DIR}/grid)

## Add centroid dependency
#add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../centroid" ${CMAKE_CURRENT_BINARY_DIR}/centroid)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../ext/pybind11" ${CMAKE_CURRENT_BINARY_DIR}/pybind11)
pybind11_add_module(example src/example.cpp)
#target_include_directories(example PUBLIC pybind11/include)

## Include library target.
#if (NOT TARGET metamatchlib)
    #add_library(metamatchlib
        #"${CMAKE_CURRENT_SOURCE_DIR}/src/metamatch/metamatch.cpp"
        #"${CMAKE_CURRENT_SOURCE_DIR}/src/metamatch/metamatch_files.cpp"
        #)
    #target_include_directories(metamatchlib PUBLIC src)
    #target_link_libraries(metamatchlib gridlib centroidlib ${CMAKE_THREAD_LIBS_INIT})
#endif()

# Include executable target.
#if (NOT TARGET metamatch)
    #add_executable(metamatch
        #"${CMAKE_CURRENT_SOURCE_DIR}/src/metamatch_command.cpp"
        #)
    #target_link_libraries(metamatch stdc++ stdc++fs metamatchlib ${CMAKE_THREAD_LIBS_INIT} gridlib centroidlib)
#endif()

# Include tests if testing library is present.
#if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} OR (${CENTROID_ADD_TESTS}))
    #unset(DOCTEST_LIBRARY_PATH CACHE)
    #find_path(
        #DOCTEST_LIBRARY_PATH
        #"doctest.h"
        #PATHS "${CMAKE_CURRENT_SOURCE_DIR}"
        #PATH_SUFFIX "${CMAKE_CURRENT_SOURCE_DIR}/../ext/doctest/doctest"
        #)
    #if (DOCTEST_LIBRARY_PATH)
        #if(NOT TARGET doctest)
            #add_library(doctest INTERFACE)
            #target_include_directories(doctest
                #INTERFACE
                #"${CMAKE_CURRENT_SOURCE_DIR}/../ext/doctest/doctest")
        #endif()
        #enable_testing()
        #if (NOT TARGET metamatch_test)
            #add_executable(
                #metamatch_test
                #tests/main.cpp
                #tests/metamatch_test.cpp
                ##tests/metamatch_runners_test.cpp
                #)
            #add_test(NAME metamatch_test COMMAND metamatch_test)
            #target_link_libraries(metamatch_test stdc++ doctest testutilslib metamatchlib)
        #endif()
    #else()
        #message("-- [${PROJECT_NAME}] Testing library not found. Ignoring tests...")
    #endif()
#endif()
