project(tapp)

cmake_minimum_required(VERSION 3.1)

# Set project properties.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set up compiler options based on platform.
if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()
if (MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

# Add required dependencies
# Threads
find_package(Threads REQUIRED)
# MIDAs
# add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../ext/midas" ${CMAKE_CURRENT_BINARY_DIR}/midas)
# PyBind
#add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../ext/pybind11" ${CMAKE_CURRENT_BINARY_DIR}/pybind11)
#pybind11_add_module(tapp src/tapp/bindings.cpp)
#target_link_libraries(tapp PRIVATE gridlib centroidlib warp2dlib metamatchlib midas stdc++ stdc++fs)

# Build TAPPs library.
add_library(tapplib
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/grid/grid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/grid/grid_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/grid/raw_data.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/grid/raw_data_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/grid/xml_reader.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/centroid/centroid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/centroid/centroid_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/warp2d/warp2d.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/metamatch/metamatch.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/metamatch/metamatch_files.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/metamatch/metamatch_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/link/link_msms.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/link/link_msms_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/base64.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/interpolation.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/search.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/serialization.cpp"
    )
target_include_directories(tapplib PUBLIC src/lib)
target_link_libraries(tapplib ${CMAKE_THREAD_LIBS_INIT})

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
