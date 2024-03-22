cmake_minimum_required(VERSION 3.1)

project(pastaq VERSION "0.11.0")

# Set project properties.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set up compiler options based on platform.
if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
endif()
if (MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

add_subdirectory("ext/zlib")
add_subdirectory("ext/eigen")

# Build pastaq library.
add_library(pastaqlib
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/centroid/centroid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/centroid/centroid_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/feature_detection/feature_detection.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/feature_detection/feature_detection_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/grid/grid.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/grid/grid_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/link/link.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/link/link_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/metamatch/metamatch.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/metamatch/metamatch_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/protein_inference/protein_inference.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/protein_inference/protein_inference_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/raw_data/raw_data.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/raw_data/raw_data_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/raw_data/xml_reader.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/base64.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/compression.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/interpolation.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/search.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/serialization.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/warp2d/warp2d.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/warp2d/warp2d_serialize.cpp"
    )

target_include_directories(pastaqlib PUBLIC src/lib)
if (WIN32)
    target_link_libraries(pastaqlib ${CMAKE_THREAD_LIBS_INIT} Eigen3::Eigen zlibstatic)
else()
    target_link_libraries(pastaqlib ${CMAKE_THREAD_LIBS_INIT} Eigen3::Eigen zlib)
endif()

# Build the python bindings.
# --------------------------
# PyBind.
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11")
pybind11_add_module(pastaq src/python-bindings/pastaq/bindings.cpp)
target_link_libraries(pastaq PRIVATE pastaqlib)

# Include tests if testing library is present.
# --------------------------------------------
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR} AND (${PASTAQ_ENABLE_TESTS}))
    unset(DOCTEST_LIBRARY_PATH CACHE)
    find_path(
        DOCTEST_LIBRARY_PATH
        "doctest.h"
        PATHS "${CMAKE_CURRENT_SOURCE_DIR}"
        PATH_SUFFIX "${CMAKE_CURRENT_SOURCE_DIR}/ext/doctest/doctest"
        )
    if (DOCTEST_LIBRARY_PATH)
        # Add the testing library.
        add_library(doctest INTERFACE)
        target_include_directories(doctest
            INTERFACE
            "${CMAKE_CURRENT_SOURCE_DIR}/ext/doctest/doctest")
        enable_testing()
        # Add tests.
        add_executable(
            pastaqlib_test
            tests/centroid_test.cpp
            tests/feature_detection_test.cpp
            tests/grid_test.cpp
            tests/main.cpp
            tests/metamatch_test.cpp
            tests/mock_stream_test.cpp
            tests/serialization_test.cpp
            tests/warp2d_test.cpp
            tests/xml_reader_test.cpp
            )
        add_test(NAME pastaqlib_test COMMAND pastaqlib_test)
        target_link_libraries(pastaqlib_test stdc++ doctest pastaqlib)
    else()
        message("-- [${PROJECT_NAME}] Testing library not found. Ignoring tests...")
    endif()
endif()
