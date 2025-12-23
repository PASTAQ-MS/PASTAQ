cmake_minimum_required(VERSION 3.14)

project(pastaq VERSION "0.11.3")

# Set project properties.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# Set paths to TDF SDK libraries and include directories
set(TDF_SDK_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/ext/third_party/timsdata_sdk/include/c")

if(WIN32)
    # Windows-specific settings
    set(TDF_SDK_LIB "${CMAKE_SOURCE_DIR}/src/python-bindings/pastaq/win64")
    find_library(TIMSDATA_LIBRARY timsdata.dll PATHS "${TDF_SDK_LIB}")
    message(STATUS "Building on Windows")
elseif(UNIX)
    # Linux-specific settings
    set(TDF_SDK_LIB "${CMAKE_SOURCE_DIR}/src/python-bindings/pastaq/linux64")
    find_library(TIMSDATA_LIBRARY libtimsdata.so PATHS "${TDF_SDK_LIB}")
    message(STATUS "Building on Linux")
else()
    message(FATAL_ERROR "Unsupported platform")
endif()


set(SQLITE_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/ext/third_party/CppSQLite")
set(SQLITE_DIR "${CMAKE_SOURCE_DIR}/ext/third_party/sqlite")

# Set up compiler options based on platform.
if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
endif()
if (MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)  # Default to Debug if not specified
endif()

# Add debug-specific flags.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Project is build in debug mode.")
    if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
    endif()
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi /Od")
    endif()
else()
    message(STATUS "Project is build in release mode.")
endif()

add_subdirectory("ext/zlib")
add_subdirectory("ext/eigen")
# add_subdirectory("ext/mstoolkit" "${CMAKE_BINARY_DIR}/mstoolkit_build")
add_subdirectory(
  "ext/mstoolkit"
  "${CMAKE_BINARY_DIR}/mstoolkit_build"
  EXCLUDE_FROM_ALL
)
# add_subdirectory(src/lib/raw_data)

# Build pastaq library.
add_library(pastaqlib STATIC
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
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/raw_data/timsdatacpp.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/raw_data/raw_data.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/raw_data/raw_data_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/raw_data/xml_reader.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/base64.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/compression.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/interpolation.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/search.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/utils/serialization.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/visualization/tims_visualization.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/warp2d/warp2d.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/lib/warp2d/warp2d_serialize.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ext/third_party/CppSQLite/CppSQLite3.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ext/third_party/sqlite/shell.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/ext/third_party/sqlite/sqlite3.c"
    )

# Try to find HDF5 (cross-platform)
find_package(HDF5 REQUIRED COMPONENTS C CXX)
message(STATUS "Found HDF5 include dirs: ${HDF5_INCLUDE_DIRS}")
message(STATUS "Found HDF5 libraries: ${HDF5_LIBRARIES}")

target_include_directories(pastaqlib PUBLIC src/lib
    ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/raw_data/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src/lib/include
    ${HDF5_INCLUDE_DIRS}
)

add_definitions(-DMZP_HDF)

# Add mstoolkit include directory
IF (UNIX)
  # for linux
  add_definitions(-DMZP_HDF -DGCC -DXML_DEV_URANDOM -D_FILE_OFFSET_BITS=64 -D_NO_THERMORAW -D_NOSQLITE)
ENDIF()
target_include_directories(pastaqlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ext/mstoolkit/include)
target_include_directories(pastaqlib PUBLIC "${SQLITE_INCLUDE_DIR}")
target_include_directories(pastaqlib PUBLIC "${SQLITE_DIR}")
target_include_directories(pastaqlib PUBLIC "${TDF_SDK_INCLUDE_DIR}")
target_include_directories(mstoolkit PUBLIC "${HDF5_INCLUDE_DIRS}")
# target_include_directories(pastaqlib PUBLIC "${HDF5_INCLUDE_DIRS}")

target_compile_definitions(mstoolkit PUBLIC MZP_HDF)

if (WIN32)
    target_link_libraries(pastaqlib PRIVATE ${CMAKE_THREAD_LIBS_INIT}
        Eigen3::Eigen 
        zlibstatic
        mstoolkit PUBLIC
        ${HDF5_LIBRARIES}
        ${TIMSDATA_LIBRARY}
    )
#    target_link_libraries(pastaqlib PUBLIC mstoolkit ${HDF5_LIBRARIES})
else()
    target_link_libraries(pastaqlib PRIVATE ${CMAKE_THREAD_LIBS_INIT} ${HDF5_LIBRARIES}
        Eigen3::Eigen 
        zlib
        mstoolkit PUBLIC
        ${HDF5_LIBRARIES}
        ${TIMSDATA_LIBRARY}
    )
#    target_link_libraries(pastaqlib PUBLIC mstoolkit ${HDF5_LIBRARIES})
endif()

add_executable(TimsDataExample "${CMAKE_CURRENT_SOURCE_DIR}/examples/cpp/timstdf.cpp")
add_executable(TimsDataSamples "${CMAKE_CURRENT_SOURCE_DIR}/examples/cpp/timsdataSamples.cpp")
add_executable(TestTimsQuery "${CMAKE_CURRENT_SOURCE_DIR}/examples/cpp/testquery.cpp")

# Link the pastaqlib library to the executable
target_link_libraries(TimsDataExample PRIVATE pastaqlib)
target_link_libraries(TimsDataSamples PRIVATE pastaqlib)
target_link_libraries(TestTimsQuery PRIVATE pastaqlib)



# Build the python bindings.
# --------------------------
# PyBind.
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11")

# Define source files for the `pastaq` module.
set(PASTAQ_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/python-bindings/pastaq/bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/python-bindings/pastaq/timsdata_pybind.cpp
)

# Define include directories.
set(PASTAQ_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/python-bindings/pastaq/include
)

pybind11_add_module(pastaq ${PASTAQ_SOURCES})

# Set include directories for the target.
target_include_directories(pastaq PRIVATE ${PASTAQ_INCLUDE_DIRS})

if(UNIX)
    set_target_properties(pastaq PROPERTIES
        INSTALL_RPATH "$ORIGIN/linux64"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
elseif(WIN32)
    set_target_properties(pastaq PROPERTIES
        INSTALL_RPATH "$ORIGIN/win64"
        BUILD_WITH_INSTALL_RPATH TRUE
    )
endif()
target_link_libraries(pastaq PRIVATE pastaqlib )

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
