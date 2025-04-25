cmake_minimum_required(VERSION 3.15)
project(SubdividerPlugin LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ----- Use VCPKG if present -----
# You should pass this from the CLI:
# export NUKE_INCLUDE_DIR=/path/Nuke/Nuke15.1v3/include
# export NUKE_LIB_DIR=/path/metaphysic/Nuke/Nuke15.1v3
# mkdir build && cd build
# cmake .. -DCMAKE_TOOLCHAIN_FILE=/path/vcpkg/scripts/buildsystems/vcpkg.cmake
# make

# ----- Set Nuke and OpenSubdiv paths -----
set(NUKE_INCLUDE_DIR $ENV{NUKE_INCLUDE_DIR} CACHE PATH "Path to Nuke headers")
set(NUKE_LIB_DIR $ENV{NUKE_LIB_DIR} CACHE PATH "Path to Nuke libraries (usually same as plugins dir)")

# ----- Include directories -----
include_directories(
    ${NUKE_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# ----- Source files -----
set(SOURCES
    Subdivider.cpp
    MeshSubdivider.cpp
)

# ----- Build shared object plugin -----
add_library(MeshSubdivider SHARED ${SOURCES})

# ----- Compile flags -----
target_compile_definitions(MeshSubdivider PRIVATE __GNUC__)
target_compile_options(MeshSubdivider PRIVATE -fPIC -Wall -Wextra)

# ----- Link to OpenSubdiv -----
find_package(OpenSubdiv REQUIRED)

target_link_libraries(MeshSubdivider PRIVATE OpenSubdiv::osdCPU_static)


# Define OpenSubdiv include path manually
# set(OPENSUBDIV_INCLUDE_DIR "${VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}/include")
set(OPENSUBDIV_INCLUDE_DIR "/home/yuri.carrara.r/vcpkg/installed/x64-linux/include")

message(STATUS ${OPENSUBDIV_INCLUDE_DIR})

# Add the include directory directly
target_include_directories(MeshSubdivider PRIVATE ${OPENSUBDIV_INCLUDE_DIR})

# ----- Link to Nuke's libDDImage if needed -----
target_link_directories(MeshSubdivider PRIVATE ${NUKE_LIB_DIR})
target_link_libraries(MeshSubdivider PRIVATE ${NUKE_LIB_DIR}/libDDImage.so)

# ----- Output plugin name and path -----
set_target_properties(MeshSubdivider PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist
    PREFIX ""   # don't prefix with 'lib'
    SUFFIX ".so"
)
