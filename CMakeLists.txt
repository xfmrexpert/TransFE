cmake_minimum_required(VERSION 3.12)
project(TransFE VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)
enable_testing()

# Include the root directory so headers referenced relatively work as expected
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Automatically gather sources from root and MeshDB
# CONFIGURE_DEPENDS ensures CMake re-runs if files are added/removed
file(GLOB ROOT_SOURCES CONFIGURE_DEPENDS "*.cpp")
file(GLOB MESHDB_SOURCES CONFIGURE_DEPENDS "MeshDB/*.cpp")

# Combine sources, deliberately excluding the "Old" directory
set(SOURCES ${ROOT_SOURCES} ${MESHDB_SOURCES})

add_executable(TransFE ${SOURCES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
