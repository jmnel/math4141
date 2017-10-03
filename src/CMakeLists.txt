cmake_minimum_required( VERSION 3.4 FATAL_ERROR )
project( Math4141 VERSION 0.1 LANGUAGES CXX )

#set( CMAKE_CXX_COMPILER_LAUNCHER ccache )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${Math4141_SOURCE_DIR}/../bin/Debug )
set( CMAKE_MODULE_PATH ${Math4141_SOURCE_DIR}/cmake/modules )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

# Run `find_package` modules here for dependencies
find_package(Boost 1.64 COMPONENTS iostreams system filesystem REQUIRED )
#find_package(PythonLibs 3.6.2 REQUIRED )
#target_include_directories( Math4141 PRIVATE ${PYTHON_INCLUDE_DIRS} )

add_subdirectory( ./simplex )
add_subdirectory( ./arcmath )
add_subdirectory( ./playground )
