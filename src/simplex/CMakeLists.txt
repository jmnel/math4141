# ---------------------------------------------------
# -- CMakeLists.txt for `arcplot` static library --
# ---------------------------------------------------
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_FLAGS "-Wall -Wextra" )
set( DEBUG_FLAGS "-Wall" "-Wextra" "-O3" )

add_executable( simplex Main.cpp )

set( DEBUG_FLAGS "-Wall" "-Wextra" "-O3" )
target_compile_options( simplex PUBLIC ${DEBUG_FLAGS} )

#include_directories( arcmath )

target_link_libraries( simplex LINK_PUBLIC ${Boost_LIBRARIES} )
target_link_libraries( simplex LINK_PUBLIC arcmath )

# Link libraries here
