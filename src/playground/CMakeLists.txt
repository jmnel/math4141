# -------------------------------------------------------
# -- CMakeLists.txt for `playground` static executable --
# -------------------------------------------------------
set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_FLAGS "-Wall -Wextra" )
set( DEBUG_FLAGS "-Wall" "-Wextra" "-O3" )

add_executable( playground Main.cpp )

set( DEBUG_FLAGS "-Wall" "-Wextra" "-O3" )
target_compile_options( playground PUBLIC ${DEBUG_FLAGS} )

target_link_libraries( playground arcmath )

# Link libraries here
