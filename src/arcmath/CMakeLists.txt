set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_FLAGS "-Wall -Wextra" )

set( SRC_FILES
    CoreMath.cpp
    MinNelderMead.cpp
    )

add_library( arcmath STATIC ${SRC_FILES} )

# Set library dependency include directories here
target_include_directories( arcmath PUBLIC . )

set( DEBUG_FLAGS "-Wall" "-Wextra" "-O3" )
target_compile_options( arcmath PUBLIC ${DEBUG_FLAGS} )

# Link library dependency libraries here

add_custom_command(
    TARGET arcmath POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/../build/compile_commands.json
    ${CMAKE_SOURCE_DIR}/../compile_commands.json )
