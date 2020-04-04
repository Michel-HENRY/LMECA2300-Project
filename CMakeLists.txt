cmake_minimum_required(VERSION 3.9)
project(ANM C)

file(GLOB_RECURSE sources "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")
add_executable(anm ${sources})
#               "${CMAKE_CURRENT_SOURCE_DIR}/src/main.c"
#               )

target_include_directories(anm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
set_target_properties(anm PROPERTIES
                      C_STANDARD 99
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

# add dependency to BOV
# set(BOV_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) #do not build examples
set(BOV_PARTICLES_VERT "${CMAKE_CURRENT_SOURCE_DIR}/shaders/particles_vert.glsl"
    CACHE STRING "" FORCE)

set(BOV_PARTICLES_GEOM "${CMAKE_CURRENT_SOURCE_DIR}/shaders/particles_geom.glsl"
    CACHE STRING "" FORCE)

set(BOV_PARTICLES_FRAG "${CMAKE_CURRENT_SOURCE_DIR}/shaders/particles_frag.glsl"
    CACHE STRING "" FORCE)

set(BOV_PARTICLES_POSTPRO "${CMAKE_CURRENT_SOURCE_DIR}/shaders/postpro.glsl"
    CACHE STRING "" FORCE)

add_subdirectory(deps/BOV)
target_link_libraries(anm 
                      PUBLIC bov)

# set anm as the startup project in visual studio
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT anm)
