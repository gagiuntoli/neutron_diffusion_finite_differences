cmake_minimum_required(VERSION 3.20)

project(Fermi VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB SOURCE_FILES src/*.cc)

set(SPARSEX_DISABLE_TESTS ON)
add_subdirectory(lib/sparsex)

file(GLOB EXEC_FILES "src/main/main_*.cc")
foreach(EXEC_FILE ${EXEC_FILES})
    get_filename_component(EXEC_NAME ${EXEC_FILE} NAME_WE)
    add_executable(${EXEC_NAME} ${EXEC_FILE})
    target_link_libraries(${EXEC_NAME} PRIVATE sparsex)
endforeach()