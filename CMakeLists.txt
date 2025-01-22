cmake_minimum_required(VERSION 3.23)
project(NeutronDiffusionFiniteDifferences)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(lib/sparsex)

file(GLOB EXECUTABLE_SOURCES "src/main/*.cc")

foreach(exe_src IN LISTS EXECUTABLE_SOURCES)
	get_filename_component(exe_name ${exe_src} NAME_WE)
	add_executable(${exe_name} ${exe_src})
	target_link_libraries(${exe_name} sparsex)
	target_include_directories(${exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib)
endforeach()
