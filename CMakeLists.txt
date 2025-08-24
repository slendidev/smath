cmake_minimum_required(VERSION 3.15)
project(SmathExamples CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(smath INTERFACE)
target_include_directories(smath INTERFACE ${CMAKE_SOURCE_DIR}/include)
add_library(smath::smath ALIAS smath)

option(BUILD_EXAMPLES "Build example programs" ON)
if(BUILD_EXAMPLES)
	file(GLOB EXAMPLE_SOURCES "${CMAKE_SOURCE_DIR}/examples/*.cpp")
	foreach(EXAMPLE_FILE ${EXAMPLE_SOURCES})
		get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
		add_executable(${EXAMPLE_NAME} ${EXAMPLE_FILE})
		target_link_libraries(${EXAMPLE_NAME} PRIVATE smath::smath)
	endforeach()
endif()
