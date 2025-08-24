cmake_minimum_required(VERSION 3.15)
project(SmathExamples CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_TESTS "Build unit tests" ON)

add_library(smath INTERFACE)
target_include_directories(smath INTERFACE ${CMAKE_SOURCE_DIR}/include)
add_library(smath::smath ALIAS smath)

if(BUILD_EXAMPLES)
	file(GLOB EXAMPLE_SOURCES "${CMAKE_SOURCE_DIR}/examples/*.cpp")
	foreach(EXAMPLE_FILE ${EXAMPLE_SOURCES})
		get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
		add_executable(${EXAMPLE_NAME} ${EXAMPLE_FILE})
		target_link_libraries(${EXAMPLE_NAME} PRIVATE smath::smath)
	endforeach()
endif()

if(BUILD_TESTS)
	enable_testing()

	include(FetchContent)
	FetchContent_Declare(
		googletest
		URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
	)
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	FetchContent_MakeAvailable(googletest)

	file(GLOB TEST_SOURCES "${CMAKE_SOURCE_DIR}/tests/*.cpp")

	add_executable(smath_tests ${TEST_SOURCES})
	target_link_libraries(smath_tests PRIVATE
		smath::smath
		GTest::gtest_main
	)

	include(GoogleTest)
	gtest_discover_tests(smath_tests)
endif()
