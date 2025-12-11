cmake_minimum_required(VERSION 3.15)
project(SmathExamples LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_EXAMPLES "Build example programs" ON)
option(BUILD_TESTS "Build unit tests" ON)

add_library(smath INTERFACE)
target_include_directories(smath
	INTERFACE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<INSTALL_INTERFACE:include>
)
add_library(smath::smath ALIAS smath)

include(CMakePackageConfigHelpers)
install(TARGETS smath
	EXPORT smathTargets
)
install(DIRECTORY include/ DESTINATION include)

install(EXPORT smathTargets
	NAMESPACE smath::
	FILE smathTargets.cmake
	DESTINATION lib/cmake/smath
)

write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/smathConfigVersion.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/smathConfig.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/smathConfig.cmake"
	INSTALL_DESTINATION "lib/cmake/smath"
)

configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/smath.pc.in"
	"${CMAKE_CURRENT_BINARY_DIR}/smath.pc"
	@ONLY
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/smath.pc"
	DESTINATION lib/pkgconfig
)

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/smathConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/smathConfigVersion.cmake"
	DESTINATION lib/cmake/smath
)

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

	find_package(GTest QUIET)

	if(NOT GTest_FOUND)
		include(FetchContent)
		FetchContent_Declare(
			googletest
			URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
		)
		set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
		FetchContent_MakeAvailable(googletest)
	endif()

	file(GLOB TEST_SOURCES "${CMAKE_SOURCE_DIR}/tests/*.cpp")

	add_executable(smath_tests ${TEST_SOURCES})
	target_link_libraries(smath_tests PRIVATE
		smath::smath
		GTest::gtest_main
	)

	include(GoogleTest)
	gtest_discover_tests(smath_tests)
endif()

