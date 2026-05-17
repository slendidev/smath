cmake_minimum_required(VERSION 3.15)
project(SmathExamples LANGUAGES CXX VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_GENERATOR MATCHES "^Visual Studio")
include(GNUInstallDirs)
include(FindPkgConfig)
endif()

option(SMATH_BUILD_EXAMPLES "Build example programs" ON)
option(SMATH_BUILD_GLFW_EXAMPLES "Build GLFW example programs" OFF)
option(SMATH_BUILD_GLUT_EXAMPLES "Build GLUT example programs" OFF)
option(SMATH_BUILD_TESTS "Build unit tests" ON)
option(SMATH_BUILD_MODULES "Enable C++20 modules support" OFF)
option(SMATH_ENABLE_LEGACY_HEADER "Install legacy include <smath.hpp> shim" OFF)
option(SMATH_INSTALL_INTEROP_HEADERS "Install interop headers" OFF)

if( NOT CMAKE_BUILD_TYPE )
        set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE )
endif()

if(NOT CMAKE_GENERATOR MATCHES "^Visual Studio")
if( NOT CMAKE_BUILD_TYPE STREQUAL Debug )
add_compile_options(-O3)
add_link_options(-Wl,--gc-sections)
add_link_options(-Wl,--print-gc-sections)
add_link_options(-Wl,-s)
else()
add_compile_options(-O0 -g)
endif()
else()
if( NOT CMAKE_BUILD_TYPE STREQUAL Debug )
add_compile_options(/O2)
else()
add_compile_options(/Od /Zi)
endif()
endif()

if (NOT CMAKE_GENERATOR MATCHES "^Visual Studio")
if ( NOT CMAKE_CROSS_COMPILING )
add_compile_options(-march=native)
if( CMAKE_GENERATOR_PLATFORM STREQUAL "x64" OR CMAKE_GENERATOR_PLATFORM STREQUAL "x86_64" OR CMAKE_GENERATOR_PLATFORM STREQUAL "AMD64" )
	add_compile_options(-mfpmath=sse)
endif()
endif()
add_compile_options(-fopenmp)
add_compile_options(-fopenmp-simd)
add_compile_options(-ftree-vectorize)
add_compile_options(-Wno-gnu-alignof-expression)
add_compile_options(-D_FILE_OFFSET_BITS=64)
add_compile_options(-fdata-sections)
add_compile_options(-ffunction-sections)
add_compile_options(-fpermissive)
add_compile_options(-Wno-array-bounds)
add_compile_options(-Wno-unused-parameter)
add_compile_options(-Wno-pedantic)
add_compile_options(-Wno-narrowing)
add_compile_options(-Wno-attributes)
else()
add_compile_options(/openmp:experimental)
endif()

if(SMATH_BUILD_MODULES)
	message(STATUS "Building smath C++ module")
	add_subdirectory(src/modules)
endif()

add_library(smath INTERFACE)
target_include_directories(smath
	INTERFACE
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
		$<$<BOOL:${SMATH_ENABLE_LEGACY_HEADER}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>>
		$<INSTALL_INTERFACE:include>
)
add_library(smath::smath ALIAS smath)

include(CMakePackageConfigHelpers)
install(TARGETS smath
	EXPORT smathTargets
)
install(FILES include/smath/smath.hpp DESTINATION include/smath)

if(SMATH_INSTALL_INTEROP_HEADERS)
	install(DIRECTORY include/smath/interop DESTINATION include/smath)
endif()

if(SMATH_ENABLE_LEGACY_HEADER)
	file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
	configure_file(
		"${CMAKE_CURRENT_SOURCE_DIR}/cmake/smath_legacy_header.hpp.in"
		"${CMAKE_CURRENT_BINARY_DIR}/include/smath.hpp"
		@ONLY
	)
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/include/smath.hpp" DESTINATION include)
endif()

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

if(SMATH_BUILD_EXAMPLES)
	file(GLOB EXAMPLE_SOURCES "${CMAKE_SOURCE_DIR}/examples/*.cpp")
	foreach(EXAMPLE_FILE ${EXAMPLE_SOURCES})
		get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
		add_executable(${EXAMPLE_NAME} ${EXAMPLE_FILE})
		target_link_libraries(${EXAMPLE_NAME} PRIVATE smath::smath)
	endforeach()
endif()

if(SMATH_BUILD_GLFW_EXAMPLES)
	find_package(OpenGL REQUIRED)
	find_package(glfw3 REQUIRED)
	file(GLOB GLFW_EXAMPLE_SOURCES "${CMAKE_SOURCE_DIR}/examples/glfw/*.cpp")
	foreach(EXAMPLE_FILE ${GLFW_EXAMPLE_SOURCES})
		get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
		add_executable("glfw_${EXAMPLE_NAME}" ${EXAMPLE_FILE})
		target_link_libraries("glfw_${EXAMPLE_NAME}" PUBLIC smath::smath glfw OpenGL OpenGL::GLU)
	endforeach()
endif()

if(SMATH_BUILD_GLUT_EXAMPLES)
	find_package(OpenGL REQUIRED)
	find_package(GLUT REQUIRED)
	file(GLOB GLUT_EXAMPLE_SOURCES "${CMAKE_SOURCE_DIR}/examples/glut/*.cpp")
	foreach(EXAMPLE_FILE ${GLUT_EXAMPLE_SOURCES})
		get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
		add_executable("glut_${EXAMPLE_NAME}" ${EXAMPLE_FILE})
		target_link_libraries("glut_${EXAMPLE_NAME}" PUBLIC smath::smath OpenGL OpenGL::GLU GLUT::GLUT)
	endforeach()
endif()

if(SMATH_BUILD_TESTS)
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
