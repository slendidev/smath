cmake_minimum_required(VERSION 3.15)
project(SmathExamples LANGUAGES CXX VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

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

string(TOUPPER "${CMAKE_SYSTEM_PROCESSOR}" SYSTEM_PROCESSOR_UPPER)
message(STATUS "OS: ${CMAKE_HOST_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR} => ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "CXX: ${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
add_compile_definitions(_FILE_OFFSET_BITS=64)
add_compile_definitions($<$<CONFIG:Debug>:DEBUG=1> $<$<NOT:$<CONFIG:Debug>>:NDEBUG=1>)
add_compile_definitions($<$<CONFIG:Debug>:_DEBUG=1> $<$<NOT:$<CONFIG:Debug>>:_NDEBUG=1>)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wno-gnu-alignof-expression)
endif()

if (CMAKE_CXX_COMPILER_FRONTEND_VARIANT MATCHES "GNU|AppleClang")
  if( NOT CMAKE_BUILD_OPTLEVEL )
    set( CMAKE_BUILD_OPTLEVEL 3 CACHE STRING "Choose the optimization level." FORCE )
  endif()
  add_compile_options(-O${CMAKE_BUILD_OPTLEVEL})
  if (NOT CMAKE_CROSSCOMPILING)
    add_compile_options(-march=native)

    if (${SYSTEM_PROCESSOR_UPPER} MATCHES "X64|X86_64|AMD64")
      add_compile_options(-mfpmath=sse)
    endif()
    if (${SYSTEM_PROCESSOR_UPPER} MATCHES "PPC|POWERPC|PPC64|PPC64LE")
      include(CheckCSourceCompiles)

      set(CMAKE_REQUIRED_FLAGS "-maltivec")
      check_c_source_compiles("#include <altivec.h> int main() { vector int v = (vector int){0}; return 0; }" HAVE_ALTIVEC)

      if(HAVE_ALTIVEC)
        add_compile_options(-maltivec)
      endif()
    endif()
  endif()
  add_compile_options(-ftree-vectorize)
  add_compile_options(-fopenmp)
  if(NOT CMAKE_CXX_COMPILER_FRONTEND_VARIANT MATCHES "AppleClang")
    add_compile_options(-fopenmp-simd)
  endif()
  add_compile_options(-fdata-sections)
  add_compile_options(-fpermissive)
  add_compile_options(-ffunction-sections)
  add_compile_options(-fno-asynchronous-unwind-tables)
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options(-g)
  else()
    # Use IPO/LTO if available
    if(NOT ${CMAKE_VERSION} VERSION_LESS "3.9.0")
      cmake_policy(SET CMP0069 NEW)
      include(CheckIPOSupported)
      check_ipo_supported(RESULT result)
      if(result)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
        message(STATUS "Using IPO/LTO.")
      endif()
    endif()
  endif()
  if (CMAKE_BUILD_TYPE MATCHES ReleaseStrip)
    add_link_options(-Wl,--gc-sections)
    add_link_options(-Wl,-s)
  endif()
  link_libraries(m)
  if (WIN32 AND MINGW)
    link_libraries(mingw32)
  endif()
else (CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
  if( NOT CMAKE_BUILD_OPTLEVEL )
    set( CMAKE_BUILD_OPTLEVEL 2 CACHE STRING "Choose the optimization level." FORCE )
  endif()
  add_compile_options(/O${CMAKE_BUILD_OPTLEVEL})
  if (NOT CMAKE_CROSSCOMPILING)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS _WIN32)
    if (CMAKE_VS_PLATFORM_NAME MATCHES "x64|Win32")
      add_compile_options(/Gv)
    endif()
  endif()
  add_compile_options(/openmp:experimental)
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options(/Zi)
    add_link_options(/DEBUG:FULL)
  endif()
  if (CMAKE_BUILD_TYPE MATCHES ReleaseStrip)
    add_link_options(/DEBUG:NONE)
  endif()
endif()

if(SMATH_BUILD_MODULES)
  message(STATUS "Building smath C++ module")
  add_subdirectory(src/modules)
endif()

add_library(smath INTERFACE)
target_include_directories(smath
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<$<BOOL:${SMATH_ENABLE_LEGACY_HEADER}>:$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}>>
  $<INSTALL_INTERFACE:include>
)
add_library(smath::smath ALIAS smath)

include(CMakePackageConfigHelpers)
install(TARGETS smath
  EXPORT smathTargets
)
install(FILES include/smath/smath.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smath)

if(SMATH_INSTALL_INTEROP_HEADERS)
  install(DIRECTORY include/smath/interop DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smath)
endif()

if(SMATH_ENABLE_LEGACY_HEADER)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/smath_legacy_header.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/smath.hpp"
    @ONLY
  )
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/smath.hpp" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

install(EXPORT smathTargets
  NAMESPACE smath::
  FILE smathTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/smath
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/smathConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/smathConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/smathConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/smath"
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/smath.pc.in"
  "${CMAKE_CURRENT_BINARY_DIR}/smath.pc"
  @ONLY
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/smath.pc"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/smathConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/smathConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/smath
)

if(SMATH_BUILD_EXAMPLES)
  file(GLOB EXAMPLE_SOURCES "${CMAKE_SOURCE_DIR}/examples/*.cpp")
  foreach(EXAMPLE_FILE ${EXAMPLE_SOURCES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_FILE})
    target_include_directories(${EXAMPLE_NAME} PRIVATE "${CMAKE_SOURCE_DIR}/examples")
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
    target_include_directories("glfw_${EXAMPLE_NAME}" PRIVATE "${CMAKE_SOURCE_DIR}/examples/glfw")
    target_link_libraries("glfw_${EXAMPLE_NAME}" PRIVATE smath::smath PUBLIC glfw OpenGL OpenGL::GLU)
  endforeach()
endif()

if(SMATH_BUILD_GLUT_EXAMPLES)
  find_package(OpenGL REQUIRED)
  find_package(GLUT REQUIRED)
  file(GLOB GLUT_EXAMPLE_SOURCES "${CMAKE_SOURCE_DIR}/examples/glut/*.cpp")
  foreach(EXAMPLE_FILE ${GLUT_EXAMPLE_SOURCES})
    get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
    add_executable("glut_${EXAMPLE_NAME}" ${EXAMPLE_FILE})
    target_include_directories("glut_${EXAMPLE_NAME}" PRIVATE "${CMAKE_SOURCE_DIR}/examples/glut")
    target_link_libraries("glut_${EXAMPLE_NAME}" PRIVATE smath::smath PUBLIC OpenGL OpenGL::GLU GLUT::GLUT)
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
