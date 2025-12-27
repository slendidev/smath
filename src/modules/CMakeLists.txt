add_library(smath_modules)

target_sources(smath_modules PUBLIC FILE_SET CXX_MODULES FILES smath.cppm)

target_link_libraries(smath_modules PUBLIC smath::smath)

target_compile_features(smath_modules PUBLIC cxx_std_20)

