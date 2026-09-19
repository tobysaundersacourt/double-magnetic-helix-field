cmake_minimum_required(VERSION 3.16)
project(AnalyticMagneticField LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Boost is used only through headers (special functions); there are no Boost
# libraries to link. BOOST_ROOT / Boost_ROOT may point either to an unpacked
# Boost tree (containing boost/) or to an installation prefix (containing
# include/boost/). Avoiding find_package(Boost) keeps this equally usable with
# system packages, Homebrew, vcpkg-style layouts, and a plain Boost archive.
find_path(BOOST_INCLUDE_DIR
    NAMES boost/math/special_functions/bessel.hpp
    HINTS ${BOOST_ROOT} ${Boost_ROOT}
    PATH_SUFFIXES include
)
if(NOT BOOST_INCLUDE_DIR)
    message(FATAL_ERROR "Boost headers were not found. Set BOOST_ROOT or Boost_ROOT to the directory containing boost/.")
endif()
message(STATUS "Boost headers: ${BOOST_INCLUDE_DIR}")

add_executable(field_solver main.cpp)
target_include_directories(field_solver PRIVATE ${BOOST_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

if(MSVC)
    target_compile_options(field_solver PRIVATE /O2 /EHsc)
    # Keep the native solver self-contained on Windows so users do not need to
    # install the matching Visual C++ runtime just for this executable.
    set_property(TARGET field_solver PROPERTY
        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
else()
    target_compile_options(field_solver PRIVATE -O3)
endif()

# Put the solver in a predictable place for single-config generators.
set_target_properties(field_solver PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
