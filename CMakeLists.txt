cmake_minimum_required(VERSION 3.20)
project(ceqc VERSION 0.0.1 LANGUAGES CXX)

# CMake/compiler dialect names are ISO C++98/11/14/17/20/23/26; there is no
# portable `CXX_STANDARD 21` value.  The project accepts CEQC_CXX_STANDARD=21
# and maps it to C++20, which is the closest standardized C++21-era mode.
set(CEQC_CXX_STANDARD "21" CACHE STRING "Requested C++ standard: 20, 21, 23")
if(CEQC_CXX_STANDARD STREQUAL "21")
  set(CMAKE_CXX_STANDARD 20)
else()
  set(CMAKE_CXX_STANDARD ${CEQC_CXX_STANDARD})
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CEQC_BUILD_TESTS "Build ceqc C++ tests" ON)
option(CEQC_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

add_library(ceqc_core
  src/app/Application.cpp
  src/cli/CommandLine.cpp
  src/io/Console.cpp
  src/io/Help.cpp
  src/rinex/RinexParser.cpp
  src/rinex/RinexWriter.cpp
  src/rinex/RinexEdit.cpp
  src/rinex/RinexSynthesis.cpp
  src/translate/TranslatorRegistry.cpp
  src/translate/RTCM3Translator.cpp
  src/translate/UBXTranslator.cpp
  src/qc/QCService.cpp
)

target_include_directories(ceqc_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(ceqc_core PUBLIC cxx_std_20)
if(MSVC)
  target_compile_options(ceqc_core PRIVATE /W4 /permissive-)
  if(CEQC_WARNINGS_AS_ERRORS)
    target_compile_options(ceqc_core PRIVATE /WX)
  endif()
else()
  target_compile_options(ceqc_core PRIVATE -Wall -Wextra -Wpedantic)
  if(CEQC_WARNINGS_AS_ERRORS)
    target_compile_options(ceqc_core PRIVATE -Werror)
  endif()
endif()

add_executable(ceqc apps/ceqc/main.cpp)
target_link_libraries(ceqc PRIVATE ceqc_core)
if(MINGW)
  target_link_options(ceqc PRIVATE -static -static-libgcc -static-libstdc++)
endif()

if(CEQC_BUILD_TESTS)
  enable_testing()
  add_executable(ceqc_tests tests/test_main.cpp)
  target_link_libraries(ceqc_tests PRIVATE ceqc_core)
  add_test(NAME ceqc_tests COMMAND ceqc_tests WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

install(TARGETS ceqc RUNTIME DESTINATION bin)
