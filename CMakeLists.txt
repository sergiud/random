cmake_minimum_required (VERSION 2.8)
project (Random)

enable_testing()

find_package (Boost 1.43.0 REQUIRED unit_test_framework)

include_directories (${Boost_INCLUDE_DIRS})

if (NOT Boost_USE_STATIC_LIBS)
  add_definitions (-DBOOST_TEST_DYN_LINK)
endif (NOT Boost_USE_STATIC_LIBS)

if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic")
endif (CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory (tests)
