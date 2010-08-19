cmake_minimum_required (VERSION 2.8)
project (Random)

find_package (Boost 1.43.0 REQUIRED unit_test_framework)

include_directories (${Boost_INCLUDE_DIRS})

if (NOT Boost_USE_STATIC_LIBS)
  add_definitions (-DBOOST_TEST_DYN_LINK)
endif ()

if (CMAKE_COMPILER_IS_GNUCXX)
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic")
endif ()

add_executable (
  welltest
  welltest.cpp
  well.hpp
)

target_link_libraries (welltest ${Boost_LIBRARIES})

