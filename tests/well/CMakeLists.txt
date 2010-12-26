cmake_minimum_required (VERSION 2.8)
project (test-well)

find_package (Boost REQUIRED unit_test_framework)

include_directories (${Boost_INCLUDE_DIRS})
include_directories (${CMAKE_SOURCE_DIR})

add_executable (
  test-well
  main.cpp
  ${CMAKE_SOURCE_DIR}/well.hpp
)

target_link_libraries (test-well ${Boost_LIBRARIES})

add_test (WELL512a test-well -t WELL512a)
add_test (WELL521a test-well -t WELL521a)
add_test (WELL521b test-well -t WELL521b)
add_test (WELL607a test-well -t WELL607a)
add_test (WELL607b test-well -t WELL607b)
add_test (WELL800a test-well -t WELL800a)
add_test (WELL800b test-well -t WELL800b)
add_test (WELL1024a test-well -t WELL1024a)
add_test (WELL1024b test-well -t WELL1024b)
add_test (WELL19937a test-well -t WELL19937a)
add_test (WELL19937b test-well -t WELL19937b)
add_test (WELL19937c test-well -t WELL19937c)
add_test (WELL21701a test-well -t WELL21701a)
add_test (WELL23209a test-well -t WELL23209a)
add_test (WELL23209b test-well -t WELL23209b)
add_test (WELL44497a test-well -t WELL44497a)
add_test (WELL44497b test-well -t WELL44497b)
