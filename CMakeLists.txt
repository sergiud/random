cmake_minimum_required (VERSION 3.15)
project (Random LANGUAGES CXX)

include (CTest)

find_package (Boost 1.43.0 REQUIRED unit_test_framework NO_MODULE)

add_executable (test_well
  tests/well/main.cpp
  well.hpp
)
target_include_directories (test_well PUBLIC
  $<BUILD_INTERFACE:${Random_SOURCE_DIR}>
)

target_link_libraries (test_well ${Boost_LIBRARIES})

add_test (WELL512a test_well -t WELL512a)
add_test (WELL521a test_well -t WELL521a)
add_test (WELL521b test_well -t WELL521b)
add_test (WELL607a test_well -t WELL607a)
add_test (WELL607b test_well -t WELL607b)
add_test (WELL800a test_well -t WELL800a)
add_test (WELL800b test_well -t WELL800b)
add_test (WELL1024a test_well -t WELL1024a)
add_test (WELL1024b test_well -t WELL1024b)
add_test (WELL19937a test_well -t WELL19937a)
add_test (WELL19937b test_well -t WELL19937b)
add_test (WELL19937c test_well -t WELL19937c)
add_test (WELL21701a test_well -t WELL21701a)
add_test (WELL23209a test_well -t WELL23209a)
add_test (WELL23209b test_well -t WELL23209b)
add_test (WELL44497a test_well -t WELL44497a)
add_test (WELL44497b test_well -t WELL44497b)
