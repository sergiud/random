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

target_link_libraries (test_well PRIVATE Boost::unit_test_framework)

add_test (maxeqdist_well19937a test_well -t maxeqdist_well19937a)
add_test (maxeqdist_well44497a test_well -t maxeqdist_well44497a)
add_test (well1024a test_well -t well1024a)
add_test (well1024b test_well -t well1024b)
add_test (well19937a test_well -t well19937a)
add_test (well19937b test_well -t well19937b)
add_test (well19937c test_well -t well19937c)
add_test (well21701a test_well -t well21701a)
add_test (well23209a test_well -t well23209a)
add_test (well23209b test_well -t well23209b)
add_test (well44497a test_well -t well44497a)
add_test (well44497b test_well -t well44497b)
add_test (well512a test_well -t well512a)
add_test (well521a test_well -t well521a)
add_test (well521b test_well -t well521b)
add_test (well607a test_well -t well607a)
add_test (well607b test_well -t well607b)
add_test (well800a test_well -t well800a)
add_test (well800b test_well -t well800b)
