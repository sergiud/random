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

enable_testing()

add_test (WELL512a welltest -t WELL512a)
add_test (WELL521a welltest -t WELL521a)
add_test (WELL521b welltest -t WELL521b)
add_test (WELL607a welltest -t WELL607a)
add_test (WELL607b welltest -t WELL607b)
add_test (WELL800a welltest -t WELL800a)
add_test (WELL800b welltest -t WELL800b)
add_test (WELL1024a welltest -t WELL1024a)
add_test (WELL1024b welltest -t WELL1024b)
add_test (WELL19937a welltest -t WELL19937a)
add_test (WELL19937b welltest -t WELL19937b)
add_test (WELL19937c welltest -t WELL19937c)
add_test (WELL21701a welltest -t WELL21701a)
add_test (WELL23209a welltest -t WELL23209a)
add_test (WELL23209b welltest -t WELL23209b)
add_test (WELL44497a welltest -t WELL44497a)
add_test (WELL44497b welltest -t WELL44497b)
