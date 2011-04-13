cmake_minimum_required (VERSION 2.8)
project (test_well)

find_package (Boost REQUIRED unit_test_framework)

include_directories (${Boost_INCLUDE_DIRS})
include_directories (${CMAKE_SOURCE_DIR})

add_executable (
  ${PROJECT_NAME}
  main.cpp
  ${CMAKE_SOURCE_DIR}/well.hpp
)

target_link_libraries (${PROJECT_NAME} ${Boost_LIBRARIES})

add_test (WELL512a ${PROJECT_NAME} -t WELL512a)
add_test (WELL521a ${PROJECT_NAME} -t WELL521a)
add_test (WELL521b ${PROJECT_NAME} -t WELL521b)
add_test (WELL607a ${PROJECT_NAME} -t WELL607a)
add_test (WELL607b ${PROJECT_NAME} -t WELL607b)
add_test (WELL800a ${PROJECT_NAME} -t WELL800a)
add_test (WELL800b ${PROJECT_NAME} -t WELL800b)
add_test (WELL1024a ${PROJECT_NAME} -t WELL1024a)
add_test (WELL1024b ${PROJECT_NAME} -t WELL1024b)
add_test (WELL19937a ${PROJECT_NAME} -t WELL19937a)
add_test (WELL19937b ${PROJECT_NAME} -t WELL19937b)
add_test (WELL19937c ${PROJECT_NAME} -t WELL19937c)
add_test (WELL21701a ${PROJECT_NAME} -t WELL21701a)
add_test (WELL23209a ${PROJECT_NAME} -t WELL23209a)
add_test (WELL23209b ${PROJECT_NAME} -t WELL23209b)
add_test (WELL44497a ${PROJECT_NAME} -t WELL44497a)
add_test (WELL44497b ${PROJECT_NAME} -t WELL44497b)
