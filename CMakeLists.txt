cmake_minimum_required(VERSION 3.10)
project(CDFtree)


enable_testing()

if(CMAKE_PROJECT_NAME STREQUAL CDFtree)
    add_subdirectory(test)
endif()

add_subdirectory(src)
add_subdirectory(third_party)

#ADD_DEFINITIONS(-Wall -Wno-deprecated-declarations -Wextra -O0 -ggdb3 -fno-inline -fno-inline-small-functions) 
#ADD_DEFINITIONS(-Wall -Wextra -O0 -ggdb3 -fno-inline -fno-inline-small-functions) 
#ADD_DEFINITIONS(-Wall -Wextra -O2)
#
#ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 
#
#find_package(Boost COMPONENTS unit_test_framework REQUIRED)
#ADD_EXECUTABLE(test test.cpp Bdistribution.h)
#TARGET_LINK_LIBRARIES(test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#
#ADD_TEST(UnitTest test)

