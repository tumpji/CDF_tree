find_package(Boost COMPONENTS unit_test_framework REQUIRED)



#ADD_DEFINITIONS(-Wall -Wno-deprecated-declarations -Wextra -O0 -ggdb3 -fno-inline -fno-inline-small-functions) 
#ADD_DEFINITIONS(-Wall -Wextra -O2 -ggdb)
ADD_DEFINITIONS(-Wall -Wextra -O0 -ggdb3 -fno-inline -fno-inline-small-functions) # -fno-implicit-templates)
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK) 

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src ${TEST_SOURCE_DIR} ${Boost_INCLUDE_DIRS})

ADD_EXECUTABLE(test_tree test_tree.cpp)
TARGET_LINK_LIBRARIES(test_tree ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
ADD_TEST(UnitTest test_tree)

# interface tests
ADD_CUSTOM_TARGET(
    prepare_python_test_libraries
    COMMAND "${CMAKE_COMMAND}" -E copy ${CMAKE_BINARY_DIR}/src/libcdftree.so ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS cdftree
    VERBATIM
    )
ADD_CUSTOM_TARGET(run_interface_test
    COMMAND ${PYTHON_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS prepare_python_test_libraries
    VERBATIM
    )

#ADD_TEST(NAME python_interface_test
#    COMMAND run_interface_test
#    VERBATIM
#    )
# speed test
#ADD_DEFINITIONS(-O3 -march=native)
ADD_EXECUTABLE(stest speed_test.cpp)


