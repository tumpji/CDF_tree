FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(PythonInterp REQUIRED)



ADD_DEFINITIONS(-Wno-missing-field-initializers )
#ADD_DEFINITIONS(-O3 -ggdb3 -Wall -Wextra -std=c++17 -fno-inline -fno-inline-small-functions -fPIC)
ADD_DEFINITIONS(-O0 -ggdb3 -Wall -Wextra -std=c++17 -fno-inline -fno-inline-small-functions -fPIC)

ADD_LIBRARY(cdftree SHARED numpy_interface.cpp) 
TARGET_INCLUDE_DIRECTORIES(cdftree PUBLIC ${CppNumpyInterface_SOURCE_DIR}/src)
TARGET_INCLUDE_DIRECTORIES(cdftree PUBLIC ${CppNumpyInterface_SOURCE_DIR}/src/utils)

TARGET_INCLUDE_DIRECTORIES(cdftree PRIVATE "/usr/lib64/python3.6/site-packages/numpy/core/include/numpy")
#/usr/lib64/python3.6/site-packages/numpy/core/include/numpy
TARGET_INCLUDE_DIRECTORIES(cdftree PRIVATE ${PYTHON_INCLUDE_DIRS})
# libs
TARGET_LINK_LIBRARIES(cdftree ${PYTHON_LIBRARY})

