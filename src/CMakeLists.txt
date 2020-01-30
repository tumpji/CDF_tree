FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(PythonInterp REQUIRED)



ADD_DEFINITIONS(-Wno-missing-field-initializers )
#ADD_DEFINITIONS(-O3 -ggdb3 -Wall -Wextra -std=c++17 -fno-inline -fno-inline-small-functions -fPIC)
ADD_DEFINITIONS(-O0 -ggdb3 -Wall -Wextra -std=c++17 -fno-inline -fno-inline-small-functions -fPIC)

ADD_LIBRARY(cdftree SHARED numpy_interface.cpp) 
TARGET_INCLUDE_DIRECTORIES(cdftree PRIVATE ${CppNumpyInterface_INCLUDE_DIRS})
# libs
TARGET_LINK_LIBRARIES(cdftree ${CppNumpyInterface_LINK_LIBRARIES})

