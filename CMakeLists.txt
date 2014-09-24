cmake_minimum_required (VERSION 2.6 FATAL_ERROR)
enable_testing()
project (DynamO)

#Check that C++11 or C++0x is supported!
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(SEND_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++0x/C++11 support. Please use a more modern C++ compiler.")
endif()

#Check that boost is available
find_package(Boost 1.46 REQUIRED COMPONENTS "program_options" "filesystem" "iostreams" "system")

if(NOT Boost_FOUND)
  message(STATUS "Cannot build, missing boost")
endif()

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/src/dynamo/)
include_directories(${PROJECT_SOURCE_DIR}/src/coil/)
include_directories(${PROJECT_SOURCE_DIR}/src/magnet)
include_directories(${PROJECT_SOURCE_DIR}/src/ext_include)

#Dynamo simulation executables
file(GLOB_RECURSE dynamo_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/dynamo/*.cpp)
add_library(dynamo OBJECT ${dynamo_SRC})

add_executable(dynarun ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarun.cpp $<TARGET_OBJECTS:dynamo>)
target_link_libraries(dynarun ${Boost_LIBRARIES})

add_executable(dynamod ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynamod.cpp $<TARGET_OBJECTS:dynamo>)
target_link_libraries(dynamod ${Boost_LIBRARIES})

add_executable(dynahist_rw ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynahist_rw.cpp $<TARGET_OBJECTS:dynamo>)
target_link_libraries(dynahist_rw ${Boost_LIBRARIES})


find_package(PkgConfig)
pkg_check_modules(GTKMM gtkmm-2.4) # look into FindPkgConfig.cmake, 
if(GTKMM_FOUND)
  #Coil visualiser
  file(GLOB_RECURSE coil_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/coil/coil/*.cpp)
  add_library(coil OBJECT ${coil_SRC})
  add_executable(dynavis ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarun.cpp $<TARGET_OBJECTS:dynamo> $<TARGET_OBJECTS:coil>)
  link_directories(dynavis ${GTKMM_LIBRARY_DIRS})
  include_directories(dynavis  ${GTKMM_INCLUDE_DIRS})
  target_link_libraries(dynavis ${Boost_LIBRARIES} gtkmm)
endif()

