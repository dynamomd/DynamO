cmake_minimum_required (VERSION 2.6 FATAL_ERROR)
enable_testing()
project (DynamO)


##########   CONFIGURATION 
#C++11/C++0x support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(SEND_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++0x/C++11 support. Please use a more modern C++ compiler.")
  endif()
endif()

# Boost support
find_package(Boost 1.46 REQUIRED COMPONENTS "program_options" "filesystem" "iostreams" "system")
if(NOT Boost_FOUND)
  message(SEND_ERROR "Cannot build, missing boost")
endif()
include_directories(${Boost_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/src/dynamo/)
include_directories(${PROJECT_SOURCE_DIR}/src/coil/)
include_directories(${PROJECT_SOURCE_DIR}/src/magnet)
include_directories(${PROJECT_SOURCE_DIR}/src/ext_include)

######### BUILD TARGETS
# DynamO files
file(GLOB_RECURSE dynamo_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/dynamo/*.cpp)
add_library(dynamo STATIC ${dynamo_SRC})

add_executable(dynarun ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarun.cpp)
target_link_libraries(dynarun ${Boost_LIBRARIES} dynamo)

add_executable(dynamod ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynamod.cpp)
target_link_libraries(dynamod ${Boost_LIBRARIES} dynamo)

add_executable(dynahist_rw ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynahist_rw.cpp)
target_link_libraries(dynahist_rw ${Boost_LIBRARIES} dynamo)

find_package(PkgConfig)
find_package(OpenGL)
find_package(GLUT)
find_package(GLEW)
pkg_check_modules(GTKMM gtkmm-2.4) # look into FindPkgConfig.cmake
pkg_check_modules(GLIB glibmm-2.4) # look into FindPkgConfig.cmake
pkg_check_modules(PNG libpng)

if(((((GTKMM_FOUND AND PNG_FOUND) AND OPENGL_FOUND) AND GLUT_FOUND) AND GLEW_FOUND) AND GLIB_FOUND)
  message(STATUS "Visualiser dependencies satisfied! Building dynavis")
  
  #Coil library
  set_target_properties(dynamo PROPERTIES COMPILE_DEFINITIONS "DYNAMO_visualizer")
  file(GLOB_RECURSE coil_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/coil/coil/*.cpp)
  add_library(coil STATIC ${coil_SRC})

  add_executable(dynavis ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarun.cpp)
  target_link_libraries(dynavis ${Boost_LIBRARIES} dynamo)
  include_directories(${GTKMM_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
  target_link_libraries(dynamo coil ${Boost_LIBRARIES} ${GTKMM_LIBRARIES} ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLEW_LIBRARIES} ${GLIB_LIBRARIES})
  add_definitions(${OPENGL_DEFINITONS} ${GLUT_DEFINITIONS} ${GLEW_DEFINITIONS})
endif()
