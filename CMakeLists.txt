cmake_minimum_required (VERSION 3.0 FATAL_ERROR)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(WIN32)
  add_definitions(-D_USE_MATH_DEFINES)
  add_definitions(-DNOMINMAX)
endif()

IF(DEFINED CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
ENDIF()

##########   OPTIONS
project(DynamO) #Project name
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

enable_testing() #Enable build of test executables and 'make test' command
include(CTest)

######################################################################
########## Packaging
######################################################################
include(SystemInformation)
include(GetGitRevisionDescription)

find_package(Python3 COMPONENTS Interpreter NumPy)

git_describe(GIT_RESULT --match initial-commit)
message(STATUS "Git description is ${GIT_RESULT}")
string(REGEX MATCH "[0-9]+" BUILD_NUMBER ${GIT_RESULT})
message(STATUS "Build number is ${BUILD_NUMBER}")
if(NOT DEFINED CPACK_GENERATOR)
  set(CPACK_GENERATOR ${SPECIFIC_SYSTEM_PREFERED_CPACK_GENERATOR})
endif()
set(PACKAGE_SUFFIX "" CACHE STRING "A suffix to add to the end of the built packages")
set(CPACK_PACKAGE_NAME "dynamomd")
set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 7)
set(CPACK_PACKAGE_VERSION_PATCH ${BUILD_NUMBER})
set(CPACK_PACKAGE_CONTACT "Marcus Bannerman <support@dynamomd.org>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Event-driven particle simulation sofware")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_SOURCE_DIR}/description.txt)
set(CPACK_PACKAGE_DESCRIPTION "This package includes the DynamO executables and scripts. DynamO is an event driven simulator capable of performing molecular and granular dynamics simulations")
if(Python3_FOUND)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS " python3")
endif()
set(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}${PACKAGE_SUFFIX})
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://dynamomd.org")
#Get Debian packages to automatically figure out the dependencies
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)

set(CPACK_WIX_UPGRADE_GUID "c70e546c-a36a-429f-a2df-f5e6cf787bb4")
set(CPACK_WIX_PRODUCT_GUID "36fe5a30-bd03-423b-822a-f505d0fa6564")


include(CPack) #Enable build of packages

install(FILES ${CMAKE_SOURCE_DIR}/copyright DESTINATION share/doc/${CPACK_PACKAGE_NAME})

##########   RELEASE MODE
if(MSVC)
  #MSVC has crazy warnings for -Wall, we'll build up support to the higher warning levels
  add_compile_options(-W1)
  #Enable the math defines (like M_PI!)
  add_definitions(-D_USE_MATH_DEFINES)
  #MSVC Crazily has "min" and "max" (yes, lowercase!) Macros defined
  #which conflict with the C++ standard algorithms! This define
  #instructs the compiler to not be an idiot and define those macros
  add_definitions(-DNOMINMAX)
else()
  add_compile_options(-Wall)
endif()


##########   DEBUG MODE
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDYNAMO_DEBUG -DMAGNET_DEBUG")
#  -D_GLIBCXX_DEBUG <- Not compatible with boost program_options

if(NOT WIN32)
  link_libraries(debug dl)
endif()

######################################################################
########## COMPILER C++11/C++0x SUPPORT TESTS
######################################################################
### First check if the compiler supports C++11 or C++0x at all
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORT_CXX11)
if(COMPILER_SUPPORT_CXX11)
  add_compile_options(-std=c++11)
else()
  check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORT_CXX0X)
  if(COMPILER_SUPPORT_CXX0X)
    add_compile_options(-std=c++0x)
    check_cxx_source_compiles("int main(int, const char**) {int array[5] = {1,2,3,4,5}; for(int& x: array) x *= 2;}" CXX0X_RANGE_BASED_FOR)
    if(NOT CXX0X_RANGE_BASED_FOR)
      message(SEND_ERROR "The compiler ${CMAKE_CXX_COMPILER} has limited C++11 support. You must install a more modern C++ compiler.")
    endif()
  else()
    message(SEND_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. You must install a more modern C++ compiler.")
  endif()
endif()

include_directories(./extern/stator/)
include_directories(./extern/)

######################################################################
# Stack trace support
######################################################################
include(CheckIncludeFiles)
check_include_files(execinfo.h STACKTRACE_AVAILABLE)
if(STACKTRACE_AVAILABLE)
  add_definitions(-DMAGNET_STACKTRACE)  
endif()

######################################################################
# Test for libbz2 (for compressed files)
######################################################################
if(NOT WIN32)
  find_package(BZip2)
  if(BZIP2_FOUND)
    include_directories(${BZIP2_INCLUDE_DIR})
    link_libraries(${BZIP2_LIBRARIES})
    add_definitions(-DDYNAMO_bzip2_support)
  endif()
endif()

######################################################################
##########  Boost support
######################################################################
if(WIN32)
  set (Boost_USE_STATIC_LIBS ON)
endif()
find_package(Boost 1.46 REQUIRED COMPONENTS program_options filesystem system unit_test_framework)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_FILESYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_IOSTREAMS_LIBRARY} ${Boost_SYSTEM_LIBRARY})
link_directories(${Boost_LIBRARY_DIRS})

######################################################################
# Determine the supported system thread library and link to it
######################################################################
find_package(Threads REQUIRED)
if(CMAKE_USE_PTHREADS_INIT)
  message(STATUS "Adding -pthread as a compile argument")
  add_compile_options(-pthread)
endif()
if(CMAKE_THREAD_LIBS_INIT)
  message(STATUS "Adding ${CMAKE_THREAD_LIBS_INIT} as a linker argument")
  link_libraries(${CMAKE_THREAD_LIBS_INIT})
endif()


######################################################################
# Test for libjudy requirements
######################################################################
check_include_files(Judy.h JUDY_H_AVAILABLE)
include(CheckLibraryExists)
check_library_exists(Judy  Judy1FreeArray "" JUDY_LIB_AVAILABLE)
include(CheckTypeSize)
check_type_size(size_t SIZEOF_SIZE_T)
set(LIBJUDY_SUPPORT FALSE)
if(SIZEOF_SIZE_T LESS 8)
  message(STATUS "libJudy support disabled (no 64 bit support).")
elseif(JUDY_H_AVAILABLE AND JUDY_LIB_AVAILABLE) 
  message(STATUS "Enabling libJudy support.")
  set(LIBJUDY_SUPPORT TRUE)
  add_definitions(-DDYNAMO_JUDY)
  link_libraries(Judy)
else()
  message(STATUS "libJudy header/library missing.")
endif()


######################################################################
# Visualiser support
######################################################################
INCLUDE (CheckIncludeFiles)
set(OpenGL_GL_PREFERENCE "LEGACY")
find_package(OpenGL)
find_package(GLUT)
find_package(GLEW)
find_package(AVCodec)
find_package(PNG)

find_package(PkgConfig)
pkg_check_modules(GTKMM gtkmm-3.0)

set(VISUALIZER_SUPPORT TRUE CACHE BOOL "Attempt to enable the visualiser")
function(visualiser_dependency varname message)
  if(${varname})
    message(STATUS "${message} found.")
  else()
    message(STATUS "${message} not found - visualiser will not be built.")
    set(VISUALIZER_SUPPORT FALSE PARENT_SCOPE)
  endif()
endfunction(visualiser_dependency)

visualiser_dependency(OPENGL_FOUND "OpenGL libraries")
visualiser_dependency(GTKMM_FOUND "GTK/GTKMM 3 library")
visualiser_dependency(PNG_FOUND "libPNG library")
visualiser_dependency(GLUT_FOUND "GLUT/freeGLUT library")
visualiser_dependency(GLEW_FOUND "GLEW library")
#visualiser_dependency(GLX_FOUND "GLX library")

find_path(CAIROMM_INCLUDE_PATH "cairomm.h" PATH_SUFFIXES "cairomm" "cairomm-1.0/cairomm")
find_path(CAIROMM_CONFIG_PATH "cairommconfig.h" PATHS "/usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}" "/usr/lib64/${CMAKE_LIBRARY_ARCHITECTURE}" PATH_SUFFIXES "cairomm" "cairomm-1.0/cairomm" "cairomm-1.0/include" "lib/cairomm-1.0/include")

if(${CAIROMM_INCLUDE_PATH} STREQUAL "CAIROMM_INCLUDE_PATH-NOTFOUND" OR ${CAIROMM_CONFIG_PATH} STREQUAL "CAIROMM_CONFIG_PATH-NOTFOUND")
  message("-- Cairomm headers not found - visualiser will not be built")
  unset(CAIROMM_INCLUDE_PATH CACHE)
  unset(CAIROMM_CONFIG_PATH CACHE)
  set(VISUALIZER_SUPPORT FALSE)
else()
  message("-- Cairomm headers found")
endif()

if(VISUALIZER_SUPPORT)
  check_include_files(cwiid.h CWIID_H_AVAILABLE)
  check_library_exists(cwiid  cwiid_open "" CWIID_LIB_AVAILABLE)
  if(CWIID_H_AVAILABLE AND CWIID_LIB_AVAILABLE)
    message(STATUS "Enabling Wiimote head tracking support")
    add_definitions(-DCOIL_wiimote)
    link_libraries(cwiid)
  endif()

  find_package(OpenVR)
  if(OPENVR_FOUND)
    message("OpenVR found!")
    add_definitions(-DCOIL_OpenVR)
    include_directories(${OPENVR_INCLUDE_DIRS})
    link_libraries(${OPENVR_LIBRARY})
  else()
    message("OpenVR not found (try setting OPENVR_ROOT_DIR to the location of the openvr git repository")
  endif()


  add_definitions(${OPENGL_DEFINITONS} ${GLUT_DEFINITIONS} ${GLEW_DEFINITIONS} ${GTKMM_DEFINITIONS})
  include_directories(${GTKMM_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS} ${CAIROMM_INCLUDE_PATH} ${CAIROMM_CONFIG_PATH} ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
  link_libraries(${GTKMM_LIBRARIES} ${CAIROMM_LIBRARIES} ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLEW_LIBRARIES})
endif()

######################################################################
########## HEADER DIRECTORIES
######################################################################
include_directories(${PROJECT_SOURCE_DIR}/src/dynamo/)
include_directories(${PROJECT_SOURCE_DIR}/src/coil/)
include_directories(${PROJECT_SOURCE_DIR}/src/magnet)
include_directories(${PROJECT_SOURCE_DIR}/src/ext_include)

######################################################################
######### TARGETS
######################################################################
### Magnet library
function(magnet_test name) #Registers a magnet unit-test
  add_executable(magnet_${name}_exe ${CMAKE_CURRENT_SOURCE_DIR}/src/magnet/tests/${name}.cpp)
  target_link_libraries(magnet_${name}_exe)
  add_test(magnet_${name} magnet_${name}_exe)
endfunction(magnet_test)

#magnet_test(scan_test : tests/scan_test.cpp magnet /opencl//OpenCL)
#magnet_test(radixsort_NVIDIA_test : tests/radixsort_NVIDIA_test.cpp magnet /opencl//OpenCL scan_test ;
#magnet_test(radixsort_AMD_test : tests/radixsort_AMD_test.cpp magnet /opencl//OpenCL scan_test ;
#magnet_test(bitonicsort_test : tests/bitonicsort_test.cpp magnet /opencl//OpenCL ;
#magnet_test(heapsort_test : tests/heapsort_test.cpp magnet /opencl//OpenCL ;
#magnet_test(sorter_test : tests/sorter_test.cpp magnet /opencl//OpenCL ;

magnet_test(threadpool_test)
#SET_TARGET_PROPERTIES(magnet_threadpool_test_exe PROPERTIES LINK_FLAGS -Wl,--no-as-needed) #Fix for a bug in gcc

target_link_libraries(magnet_threadpool_test_exe ${CMAKE_THREAD_LIBS_INIT})
magnet_test(cubic_quartic_test)
magnet_test(vector_test)
magnet_test(quaternion_test)
magnet_test(dilate_test)
magnet_test(splinetest)
magnet_test(plane_intersection)
magnet_test(triangle_intersection)
magnet_test(intersection_genalg)
magnet_test(offcenterspheres)
magnet_test(stack_vector_test)

if(JUDY_SUPPORT)
  magnet_test(judy_test)
endif(JUDY_SUPPORT)

if(AVCodec_FOUND)
  set(CMAKE_REQUIRED_INCLUDES ${AVCodec_INCLUDE_DIRS})
  check_cxx_source_compiles("#include <libavcodec/avcodec.h>
int main() { &avcodec_send_frame; return 0; }" AVCodec_Send_Frame_Avail)
  if(${AVCodec_Send_Frame_Avail})
    message(STATUS "Found suitable libAVCodec version")
    link_libraries(${AVCodec_LIBRARIES})
    include_directories(${AVCodec_INCLUDE_DIRS})
    add_definitions(-DMAGNET_FFMPEG_SUPPORT)
  else()
    message(STATUS "libAVCodec version is out of date or incompatible - recording will be disabled")
  endif()
endif()

### Coil
# configuration
if(VISUALIZER_SUPPORT)
  message(STATUS "Visualiser dependencies satisfied.")
  
  # Coil library
  file(GLOB_RECURSE coil_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/coil/coil/*.cpp)
  add_library(coil STATIC ${coil_SRC})
  link_libraries(coil)
  add_definitions(-DDYNAMO_visualizer)
  add_executable(coilvol ${CMAKE_CURRENT_SOURCE_DIR}/src/coil/programs/coilvol.cpp)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/coilvol DESTINATION bin)
endif()

### DynamO
file(GLOB_RECURSE dynamo_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/dynamo/*.cpp)
add_library(dynamo STATIC ${dynamo_SRC})
link_libraries(dynamo)

function(dynamo_exe name) #Registers a dynamo executable given the source file name
  add_executable(${name} ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/${name}.cpp)
  install(PROGRAMS $<TARGET_FILE:${name}> DESTINATION bin)
endfunction(dynamo_exe)
dynamo_exe(dynarun)
dynamo_exe(dynamod)
dynamo_exe(dynahist_rw)
dynamo_exe(dynapotential)
#dynamo_exe(dynacollide)
if(VISUALIZER_SUPPORT)
  #Can't use dynamo_exe here, as we just need to compile "dynarun.cpp" differently
  add_executable(dynavis ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarun.cpp)
  set_target_properties(dynavis PROPERTIES COMPILE_DEFINITIONS "DYNAMO_loadvisualiser")
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/dynavis DESTINATION bin)
endif()

# scripts
if(Python3_FOUND)
  if(Python3_NumPy_FOUND)
    install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynatransport DESTINATION bin)
    install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarmsd      DESTINATION bin)
    message(STATUS "Python NumPy found, installing dynatransport, dynamaprmsd, and dynarmsd")
  else()
    message(WARNING "Python NumPy library not found, disabling install of dynatransport, dynamaprmsd, and dynarmsd")
  endif()

  install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynamo2xyz    DESTINATION bin)
  install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynamaprmsd   DESTINATION bin)
else()
  message(WARNING "Python 3 not found, cannot install all DynamO tools.")
endif()


# unit tests
function(dynamo_test name) #Registers a unit test of DynamO
  add_executable(dynamo_${name}_exe ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/tests/${name}.cpp)
  add_test(dynamo_${name} dynamo_${name}_exe)
endfunction(dynamo_test)

dynamo_test(hardsphere_test)
dynamo_test(shearing_test)
dynamo_test(binaryhardsphere_test)
dynamo_test(squarewell_test)
dynamo_test(2dstepped_potential_test)
dynamo_test(infmass_spheres_test)
dynamo_test(lines_test)
dynamo_test(static_spheres_test)
dynamo_test(gravityplate_test)
dynamo_test(polymer_test)
dynamo_test(swingspheres_test)
dynamo_test(squarewellwall_test)
dynamo_test(thermalisedwalls_test)
dynamo_test(event_sorters_test)


if(Python3_Interpreter_FOUND)
  add_test(NAME dynamo_replica_exchange
    COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/tests/replex_test.py
    --dynarun=$<TARGET_FILE:dynarun>
    --dynamod=$<TARGET_FILE:dynamod>
    --dynahist_rw=$<TARGET_FILE:dynahist_rw>)

  add_test(NAME dynamo_multicanonical_cmap
    COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/tests/multicanonical_cmap_test.py
    --dynarun=$<TARGET_FILE:dynarun>
    --dynamod=$<TARGET_FILE:dynamod>
    --dynahist_rw=$<TARGET_FILE:dynahist_rw>)
  
  if(Python3_NumPy_FOUND)
    add_test(NAME dynamo_dynatransport
      COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/tests/dynatransport_test.py
      --dynarun=$<TARGET_FILE:dynarun>
      --dynamod=$<TARGET_FILE:dynamod>
      --dynahist_rw=$<TARGET_FILE:dynahist_rw>
      --dynatransport=${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynatransport
      --python=${Python3_EXECUTABLE}
      )
  endif()
endif()
