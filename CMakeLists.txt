cmake_minimum_required (VERSION 2.6 FATAL_ERROR)
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()

##########   OPTIONS
project(DynamO) #Project name
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

enable_testing() #Enable build of test executables and 'make test' command
include(CTest)

set(CPACK_GENERATOR DEB RPM)
set(CPACK_PACKAGE_NAME "dynamo")
set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 6)
set(CPACK_PACKAGE_VERSION_PATCH 0)
set(CPACK_PACKAGE_CONTACT "support@dynamomd.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "An event-driven particle simulation tool.")
include(CPack) #Enable build of packages

##########   RELEASE MODE
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

##########   DEBUG MODE
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDYNAMO_DEBUG -DMAGNET_DEBUG -Wall -g")

link_libraries(debug dl)

########## COMPILER SUPPORT TESTS
#C++11/C++0x support
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORT_CXX11)
if(COMPILER_SUPPORT_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORT_CXX0X)
  if(COMPILER_SUPPORT_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    check_cxx_source_compiles("int main(int, const char**) {int array[5] = {1,2,3,4,5}; for(int& x: array) x *= 2;}" CXX0X_RANGE_BASED_FOR)
    if(NOT CXX0X_RANGE_BASED_FOR)
      message(SEND_ERROR "The compiler ${CMAKE_CXX_COMPILER} has limited C++11 support. You must install a more modern C++ compiler.")
    endif()
  else()
    message(SEND_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. You must install a more modern C++ compiler.")
  endif()
endif()

check_cxx_source_compiles("template <typename T> struct wrapped{ typedef T type;};
template <typename T> typename T::type unwrap1(T a) { return typename T::type();}
int unwrap(int a) { return 1; }
template <typename T> auto unwrap(T t) -> decltype(unwrap(unwrap1(t))) { return unwrap(unwrap1(t)); }
int main() { unwrap(wrapped<wrapped<int>>()); }" CXX11_TYPE_MANGLING)
if(NOT CXX11_TYPE_MANGLING)
  message(WARNING "Although the compiler advertises support for C++0x/C++11, this support is limited (gcc-4.7 or lower). Many advanced potentials will be disabled.")
endif()


# Boost support
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.46 COMPONENTS "program_options" "filesystem" "iostreams" "system" "unit_test_framework")
if(NOT Boost_FOUND)
  message(STATUS "Cannot find a static installation of boost, trying a dynamic version")
  set(Boost_USE_STATIC_LIBS OFF)
  find_package(Boost 1.46 REQUIRED COMPONENTS "program_options" "filesystem" "iostreams" "system" "unit_test_framework")
  if(NOT Boost_FOUND)
    message(SEND_ERROR "Cannot find a static installation of boost, trying a dynamic version")
  endif()
endif()
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_FILESYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_IOSTREAMS_LIBRARY} ${Boost_SYSTEM_LIBRARY})

#BZip2 support required for boost iostreams
find_package(BZip2 REQUIRED)
include_directories(${BZIP2_INCLUDE_DIR})
link_libraries(${BZIP2_LIBRARIES})

# Determine the supported system thread library and link to it
find_package(Threads)
link_libraries(${CMAKE_THREAD_LIBS_INIT})

# Test for libjudy requirements
include(CheckIncludeFiles)
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

# Visualiser support
find_package(PkgConfig)
find_package(OpenGL)
find_package(GLUT)
find_package(GLEW)
pkg_check_modules(GTKMM gtkmm-2.4)
pkg_check_modules(PNG libpng)
set(VISUALIZER_SUPPORT FALSE)
if(OPENGL_FOUND)
  message(STATUS "OpenGL found.")
  if(GTKMM_FOUND)
    message(STATUS "GTKMM found.")
    if(PNG_FOUND)
      message(STATUS "libPNG found.")
      if(GTKMM_FOUND)
	message(STATUS "GTKmm found.")
	if(GLUT_FOUND)
	  message(STATUS "GLUT found.")
	  if(GLEW_FOUND)
	    message(STATUS "GLEW found.")
	    include_directories(${GTKMM_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
	    add_definitions(${OPENGL_DEFINITONS} ${GLUT_DEFINITIONS} ${GLEW_DEFINITIONS})
	    link_libraries(${GTKMM_LIBRARIES} ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLEW_LIBRARIES})
	    set(VISUALIZER_SUPPORT TRUE)
	  endif()
	endif()
      endif()
    endif()
  endif()
endif()

########## HEADER DIRECTORIES
include_directories(${PROJECT_SOURCE_DIR}/src/dynamo/)
include_directories(${PROJECT_SOURCE_DIR}/src/coil/)
include_directories(${PROJECT_SOURCE_DIR}/src/magnet)
include_directories(${PROJECT_SOURCE_DIR}/src/ext_include)

######### TARGETS
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
SET_TARGET_PROPERTIES(magnet_threadpool_test_exe PROPERTIES LINK_FLAGS -Wl,--no-as-needed) #Fix for a bug in gcc

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

#### CAS library tests ####
if(CXX11_TYPE_MANGLING)
  magnet_test(polynomial_test)
  magnet_test(symbolic_test)
  magnet_test(stable_alg_test)
  magnet_test(numeric_test)
  magnet_test(symbolic_performance)  
endif()

if(JUDY_SUPPORT)
  magnet_test(judy_test)
endif(JUDY_SUPPORT)

#Test for AV
find_package(AVCodec)
if(AVCodec_FOUND)
  #Note, avcodec has multiple headers, perhaps set up a config.h file?
  link_libraries(${AVCodec_LIBRARIES})
  include_directories(${AVCodec_INCLUDE_DIRS})
  add_definitions(-DMAGNET_FFMPEG_SUPPORT)
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
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${name} DESTINATION bin)
endfunction(dynamo_exe)
dynamo_exe(dynarun)
dynamo_exe(dynamod)
dynamo_exe(dynahist_rw)
#dynamo_exe(dynacollide)
#dynamo_exe(dynapotential)
if(VISUALIZER_SUPPORT)
  add_executable(dynavis ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarun.cpp)
  set_target_properties(dynavis PROPERTIES COMPILE_DEFINITIONS "DYNAMO_loadvisualiser")
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/dynavis DESTINATION bin)
endif()

# scripts
find_package(PythonInterp 2.7)
if(PYTHONINTERP_FOUND)
  install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynatransport DESTINATION bin)
  install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynamo2xyz    DESTINATION bin)
  install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarmsd      DESTINATION bin)
  install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynamaprmsd   DESTINATION bin)
else()
  message(WARNING "Python 2.7 not found, cannot install all DynamO tools.")
endif()

# unit tests
function(dynamo_test name) #Registers a unit test of DynamO
  add_executable(dynamo_${name}_exe ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/tests/${name}.cpp)
  set_target_properties(dynamo_${name}_exe PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK")
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
