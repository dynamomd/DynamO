cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

##########   OPTIONS
enable_testing() #Enable build of test executables and 'make test' command
project(DynamO) #Project name

##########   RELEASE MODE
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

##########   DEBUG MODE
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDYNAMO_DEBUG -DMAGNET_DEBUG -Wall -fstack-protector-all")
link_libraries(debug dl)

########## COMPILER SUPPORT TESTS
#C++11/C++0x support
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(SEND_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++0x/C++11 support. Please use a more modern C++ compiler.")
  endif()
endif()

# Boost support
find_package(Boost 1.46 REQUIRED COMPONENTS "program_options" "filesystem" "iostreams" "system" "unit_test_framework")
if(NOT Boost_FOUND)
  message(SEND_ERROR "Cannot build DynamO as the required boost libraries are missing.")
endif()
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_FILESYSTEM_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_IOSTREAMS_LIBRARY} ${Boost_SYSTEM_LIBRARY})

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
  add_definitions(-DYNAMO_JUDY)
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
  target_link_libraries(magnet_${name}_exe ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  set_target_properties(magnet_${name}_exe PROPERTIES COMPILE_DEFINITIONS "BOOST_TEST_DYN_LINK")
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
if(JUDY_SUPPORT)
  magnet_test(judy_test)
endif(JUDY_SUPPORT)

### Coil
# configuration
if(VISUALIZER_SUPPORT)
  message(STATUS "Visualiser dependencies satisfied.")
  # Coil library
  file(GLOB_RECURSE coil_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/coil/coil/*.cpp)
  add_library(coil STATIC ${coil_SRC})
  link_libraries(coil)
  add_definitions(-DDYNAMO_visualiser)
endif()

### DynamO
file(GLOB_RECURSE dynamo_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/dynamo/*.cpp)
add_library(dynamo OBJECT ${dynamo_SRC})

function(dynamo_exe name) #Registers a dynamo executable given the source file name
  add_executable(${name} ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/${name}.cpp $<TARGET_OBJECTS:dynamo>)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${name} DESTINATION bin)
endfunction(dynamo_exe)
dynamo_exe(dynarun)
dynamo_exe(dynamod)
dynamo_exe(dynahist_rw)
#dynamo_exe(dynacollide)
#dynamo_exe(dynapotential)
if(VISUALIZER_SUPPORT)
  add_executable(dynavis ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarun.cpp $<TARGET_OBJECTS:dynamo>)
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
  add_executable(dynamo_${name}_exe ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/tests/${name}.cpp $<TARGET_OBJECTS:dynamo>)
  target_link_libraries(dynamo_${name}_exe ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
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
