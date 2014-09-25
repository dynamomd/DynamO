cmake_minimum_required (VERSION 2.6 FATAL_ERROR)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CheckTypeSize)
include(CheckCXXCompilerFlag)

##########   OPTIONS
enable_testing() #Enable build of test executables and 'make test' command
project (DynamO) #Project name

########## SUPPORT TESTING
#C++11/C++0x support
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

#Boost support
find_package(Boost 1.46 REQUIRED COMPONENTS "program_options" "filesystem" "iostreams" "system" "unit_test_framework")
if(NOT Boost_FOUND)
  message(SEND_ERROR "Cannot build any programs, boost library is missing.")
endif()
include_directories(${Boost_INCLUDE_DIRS})

find_package(PkgConfig)
find_package(OpenGL)
find_package(GLUT)
find_package(GLEW)
pkg_check_modules(GTKMM gtkmm-2.4)
pkg_check_modules(PNG libpng)

find_package(PythonInterp 2.7)

check_include_files(Judy.h JUDY_H_AVAILABLE)
check_library_exists(Judy  Judy1FreeArray "" JUDY_LIB_AVAILABLE)
CHECK_TYPE_SIZE(size_t SIZEOF_SIZE_T)

find_package( Threads )
########## HEADER DIRECTORIES
include_directories(${PROJECT_SOURCE_DIR}/src/dynamo/)
include_directories(${PROJECT_SOURCE_DIR}/src/coil/)
include_directories(${PROJECT_SOURCE_DIR}/src/magnet)
include_directories(${PROJECT_SOURCE_DIR}/src/ext_include)

######### TARGETS
### Magnet library
function(magnet_test name) #Registers a magnet unit-test
  add_executable(magnet_${name}_exe ${CMAKE_CURRENT_SOURCE_DIR}/src/magnet/tests/${name}.cpp)
  target_link_libraries(magnet_${name}_exe ${Boost_LIBRARIES})
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
if((JUDY_H_AVAILABLE AND JUDY_LIB_AVAILABLE) AND (SIZEOF_SIZE_T LESS 8))
  magnet_test(judy_test)
  target_link_libraries(magnet_judy_test_exe Judy)
endif()
magnet_test(plane_intersection)
magnet_test(triangle_intersection)
magnet_test(intersection_genalg)
magnet_test(offcenterspheres)

### DynamO
file(GLOB_RECURSE dynamo_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/dynamo/*.cpp)
add_library(dynamo STATIC ${dynamo_SRC})
target_link_libraries(dynamo ${Boost_LIBRARIES})

if(SIZEOF_SIZE_T LESS 8)
  message(STATUS "libJudy support disabled (no 64 bit support).")
elseif(JUDY_H_AVAILABLE AND JUDY_LIB_AVAILABLE) 
  message(STATUS "Enabling libJudy support.")
  add_definitions(-DYNAMO_JUDY)
  target_link_libraries(dynamo Judy)
else()
  message(STATUS "libJudy header/library missing.")
endif()

# executables
function(dynamo_exe name) #Registers a dynamo executable given the source file name
  add_executable(${name} ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/${name}.cpp)
  target_link_libraries(${name} dynamo)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${name} DESTINATION bin)
endfunction(dynamo_exe)

dynamo_exe(dynarun)
dynamo_exe(dynamod)
dynamo_exe(dynahist_rw)
#dynamo_exe(dynacollide)
#dynamo_exe(dynapotential)

# scripts
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
  target_link_libraries(dynamo_${name}_exe dynamo)
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

### Coil
# configuration
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
	    message(STATUS "Visualiser dependencies satisfied. Enabling dynavis build.")

	    # further configuration
	    include_directories(${GTKMM_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
	    add_definitions(${OPENGL_DEFINITONS} ${GLUT_DEFINITIONS} ${GLEW_DEFINITIONS})
	    # Coil library
	    set_target_properties(dynamo PROPERTIES COMPILE_DEFINITIONS "DYNAMO_visualizer")
	    file(GLOB_RECURSE coil_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/coil/coil/*.cpp)
	    add_library(coil STATIC ${coil_SRC})
	    target_link_libraries(coil ${GTKMM_LIBRARIES} ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLEW_LIBRARIES})

	    # Add coil as dependency to Dynamo
	    target_link_libraries(dynamo coil)
	    add_definitions(-DDYNAMO_visualiser)

	    # dynavis
	    add_executable(dynavis ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarun.cpp)
	    target_link_libraries(dynavis dynamo)  
	  endif()    
	endif()    
      endif()    
    endif()    
  endif()    
endif()    
