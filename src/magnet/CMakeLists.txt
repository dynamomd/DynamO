set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMAGNET_DEBUG")

find_package(Boost 1.70.0 REQUIRED unit_test_framework)


add_library(magnet INTERFACE)
target_include_directories(magnet INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${CMAKE_CURRENT_SOURCE_DIR}/../../extern
)

######################################################################
# Stack trace support
######################################################################
include(CheckIncludeFiles)
check_include_files(execinfo.h STACKTRACE_AVAILABLE)
if(STACKTRACE_AVAILABLE)
  target_compile_definitions(magnet INTERFACE MAGNET_STACKTRACE)
endif()

######################################################################
# Determine the supported system thread library and link to it when needed
######################################################################
find_package(Threads REQUIRED)

if(CMAKE_USE_PTHREADS_INIT)
  message(STATUS "Adding -pthread as a compile argument")
  target_compile_options(magnet INTERFACE -pthread)
endif()

if(CMAKE_THREAD_LIBS_INIT)
  message(STATUS "Adding ${CMAKE_THREAD_LIBS_INIT} as a linker argument")
  target_link_libraries(magnet INTERFACE ${CMAKE_THREAD_LIBS_INIT})
endif()

function(magnet_test name) #Registers a magnet unit-test
  add_executable(magnet_${name}_exe ${CMAKE_CURRENT_SOURCE_DIR}/tests/${name}.cpp)
  target_link_libraries(magnet_${name}_exe magnet)
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
