######################################################################
# Visualizer support
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

set(VISUALIZER_SUPPORT TRUE CACHE BOOL "Visualizer build supported")
function(visualizer_dependency varname message)
  if(${varname})
    message(STATUS "${message} found.")
  else()
    message(STATUS "${message} not found - visualizer cannot be built.")
    set(VISUALIZER_SUPPORT FALSE PARENT_SCOPE)
  endif()
endfunction(visualizer_dependency)

visualizer_dependency(OPENGL_FOUND "OpenGL libraries")
visualizer_dependency(GTKMM_FOUND "GTK/GTKMM 3 library")
visualizer_dependency(PNG_FOUND "libPNG library")
visualizer_dependency(GLUT_FOUND "GLUT/freeGLUT library")
visualizer_dependency(GLEW_FOUND "GLEW library")
#visualizer_dependency(GLX_FOUND "GLX library")

find_path(CAIROMM_INCLUDE_PATH "cairomm.h" PATH_SUFFIXES "cairomm" "cairomm-1.0/cairomm")
find_path(CAIROMM_CONFIG_PATH "cairommconfig.h" PATHS "/usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}" "/usr/lib64/${CMAKE_LIBRARY_ARCHITECTURE}" PATH_SUFFIXES "cairomm" "cairomm-1.0/cairomm" "cairomm-1.0/include" "lib/cairomm-1.0/include")

if(${CAIROMM_INCLUDE_PATH} STREQUAL "CAIROMM_INCLUDE_PATH-NOTFOUND" OR ${CAIROMM_CONFIG_PATH} STREQUAL "CAIROMM_CONFIG_PATH-NOTFOUND")
  message("-- Cairomm headers not found - visualizer cannot be built")
  unset(CAIROMM_INCLUDE_PATH CACHE)
  unset(CAIROMM_CONFIG_PATH CACHE)
  set(VISUALIZER_SUPPORT FALSE)
else()
  message("-- Cairomm headers found")
endif()

##########################################################################
### Coil
##########################################################################
set(Coil_FOUND (VISUALIZER_SUPPORT AND VISUALIZER_ENABLED) CACHE BOOL "Coil library available")

if(VISUALIZER_SUPPORT AND VISUALIZER_ENABLED)
  message(STATUS "Visualizer dependencies satisfied.")

  # Coil library
  file(GLOB_RECURSE coil_SRC ${CMAKE_CURRENT_SOURCE_DIR}/coil/*.cpp)
  add_library(coil STATIC ${coil_SRC})
  target_link_libraries(coil PUBLIC magnet Boost::program_options ${GTKMM_LIBRARIES} ${CAIROMM_LIBRARIES} ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLEW_LIBRARIES})
  target_include_directories(coil PUBLIC ${PROJECT_SOURCE_DIR}/src/coil/)
  target_include_directories(coil PUBLIC ${PROJECT_SOURCE_DIR}/extern/stator/)
  target_include_directories(coil PUBLIC ${GTKMM_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS} ${CAIROMM_INCLUDE_PATH} ${CAIROMM_CONFIG_PATH} ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS})
  target_compile_definitions(coil PUBLIC ${OPENGL_DEFINITONS} ${GLUT_DEFINITIONS} ${GLEW_DEFINITIONS} ${GTKMM_DEFINITIONS})

  # cwiid support
  check_include_files(cwiid.h CWIID_H_AVAILABLE)
  check_library_exists(cwiid  cwiid_open "" CWIID_LIB_AVAILABLE)
  if(CWIID_H_AVAILABLE AND CWIID_LIB_AVAILABLE)
    message(STATUS "Enabling Wiimote head tracking support")
    target_compile_definitions(coil PRIVATE COIL_wiimote)
    target_link_libraries(coil cwiid)
  endif()

  # OpenVR support
  find_package(OpenVR)
  if(OPENVR_FOUND)
    message("OpenVR found!")
    target_add_definitions(coil PRIVATE COIL_OpenVR)
    target_include_directories(coil ${OPENVR_INCLUDE_DIRS})
    target_link_libraries(coil ${OPENVR_LIBRARY})
  else()
    message("OpenVR not found (try setting OPENVR_ROOT_DIR to the location of the openvr git repository")
  endif()
  
  # libAVCodec support
  if(AVCodec_FOUND)
    set(CMAKE_REQUIRED_INCLUDES ${AVCodec_INCLUDE_DIRS})
    include(CheckCXXCompilerFlag)
    check_cxx_source_compiles("#include <libavcodec/avcodec.h>
  int main() { &avcodec_send_frame; return 0; }" AVCodec_Send_Frame_Avail)
    if(${AVCodec_Send_Frame_Avail})
      message(STATUS "Found suitable libAVCodec version")
      target_link_libraries(coil PRIVATE ${AVCodec_LIBRARIES})
      target_include_directories(coil PRIVATE ${AVCodec_INCLUDE_DIRS})
      target_compile_definitions(coil PRIVATE MAGNET_FFMPEG_SUPPORT)
    else()
      message(STATUS "libAVCodec version is out of date or incompatible - recording will be disabled")
    endif()
  endif()

  # Add coilvol program
  add_executable(coilvol ${CMAKE_CURRENT_SOURCE_DIR}/programs/coilvol.cpp)
  target_link_libraries(coilvol PUBLIC coil Boost::program_options Boost::filesystem Boost::system)
  install(TARGETS coilvol RUNTIME DESTINATION bin)
endif()
