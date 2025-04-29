cmake_minimum_required(VERSION 3.28)

if(POLICY CMP0167) # Must be before project otherwise VCPKG hooks trigger the warning
  cmake_policy(SET CMP0167 NEW) #  Use the new boost configuration
endif()


if(NOT MSVC)
  # To make this compile on NixOS, we need to use the static linked MSVC
  # But this is not available on windows
  set(VCPKG_BOOTSTRAP_OPTIONS "-musl")
endif()
set(CMAKE_TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/extern/vcpkg/scripts/buildsystems/vcpkg.cmake)

project(DynamO) #Project name
enable_testing() #Enable build of test executables and 'make test' command
include(CTest)

# Include our custom CMAKE modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Define any build options
option(VISUALIZER_ENABLED "Enable the visualizer" ON)

# if SKBUILD_SCRIPTS_DIR is not set, set it to a non-absolute path (otherwise windows builds fail)
if(NOT DEFINED SKBUILD_SCRIPTS_DIR)
  set(SKBUILD_SCRIPTS_DIR "bin")
endif()

######################################################################
########## COMPILER FLAGS
######################################################################
# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Check the build type has been configured
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Fix some compiler flags for MSVC
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
  link_libraries(debug dl)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # GCC incorrectly determines array accesses are out of bounds
  # when using std::sort on std::vector. This is a known issue
  # with GCC and std::sort. See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=107986
  add_compile_options(-Wno-array-bounds)
endif()

# Set some compiler flags for debug mode
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDYNAMO_DEBUG")

######################################################################
########## Git information and check submodules are installed
######################################################################
find_package(Git QUIET)

include(GetGitRevisionDescription)
git_describe(GIT_RESULT --match initial-commit)
message(STATUS "Git description is ${GIT_RESULT}")
string(REGEX MATCH "[0-9]+" BUILD_NUMBER ${GIT_RESULT})
message(STATUS "Build number is ${BUILD_NUMBER}")

########## Check we have git submodules installed!
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/stator/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

######################################################################
######### Dependencies
######################################################################

#set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.70.0 REQUIRED system filesystem program_options unit_test_framework)

add_subdirectory(src/magnet)
find_package(Eigen3)
add_subdirectory(src/coil)

######################################################################
### DynamO
######################################################################

file(GLOB_RECURSE dynamo_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/dynamo/*.cpp)
add_library(dynamo STATIC ${dynamo_SRC})
target_link_libraries(dynamo PUBLIC magnet Boost::program_options Boost::system Boost::filesystem Eigen3::Eigen)
target_include_directories(dynamo PUBLIC ${PROJECT_SOURCE_DIR}/src/dynamo/)

message(STATUS "Coil_FOUND: ${Coil_FOUND}")
if (Coil_FOUND)
  # This needs to be PUBLIC as anything using the headers (i.e. unit tests)
  # needs to know about coil being enabled, otherwise we get memory layouts of
  # locals being different in different translation units
  target_compile_definitions(dynamo PUBLIC DYNAMO_visualizer)
  target_link_libraries(dynamo PUBLIC coil)
  # Add the dynavis program
  add_executable(dynavis ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarun.cpp)
  # This define is only used in dynarun.cpp
  target_compile_definitions(dynavis PRIVATE DYNAMO_loadvisualizer)
  target_link_libraries(dynavis dynamo)
  install(TARGETS dynavis DESTINATION "${SKBUILD_SCRIPTS_DIR}")
endif()

# Test for libjudy requirements
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
  target_compile_definitions(dynamo PRIVATE DYNAMO_JUDY)
  target_link_libraries(dynamo PRIVATE Judy)
  magnet_test(judy_test)
else()
  message(STATUS "libJudy header/library missing.")
endif()

# Test for Bzip2 requirements
find_package(BZip2)
if(BZIP2_FOUND)
  target_include_directories(dynamo PRIVATE ${BZIP2_INCLUDE_DIR})
  target_link_libraries(dynamo PRIVATE ${BZIP2_LIBRARIES})
  target_compile_definitions(dynamo PUBLIC DYNAMO_bzip2_support=1)
  message(STATUS "libbz2 found - compressed file support enabled")
else()
  message(WARNING "libbz2 not found - compressed file support disabled")
endif()

function(dynamo_exe name) #Registers a dynamo executable given the source file name
  add_executable(${name} ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/${name}.cpp)
  target_link_libraries(${name} dynamo)
  install(TARGETS ${name} DESTINATION "${SKBUILD_SCRIPTS_DIR}")
endfunction(dynamo_exe)

dynamo_exe(dynarun)
dynamo_exe(dynamod)
dynamo_exe(dynahist_rw)
dynamo_exe(dynapotential)

# scripts
find_package(Python3 COMPONENTS Interpreter NumPy)
if(Python3_FOUND)
  if(Python3_NumPy_FOUND)
    install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynatransport DESTINATION "${SKBUILD_SCRIPTS_DIR}")
    install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynarmsd      DESTINATION "${SKBUILD_SCRIPTS_DIR}")
    message(STATUS "Python NumPy found, installing dynatransport, dynamaprmsd, and dynarmsd")
  else()
    message(WARNING "Python NumPy library not found, disabling install of dynatransport, dynamaprmsd, and dynarmsd")
  endif()

  install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynamo2xyz    DESTINATION "${SKBUILD_SCRIPTS_DIR}")
  install(PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/programs/dynamaprmsd   DESTINATION "${SKBUILD_SCRIPTS_DIR}")
else()
  message(WARNING "Python 3 not found, cannot install all DynamO tools.")
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

#######################################################################
# Checking and applying formatting to source files
#######################################################################
# Check if clang-format is available
find_program(CLANG_FORMAT_EXECUTABLE NAMES clang-format)
if(CLANG_FORMAT_EXECUTABLE)
  message(STATUS "clang-format found: ${CLANG_FORMAT_EXECUTABLE}")
else()
  message(WARNING "clang-format not found. Skipping formatting checks.")
endif()

file(GLOB_RECURSE dynamo_SRC_hpp ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/*.hpp)
file(GLOB_RECURSE dynamo_SRC_cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/dynamo/*.cpp)
file(GLOB_RECURSE coil_SRC_hpp ${CMAKE_CURRENT_SOURCE_DIR}/src/coil/*.hpp)
file(GLOB_RECURSE coil_SRC_cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/coil/*.cpp)
file(GLOB_RECURSE magnet_SRC_hpp ${CMAKE_CURRENT_SOURCE_DIR}/src/magnet/*.hpp)
file(GLOB_RECURSE magnet_SRC_cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/magnet/*.cpp)

set(
  CPP_SOURCE_FILES
  ${dynamo_SRC_hpp} ${dynamo_SRC_cpp}
  ${coil_SRC_hpp} ${coil_SRC_cpp}
  ${magnet_SRC_hpp} 
  ${magnet_SRC_cpp}
)

add_custom_target(
    check-src-format
    COMMAND ${CLANG_FORMAT_EXECUTABLE} --dry-run ${CPP_SOURCE_FILES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "checking format of source files..."
)

add_custom_target(
    src-format
    COMMAND clang-format -i ${CPP_SOURCE_FILES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Formatting all source files in src..."
)

########################################################################
########## CPack
########################################################################
if(NOT DEFINED CPACK_GENERATOR)
  set(CPACK_GENERATOR ${SPECIFIC_SYSTEM_PREFERED_CPACK_GENERATOR})
endif()
set(PACKAGE_SUFFIX "" CACHE STRING "A suffix to add to the end of the built packages")
set(CPACK_PACKAGE_NAME "dynamomd")
set(CPACK_PACKAGE_VERSION_MAJOR 1)
set(CPACK_PACKAGE_VERSION_MINOR 7)
set(CPACK_PACKAGE_VERSION_PATCH ${BUILD_NUMBER})
set(CPACK_PACKAGE_CONTACT "Marcus Bannerman <support@dynamomd.com>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Event-driven particle simulation sofware")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_SOURCE_DIR}/description.txt)
set(CPACK_PACKAGE_DESCRIPTION "This package includes the DynamO executables and scripts. DynamO is an event driven simulator capable of performing molecular and granular dynamics simulations")

if(Python3_FOUND)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS " python3")
endif()

set(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}${PACKAGE_SUFFIX})
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "http://dynamomd.com")
#Get Debian packages to automatically figure out the dependencies
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_WIX_UPGRADE_GUID "c70e546c-a36a-429f-a2df-f5e6cf787bb4")
set(CPACK_WIX_PRODUCT_GUID "36fe5a30-bd03-423b-822a-f505d0fa6564")
include(CPack) #Enable build of packages
install(FILES ${CMAKE_SOURCE_DIR}/copyright DESTINATION share/doc/${CPACK_PACKAGE_NAME})
