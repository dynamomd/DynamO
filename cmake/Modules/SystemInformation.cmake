# define a set of string with may-be useful readable name
# this file is meant to be included in a CMakeLists.txt
# not as a standalone CMake script
set(SPECIFIC_COMPILER_NAME "")
set(SPECIFIC_SYSTEM_VERSION_NAME "")
set(SPECIFIC_SYSTEM_PREFERED_CPACK_GENERATOR "")
set(DISTRO_ID "")
set(DISTRO_RELEASE "")

function (read_release valuename FROM filename INTO varname)
  file (STRINGS ${filename} _distrib REGEX "^${valuename}=")
  string (REGEX REPLACE "^${valuename}=\"?\(.*\)" "\\1" ${varname} "${_distrib}")
  # remove trailing quote that got globbed by the wildcard (greedy match)
  string (REGEX REPLACE "\"$" "" ${varname} "${${varname}}")
  set (${varname} "${${varname}}" PARENT_SCOPE)
endfunction (read_release valuename FROM filename INTO varname)

# In the WIN32 case try to guess a "readable system name"
if(WIN32)
  set(SPECIFIC_SYSTEM_PREFERED_PACKAGE "NSIS")
  # information taken from
  # http://www.codeguru.com/cpp/w-p/system/systeminformation/article.php/c8973/
  # Win9x series
  if(CMAKE_SYSTEM_VERSION MATCHES "4.0")
    set(SPECIFIC_SYSTEM_VERSION_NAME "Win95")
  endif(CMAKE_SYSTEM_VERSION MATCHES "4.0")
  if(CMAKE_SYSTEM_VERSION MATCHES "4.10")
    set(SPECIFIC_SYSTEM_VERSION_NAME "Win98")
  endif(CMAKE_SYSTEM_VERSION MATCHES "4.10")
  if(CMAKE_SYSTEM_VERSION MATCHES "4.90")
    set(SPECIFIC_SYSTEM_VERSION_NAME "WinME")
  endif(CMAKE_SYSTEM_VERSION MATCHES "4.90")

  # WinNTyyy series
  if(CMAKE_SYSTEM_VERSION MATCHES "3.0")
    set(SPECIFIC_SYSTEM_VERSION_NAME "WinNT351")
  endif(CMAKE_SYSTEM_VERSION MATCHES "3.0")
  if(CMAKE_SYSTEM_VERSION MATCHES "4.1")
    set(SPECIFIC_SYSTEM_VERSION_NAME "WinNT4")
  endif(CMAKE_SYSTEM_VERSION MATCHES "4.1")

  # Win2000/XP series
  if(CMAKE_SYSTEM_VERSION MATCHES "5.0")
    set(SPECIFIC_SYSTEM_VERSION_NAME "Win2000")
  endif(CMAKE_SYSTEM_VERSION MATCHES "5.0")
  if(CMAKE_SYSTEM_VERSION MATCHES "5.1")
    set(SPECIFIC_SYSTEM_VERSION_NAME "WinXP")
  endif(CMAKE_SYSTEM_VERSION MATCHES "5.1")
  if(CMAKE_SYSTEM_VERSION MATCHES "5.2")
    set(SPECIFIC_SYSTEM_VERSION_NAME "Win2003")
  endif(CMAKE_SYSTEM_VERSION MATCHES "5.2")

  # WinVista/7/8 series
  if(CMAKE_SYSTEM_VERSION MATCHES "6.0")
    set(SPECIFIC_SYSTEM_VERSION_NAME "WinVISTA")
  endif(CMAKE_SYSTEM_VERSION MATCHES "6.0")
  if(CMAKE_SYSTEM_VERSION MATCHES "6.1")
    set(SPECIFIC_SYSTEM_VERSION_NAME "Win7")
  endif(CMAKE_SYSTEM_VERSION MATCHES "6.1")
  if(CMAKE_SYSTEM_VERSION MATCHES "6.2")
    set(SPECIFIC_SYSTEM_VERSION_NAME "Win8")
  endif(CMAKE_SYSTEM_VERSION MATCHES "6.2")
  
  # Cross-compiling
  if (CMAKE_CROSSCOMPILING)
     set(SPECIFIC_SYSTEM_VERSION_NAME "cross_${CMAKE_HOST_SYSTEM_NAME}")
  endif()

  # Compilers
  # taken from http://sourceforge.net/p/predef/wiki/Compilers/
  if(MSVC)
    set(SPECIFIC_COMPILER_NAME "MSVC-Unknown-${MSVC_VERSION}")
    if(MSVC_VERSION EQUAL 1200)
      set(SPECIFIC_COMPILER_NAME "MSVC-6.0")
    endif(MSVC_VERSION EQUAL 1200)
    if(MSVC_VERSION EQUAL 1300)
      set(SPECIFIC_COMPILER_NAME "MSVC-7.0")
    endif(MSVC_VERSION EQUAL 1300)
    if(MSVC_VERSION EQUAL 1310)
      set(SPECIFIC_COMPILER_NAME "MSVC-7.1-2003") #Visual Studio 2003
    endif(MSVC_VERSION EQUAL 1310)
    if(MSVC_VERSION EQUAL 1400)
      set(SPECIFIC_COMPILER_NAME "MSVC-8.0-2005") #Visual Studio 2005
    endif(MSVC_VERSION EQUAL 1400)
    if(MSVC_VERSION EQUAL 1500)
      set(SPECIFIC_COMPILER_NAME "MSVC-9.0-2008") #Visual Studio 2008
    endif(MSVC_VERSION EQUAL 1500)
    if(MSVC_VERSION EQUAL 1600)
      set(SPECIFIC_COMPILER_NAME "MSVC-10.0-2010") #Visual Studio 2010
    endif(MSVC_VERSION EQUAL 1600)
    if(MSVC_VERSION EQUAL 1700)
      set(SPECIFIC_COMPILER_NAME "MSVC-11.0-2012") #Visual Studio 2012
    endif(MSVC_VERSION EQUAL 1700)
    if(MSVC_VERSION EQUAL 1800)
      set(SPECIFIC_COMPILER_NAME "MSVC-12.0-2013") #Visual Studio 2013
    endif(MSVC_VERSION EQUAL 1800)
  endif(MSVC)
  if(MINGW)
    set(SPECIFIC_COMPILER_NAME "MinGW")
  endif(MINGW)
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    set(SPECIFIC_SYSTEM_VERSION_NAME "${SPECIFIC_SYSTEM_VERSION_NAME}-x86_64")
  endif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
endif(WIN32)

# In the Linux case try to guess the distro name/type
# using either lsb_release program or fallback
# to the content of the /etc/issue file
if(UNIX)
  if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(SPECIFIC_SYSTEM_VERSION_NAME "${CMAKE_SYSTEM_NAME}")
    set(SPECIFIC_SYSTEM_PREFERED_CPACK_GENERATOR "TGZ")

    if (EXISTS "/etc/os-release")
      #We prefer os-release as it is well formatted and
      #standardised. Also, /etc/issue contains no version info on
      #modern distros (e.g. CentOS7)
      read_release(ID FROM /etc/os-release INTO DISTRO_ID)
      read_release(VERSION_ID FROM /etc/os-release INTO DISTRO_RELEASE)
    elseif (EXISTS "/etc/issue")
      #If /etc/os-release does not exist, we fall back to parsing
      #/etc/issue, which is pretty popular.
      set(LINUX_NAME "")
      file(READ "/etc/issue" LINUX_ISSUE)
      #CentOS case
      if(LINUX_ISSUE MATCHES "CentOS")
        string(REGEX MATCH "release ([0-9\\.]+)" CENTOS "${LINUX_ISSUE}")
        set(DISTRO_ID "CentOS")
        set(DISTRO_RELEASE "${CMAKE_MATCH_1}")
        # FIXME can we find that in /etc/issue
        set(DISTRO_CODENAME "")
      endif(LINUX_ISSUE MATCHES "CentOS")
      #Scientific Linux
      if(LINUX_ISSUE MATCHES "Scientific")
        string(REGEX MATCH "release ([0-9\\.]+)" SCIENTIFIC "${LINUX_ISSUE}")
        set(DISTRO_ID "ScientificLinux")
        set(DISTRO_RELEASE "${CMAKE_MATCH_1}")
        # FIXME can we find that in /etc/issue
        set(DISTRO_CODENAME "")
      endif(LINUX_ISSUE MATCHES "Scientific")
      # Fedora case
      if(LINUX_ISSUE MATCHES "Fedora")
        string(REGEX MATCH "release ([0-9]+)" FEDORA "${LINUX_ISSUE}")
        set(DISTRO_ID "Fedora")
        set(DISTRO_RELEASE "${CMAKE_MATCH_1}")
        # FIXME can we find that in /etc/issue
        set(DISTRO_CODENAME "")
      endif(LINUX_ISSUE MATCHES "Fedora")
      # Ubuntu case
      if(LINUX_ISSUE MATCHES "Ubuntu")
        string(REGEX MATCH "buntu ([0-9]+\\.[0-9]+)" UBUNTU "${LINUX_ISSUE}")
        set(DISTRO_ID "Ubuntu")
        set(DISTRO_RELEASE "${CMAKE_MATCH_1}")
        # FIXME can we find that in /etc/issue
        set(DISTRO_CODENAME "")
      endif(LINUX_ISSUE MATCHES "Ubuntu")
      # Debian case
      if(LINUX_ISSUE MATCHES "Debian")
        string(REGEX MATCH "Debian .*ux ([0-9]+\\.[0-9]+)"
               DEBIAN "${LINUX_ISSUE}")
        set(DISTRO_ID "Debian")
        set(DISTRO_RELEASE "${CMAKE_MATCH_1}")
        set(DISTRO_CODENAME "")
      endif(LINUX_ISSUE MATCHES "Debian")
      # Open SuSE case
      if(LINUX_ISSUE MATCHES "SUSE")
        string(REGEX MATCH "SUSE ([0-9]+\\.[0-9]+)" SUSE "${LINUX_ISSUE}")
        set(DISTRO_ID "SUSE")
        set(DISTRO_RELEASE "${CMAKE_MATCH_1}")
        set(DISTRO_CODENAME "")
      endif(LINUX_ISSUE MATCHES "SUSE")
      # Mandriva case
      # TODO
    endif(EXISTS "/etc/os-release")
    string(TOLOWER ${DISTRO_ID} DISTRO_ID)
    # Now mangle some names
    set(LINUX_NAME "${DISTRO_ID}_${DISTRO_RELEASE}")
    if(DISTRO_ID MATCHES "fedora|mandriva|suse|opensuse|centos|scientific")
      set(SPECIFIC_SYSTEM_PREFERED_CPACK_GENERATOR "RPM")
    endif(DISTRO_ID MATCHES "fedora|mandriva|suse|opensuse|centos|scientific")
    if(DISTRO_ID MATCHES "debian|ubuntu")
      set(SPECIFIC_SYSTEM_PREFERED_CPACK_GENERATOR "DEB")
    endif(DISTRO_ID MATCHES "debian|ubuntu")
    if(LINUX_NAME)
      set(SPECIFIC_SYSTEM_VERSION_NAME "${CMAKE_SYSTEM_NAME}-${LINUX_NAME}")
    endif(LINUX_NAME)
  endif(CMAKE_SYSTEM_NAME MATCHES "Linux")
  set(SPECIFIC_SYSTEM_VERSION_NAME
     "${SPECIFIC_SYSTEM_VERSION_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
  set(SPECIFIC_COMPILER_NAME "")
endif(UNIX)
