# Install script for directory: /Users/alessiapaccagnella/desktop/gp20-alessiapacca/libigl/external/embree/common

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/alessiapaccagnella/desktop/gp20-alessiapacca/libigl/tutorial/embree/common/sys/cmake_install.cmake")
  include("/Users/alessiapaccagnella/desktop/gp20-alessiapacca/libigl/tutorial/embree/common/math/cmake_install.cmake")
  include("/Users/alessiapaccagnella/desktop/gp20-alessiapacca/libigl/tutorial/embree/common/simd/cmake_install.cmake")
  include("/Users/alessiapaccagnella/desktop/gp20-alessiapacca/libigl/tutorial/embree/common/lexers/cmake_install.cmake")
  include("/Users/alessiapaccagnella/desktop/gp20-alessiapacca/libigl/tutorial/embree/common/tasking/cmake_install.cmake")
  include("/Users/alessiapaccagnella/desktop/gp20-alessiapacca/libigl/tutorial/embree/common/algorithms/cmake_install.cmake")

endif()

