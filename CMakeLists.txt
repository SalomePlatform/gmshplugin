# Copyright (C) 2012-2015  ALNEOS
# Copyright (C) 2016-2021  EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.alneos.com/ or email : contact@alneos.fr
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8 FATAL_ERROR)

# Project name
# ============
# original
PROJECT(SalomeGMSHPLUGIN C CXX)
# upper case
STRING(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC)

# Common CMake macros
# ===================
SET(CONFIGURATION_ROOT_DIR $ENV{CONFIGURATION_ROOT_DIR} CACHE PATH "Path to the Salome CMake configuration files")
IF(EXISTS ${CONFIGURATION_ROOT_DIR})
  LIST(APPEND CMAKE_MODULE_PATH "${CONFIGURATION_ROOT_DIR}/cmake")
  INCLUDE(SalomeMacros NO_POLICY_SCOPE)
ELSE()
  MESSAGE(FATAL_ERROR "We absolutely need the Salome CMake configuration files, please define CONFIGURATION_ROOT_DIR !")
ENDIF()

# Versioning
# ===========
SALOME_SETUP_VERSION(9.9.0)
MESSAGE(STATUS "Building ${PROJECT_NAME_UC} ${${PROJECT_NAME_UC}_VERSION} from \"${${PROJECT_NAME_UC}_GIT_SHA1}\"")

# Find KERNEL
# ===========
SET(KERNEL_ROOT_DIR $ENV{KERNEL_ROOT_DIR} CACHE PATH "Path to the Salome KERNEL")
IF(EXISTS ${KERNEL_ROOT_DIR})
  FIND_PACKAGE(SalomeKERNEL REQUIRED)
  KERNEL_WITH_CORBA() # check whether KERNEL built with CORBA
ELSE(EXISTS ${KERNEL_ROOT_DIR})
  MESSAGE(FATAL_ERROR "We absolutely need a Salome KERNEL, please define KERNEL_ROOT_DIR")
ENDIF(EXISTS ${KERNEL_ROOT_DIR})

# Platform setup
# ==============
INCLUDE(SalomeSetupPlatform)   # From KERNEL
# Always build libraries as shared objects:
SET(BUILD_SHARED_LIBS TRUE)

# User options 
# ============
OPTION(SALOME_BUILD_TESTS "Build SALOME tests" ON)
OPTION(SALOME_BUILD_DOC "Generate SALOME GMSHPLUGIN documentation" ON)

IF(SALOME_BUILD_TESTS)
  ENABLE_TESTING()
ENDIF()

# Advanced options:
OPTION(SALOME_BUILD_GUI "Enable GUI" ON)

##
## From KERNEL:
##

# Various
# Python
FIND_PACKAGE(SalomePythonInterp REQUIRED)
FIND_PACKAGE(SalomePythonLibs REQUIRED)
# Boost
FIND_PACKAGE(SalomeBoost REQUIRED)
# CORBA
FIND_PACKAGE(SalomeOmniORB REQUIRED)
FIND_PACKAGE(SalomeOmniORBPy REQUIRED)

# Other KERNEL optionals:  
IF(SALOME_BUILD_DOC)
  FIND_PACKAGE(SalomeDoxygen)
  SALOME_LOG_OPTIONAL_PACKAGE(Doxygen SALOME_BUILD_DOC)
ENDIF()

# Find GUI (optional)
# ===========
IF(SALOME_BUILD_GUI)
  SET(GUI_ROOT_DIR $ENV{GUI_ROOT_DIR} CACHE PATH "Path to the Salome GUI")
  IF(EXISTS ${GUI_ROOT_DIR})
    FIND_PACKAGE(SalomeGUI)
    SALOME_GUI_WITH_CORBA() # check whether GUI built with CORBA
    ADD_DEFINITIONS(${GUI_DEFINITIONS})
    INCLUDE_DIRECTORIES(${GUI_INCLUDE_DIRS})
  ELSE(EXISTS ${GUI_ROOT_DIR})
    MESSAGE(FATAL_ERROR "We absolutely need a Salome GUI, please define GUI_ROOT_DIR")
  ENDIF(EXISTS ${GUI_ROOT_DIR})
  SALOME_LOG_OPTIONAL_PACKAGE(SalomeGUI SALOME_BUILD_GUI)
ENDIF(SALOME_BUILD_GUI)

##
## From GUI:
##

IF(SALOME_BUILD_GUI)
  # Qt
  FIND_PACKAGE(SalomeQt5 REQUIRED)
ENDIF(SALOME_BUILD_GUI)
# OCCT
FIND_PACKAGE(SalomeOpenCASCADE REQUIRED)
# VTK
FIND_PACKAGE(SalomeVTK REQUIRED)

# Find GEOM
# ===========

SET(GEOM_ROOT_DIR $ENV{GEOM_ROOT_DIR} CACHE PATH "Path to the Salome GEOM")
IF(EXISTS ${GEOM_ROOT_DIR})
  FIND_PACKAGE(SalomeGEOM REQUIRED)
ELSE(EXISTS ${GEOM_ROOT_DIR})
  MESSAGE(FATAL_ERROR "We absolutely need a Salome GEOM, please define GEOM_ROOT_DIR")
ENDIF(EXISTS ${GEOM_ROOT_DIR})

# Find SMESH
# ===========

SET(SMESH_ROOT_DIR $ENV{SMESH_ROOT_DIR} CACHE PATH "Path to the Salome SMESH")
IF(EXISTS ${SMESH_ROOT_DIR})
  FIND_PACKAGE(SalomeSMESH REQUIRED)
  ADD_DEFINITIONS(${SMESH_DEFINITIONS})
  INCLUDE_DIRECTORIES(${SMESH_INCLUDE_DIRS})
  IF(SALOME_BUILD_GUI)
    IF(NOT SALOME_SMESH_BUILD_GUI)
      MESSAGE(FATAL_ERROR "We absolutely need a Salome SMESH with GUI, please set SALOME_BUILD_GUI=ON in SMESH module")
    ENDIF(NOT SALOME_SMESH_BUILD_GUI)
  ENDIF(SALOME_BUILD_GUI)
  IF(SALOME_BUILD_DOC)
    IF(NOT SALOME_SMESH_BUILD_DOC)
      MESSAGE(FATAL_ERROR "We absolutely need a Salome SMESH documentation to build GMSHPLUGIN documentation. "
			  "So, you should either set SALOME_BUILD_DOC=OFF in GMSHPLUGIN or set SALOME_BUILD_DOC=ON in SMESH module.")
    ENDIF(NOT SALOME_SMESH_BUILD_DOC)
  ENDIF(SALOME_BUILD_DOC)
ELSE(EXISTS ${SMESH_ROOT_DIR})
  MESSAGE(FATAL_ERROR "We absolutely need a Salome SMESH, please define SMESH_ROOT_DIR")
ENDIF(EXISTS ${SMESH_ROOT_DIR})

##
## GMSHPLUGIN specifics
##
FIND_PACKAGE(SalomeGMSH REQUIRED)

# Detection summary:
SALOME_PACKAGE_REPORT_AND_CHECK()

# Directories
# (default values taken from KERNEL)
# ===========
SET(SALOME_INSTALL_BINS "${SALOME_INSTALL_BINS}" CACHE PATH "Install path: SALOME binaries")
SET(SALOME_INSTALL_LIBS "${SALOME_INSTALL_LIBS}" CACHE PATH "Install path: SALOME libs")
SET(SALOME_INSTALL_IDLS "${SALOME_INSTALL_IDLS}" CACHE PATH "Install path: SALOME IDL files")
SET(SALOME_INSTALL_HEADERS "${SALOME_INSTALL_HEADERS}" CACHE PATH "Install path: SALOME headers")
SET(SALOME_INSTALL_SCRIPT_SCRIPTS "${SALOME_INSTALL_SCRIPT_SCRIPTS}" CACHE PATH 
   "Install path: SALOME scripts")
SET(SALOME_INSTALL_SCRIPT_DATA "${SALOME_INSTALL_SCRIPT_DATA}" CACHE PATH 
   "Install path: SALOME script data")
SET(SALOME_INSTALL_SCRIPT_PYTHON "${SALOME_INSTALL_SCRIPT_PYTHON}" CACHE PATH 
   "Install path: SALOME Python scripts")
SET(SALOME_INSTALL_APPLISKEL_SCRIPTS "${SALOME_INSTALL_APPLISKEL_SCRIPTS}" CACHE PATH 
   "Install path: SALOME application skeleton - scripts")
SET(SALOME_INSTALL_APPLISKEL_PYTHON "${SALOME_INSTALL_APPLISKEL_PYTHON}" CACHE PATH 
   "Install path: SALOME application skeleton - Python")
SET(SALOME_INSTALL_PYTHON "${SALOME_INSTALL_PYTHON}" CACHE PATH "Install path: SALOME Python stuff")
SET(SALOME_INSTALL_PYTHON_SHARED "${SALOME_INSTALL_PYTHON_SHARED}" CACHE PATH 
   "Install path: SALOME Python shared modules")
SET(SALOME_INSTALL_CMAKE_LOCAL "${SALOME_INSTALL_CMAKE_LOCAL}" CACHE PATH 
    "Install path: local SALOME CMake files") 
SET(SALOME_INSTALL_AMCONFIG_LOCAL "${SALOME_INSTALL_AMCONFIG_LOCAL}" CACHE PATH
  "Install path: local SALOME config files (obsolete, to be removed)")
SET(SALOME_INSTALL_RES "${SALOME_INSTALL_RES}" CACHE PATH "Install path: SALOME resources")
SET(SALOME_INSTALL_DOC "${SALOME_INSTALL_DOC}" CACHE PATH "Install path: SALOME documentation")

# Specific to GMSHPLUGIN:
SET(SALOME_GMSHPLUGIN_INSTALL_RES_DATA "${SALOME_INSTALL_RES}/gmshplugin" CACHE PATH 
    "Install path: SALOME GMSHPLUGIN specific data")
SET(SALOME_GMSHPLUGIN_INSTALL_TESTS "${SALOME_INSTALL_BINS}/test")

MARK_AS_ADVANCED(SALOME_INSTALL_BINS SALOME_INSTALL_LIBS SALOME_INSTALL_IDLS SALOME_INSTALL_HEADERS)
MARK_AS_ADVANCED(SALOME_INSTALL_SCRIPT_SCRIPTS SALOME_INSTALL_SCRIPT_DATA SALOME_INSTALL_SCRIPT_PYTHON)
MARK_AS_ADVANCED(SALOME_INSTALL_APPLISKEL_SCRIPTS  SALOME_INSTALL_APPLISKEL_PYTHON SALOME_INSTALL_CMAKE_LOCAL SALOME_INSTALL_RES)
MARK_AS_ADVANCED(SALOME_INSTALL_PYTHON SALOME_INSTALL_PYTHON_SHARED)
MARK_AS_ADVANCED(SALOME_INSTALL_AMCONFIG_LOCAL SALOME_INSTALL_DOC)
MARK_AS_ADVANCED(SALOME_GMSHPLUGIN_INSTALL_RES_DATA SALOME_GMSHPLUGIN_INSTALL_TESTS)

# Accumulate environment variables for GMSHPLUGIN module
SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH NOCHECK ${CMAKE_INSTALL_PREFIX}/${SALOME_INSTALL_BINS}
                                         ${CMAKE_INSTALL_PREFIX}/${SALOME_INSTALL_PYTHON})
SALOME_ACCUMULATE_ENVIRONMENT(LD_LIBRARY_PATH NOCHECK ${CMAKE_INSTALL_PREFIX}/${SALOME_INSTALL_LIBS})
# Sources 
# ========

ADD_SUBDIRECTORY(idl)
ADD_SUBDIRECTORY(resources)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(bin)
IF(SALOME_BUILD_TESTS)
  ADD_SUBDIRECTORY(tests)
ENDIF(SALOME_BUILD_TESTS)
IF(SALOME_BUILD_DOC)
  ADD_SUBDIRECTORY(doc)
ENDIF(SALOME_BUILD_DOC)

# Header configuration
# ====================
SALOME_CONFIGURE_FILE(GMSHPLUGIN_version.h.in GMSHPLUGIN_version.h INSTALL ${SALOME_INSTALL_HEADERS})

# Configuration export
# (here only the level 1 prerequisites are exposed)
# ====================
INCLUDE(CMakePackageConfigHelpers)

# List of targets in this project we want to make visible to the rest of the world.
# They all have to be INSTALL'd with the option "EXPORT ${PROJECT_NAME}TargetGroup"
SET(_${PROJECT_NAME}_exposed_targets 
  GMSHEngine SalomeIDLGMSHPLUGIN
)

IF(SALOME_BUILD_GUI)
  LIST(APPEND _${PROJECT_NAME}_exposed_targets 
    GMSHPluginGUI
  )
ENDIF(SALOME_BUILD_GUI)

# Add all targets to the build-tree export set
EXPORT(TARGETS ${_${PROJECT_NAME}_exposed_targets}
  FILE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake)

# Create the configuration files:
#   - in the build tree:

# Ensure the variables are always defined for the configure:
SET(SMESH_ROOT_DIR "${SMESH_ROOT_DIR}")
SET(MESHGEMS_ROOT_DIR "${MESHGEMS_ROOT_DIR}")
 
SET(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")

# Build variables that will be expanded when configuring Salome<MODULE>Config.cmake:
SALOME_CONFIGURE_PREPARE(MESHGEMS)

CONFIGURE_PACKAGE_CONFIG_FILE(${PROJECT_NAME}Config.cmake.in 
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION "${SALOME_INSTALL_CMAKE_LOCAL}"
    PATH_VARS CONF_INCLUDE_DIRS SALOME_INSTALL_CMAKE_LOCAL CMAKE_INSTALL_PREFIX
    SMESH_ROOT_DIR MESHGEMS_ROOT_DIR)

WRITE_BASIC_PACKAGE_VERSION_FILE(${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${${PROJECT_NAME_UC}_VERSION}
    COMPATIBILITY AnyNewerVersion)
  
# Install the CMake configuration files:
INSTALL(FILES
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "${SALOME_INSTALL_CMAKE_LOCAL}")

# Install the export set for use with the install-tree
INSTALL(EXPORT ${PROJECT_NAME}TargetGroup DESTINATION "${SALOME_INSTALL_CMAKE_LOCAL}" 
  FILE ${PROJECT_NAME}Targets.cmake)
