# Copyright (C) 2012-2015  ALNEOS
# Copyright (C) 2016-2024  EDF
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

INCLUDE(UseQtExt)

# --- options ---
# additional include directories
INCLUDE_DIRECTORIES(
  ${QT_INCLUDES}
  ${OpenCASCADE_INCLUDE_DIR}
  ${KERNEL_INCLUDE_DIRS}
  ${GUI_INCLUDE_DIRS}
  ${GEOM_INCLUDE_DIRS}
  ${SMESH_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${OMNIORB_INCLUDE_DIR}
  ${GMSH_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${PROJECT_BINARY_DIR}/idl
  ${PROJECT_SOURCE_DIR}/src/GMSHPlugin
)

# additional preprocessor / compiler flags
ADD_DEFINITIONS(
  ${QT_DEFINITIONS}
  ${OMNIORB_DEFINITIONS}
  ${OpenCASCADE_DEFINITIONS}
  ${BOOST_DEFINITIONS}
)

# libraries to link to
SET(_link_LIBRARIES
  ${QT_LIBRARIES}
  ${KERNEL_SALOMELocalTrace}
  ${GEOM_GEOM}
  ${SMESH_SMESH}
  ${SMESH_PluginUtils}
  ${SMESH_StdMeshersGUI}
  ${SMESH_SMESHFiltersSelection}
  ${GUI_suit}
  ${GUI_qtx}
  ${GUI_SalomeApp}
  ${OpenCASCADE_FoundationClasses_LIBRARIES}
  SalomeIDLGMSHPLUGIN
)

# --- headers ---

# header files / to be processed by moc
SET(_moc_HEADERS
  GMSHPluginGUI_HypothesisCreator.h
)

# header files / no moc processed
SET(_other_HEADERS
  GMSHPluginGUI.h
)

SET(GMSHPluginGUI_HEADERS ${_moc_HEADERS} ${_other_HEADERS})

# --- sources ---

# sources / moc wrappings
QT_WRAP_MOC(_moc_SOURCES ${_moc_HEADERS})

# sources / static
SET(_other_SOURCES
  GMSHPluginGUI.cxx
  GMSHPluginGUI_HypothesisCreator.cxx
)

# --- resources ---

# resource files / to be processed by lrelease
SET(_ts_RESOURCES
  GMSHPlugin_images.ts
  GMSHPlugin_msg_en.ts
  GMSHPlugin_msg_fr.ts
)

# sources / to compile
SET(GMSHPluginGUI_SOURCES ${_other_SOURCES} ${_moc_SOURCES})

# --- rules ---

ADD_LIBRARY(GMSHPluginGUI ${GMSHPluginGUI_SOURCES})
TARGET_LINK_LIBRARIES(GMSHPluginGUI ${_link_LIBRARIES} )
INSTALL(TARGETS GMSHPluginGUI EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${SALOME_INSTALL_LIBS})

INSTALL(FILES ${GMSHPluginGUI_HEADERS} DESTINATION ${SALOME_INSTALL_HEADERS})
QT_INSTALL_TS_RESOURCES("${_ts_RESOURCES}" "${SALOME_GMSHPLUGIN_INSTALL_RES_DATA}")
