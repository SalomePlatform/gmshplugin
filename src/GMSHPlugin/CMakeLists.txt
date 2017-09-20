# Copyright (C) 2012-2015  ALNEOS
# Copyright (C) 2016  EDF R&D
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
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

# --- options ---
# additional include directories
INCLUDE_DIRECTORIES(
  ${KERNEL_INCLUDE_DIRS}
  ${CAS_INCLUDE_DIRS}
  ${GEOM_INCLUDE_DIRS}
  ${SMESH_INCLUDE_DIRS}
  ${VTK_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${OMNIORB_INCLUDE_DIR}
  ${GMSH_INCLUDE_DIRS}
  ${PROJECT_BINARY_DIR}/idl
)

# additional preprocessor / compiler flags
ADD_DEFINITIONS(  
  ${OMNIORB_DEFINITIONS}
  ${CAS_DEFINITIONS}
  ${BOOST_DEFINITIONS}
)

# libraries to link to
SET(_link_LIBRARIES
  ${CAS_TKernel}
  ${CAS_TKMath}
  ${CAS_TKBRep}
  ${CAS_TKGeomBase}
  ${CAS_TKG2d}
  ${CAS_TKG3d}
  ${CAS_TKTopAlgo}
  ${CAS_TKGeomAlgo}
  ${CAS_TKCDF}
  ${GMSH_LIBRARIES}
  ${KERNEL_SALOMELocalTrace}
  ${KERNEL_SALOMEBasics}
  ${KERNEL_SalomeGenericObj}
  ${KERNEL_SalomeNS}
  ${KERNEL_OpUtil}
  ${SMESH_SMDS}
  ${SMESH_SMESHDS}
  ${SMESH_SMESHimpl}
  ${SMESH_SMESHEngine}
  ${SMESH_MeshDriverGMF}
  SalomeIDLGMSHPLUGIN
)

# --- headers ---

# header files
SET(GMSHEngine_HEADERS
  GMSHPlugin_GMSH_2D.hxx
  GMSHPlugin_GMSH_2D_i.hxx
  GMSHPlugin_GMSH.hxx
  GMSHPlugin_GMSH_i.hxx
  GMSHPlugin_Hypothesis.hxx
  GMSHPlugin_Hypothesis_i.hxx
  GMSHPlugin_Hypothesis_2D.hxx
  GMSHPlugin_Hypothesis_2D_i.hxx
  GMSHPlugin_Mesher.hxx
  GMSHPlugin_Defs.hxx
)

# --- sources ---

# sources / static
SET(GMSHEngine_SOURCES
  GMSHPlugin_GMSH_2D.cxx
  GMSHPlugin_GMSH_2D_i.cxx
  GMSHPlugin_GMSH.cxx
  GMSHPlugin_GMSH_i.cxx
  GMSHPlugin_Hypothesis.cxx
  GMSHPlugin_Hypothesis_i.cxx
  GMSHPlugin_Hypothesis_2D.cxx
  GMSHPlugin_Hypothesis_2D_i.cxx
  GMSHPlugin_Mesher.cxx
  GMSHPlugin_i.cxx
)

# --- scripts ---

# scripts / static
SET(_bin_SCRIPTS
  __init__.py
  GMSHPluginBuilder.py
)

# --- rules ---

ADD_LIBRARY(GMSHEngine ${GMSHEngine_SOURCES})
TARGET_LINK_LIBRARIES(GMSHEngine ${_link_LIBRARIES} )
INSTALL(TARGETS GMSHEngine EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${SALOME_INSTALL_LIBS})

INSTALL(FILES ${GMSHEngine_HEADERS} DESTINATION ${SALOME_INSTALL_HEADERS})

SALOME_INSTALL_SCRIPTS("${_bin_SCRIPTS}" ${SALOME_INSTALL_PYTHON}/salome/GMSHPlugin)