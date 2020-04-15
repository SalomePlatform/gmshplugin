# Copyright (C) 2012-2015  ALNEOS
# Copyright (C) 2016-2020  EDF R&D
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

INCLUDE(UseOmniORB)  # Provided by KERNEL

INCLUDE_DIRECTORIES(
  ${OMNIORB_INCLUDE_DIR}
  ${KERNEL_INCLUDE_DIRS}
  ${GEOM_INCLUDE_DIRS}
  ${SMESH_INCLUDE_DIRS}
  ${PROJECT_BINARY_DIR}/idl
)

SET(SalomeIDLGMSHPLUGIN_IDLSOURCES
  GMSHPlugin_Algorithm.idl
)

SET(_idl_include_dirs
  ${KERNEL_ROOT_DIR}/idl/salome
  ${GEOM_ROOT_DIR}/idl/salome
  ${SMESH_ROOT_DIR}/idl/salome
)

SET(_idl_link_flags
  ${KERNEL_SalomeIDLKernel}
  ${GEOM_SalomeIDLGEOM}
  ${SMESH_SalomeIDLSMESH}
)

OMNIORB_ADD_MODULE(SalomeIDLGMSHPLUGIN "${SalomeIDLGMSHPLUGIN_IDLSOURCES}" "${_idl_include_dirs}" "${_idl_link_flags}")
INSTALL(TARGETS SalomeIDLGMSHPLUGIN EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${SALOME_INSTALL_LIBS})
