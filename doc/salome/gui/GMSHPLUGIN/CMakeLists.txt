# Copyright (C) 2012-2015  ALNEOS
# Copyright (C) 2016-2019  EDF R&D
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

SALOME_CONFIGURE_FILE(doxyfile.in doxyfile)
SALOME_CONFIGURE_FILE(doxyfile_py.in doxyfile_py)
SALOME_CONFIGURE_FILE(static/header.html.in ${CMAKE_CURRENT_BINARY_DIR}/static/header.html)
SALOME_CONFIGURE_FILE(static/header_py.html.in ${CMAKE_CURRENT_BINARY_DIR}/static/header_py.html)

SET(DOC_SMESH_MeshersList GMSHPlugin)
SET(smesh_file "$(SMESH_ROOT_DIR)/bin/salome/collect_mesh_methods.py")
SET(kernel_file "$ENV{KERNEL_ROOT_DIR}/bin/salome/prepare_generating_doc.py")

SALOME_ACCUMULATE_ENVIRONMENT(SMESH_MeshersList NOCHECK ${DOC_SMESH_MeshersList})

SET(_cmd_options ${smesh_file}  -d -o tmp/smeshBuilder.py GMSHPlugin)
SALOME_GENERATE_ENVIRONMENT_SCRIPT(_cmd env_script "${PYTHON_EXECUTABLE}" "${_cmd_options}")            

ADD_CUSTOM_TARGET(usr_docs ${CMAKE_COMMAND} -E make_directory tmp
  COMMAND ${PYTHON_EXECUTABLE} ${kernel_file} -o tmp/GMSHPluginBuilder.py ${CMAKE_SOURCE_DIR}/src/GMSHPlugin/GMSHPluginBuilder.py
  COMMAND ${PYTHON_EXECUTABLE} ${kernel_file} -o tmp/smesh_algorithm.py ${SMESH_ROOT_DIR}/${SALOME_INSTALL_PYTHON}/salome/smesh/smesh_algorithm.py
  COMMAND ${_cmd}
  COMMAND ${DOXYGEN_EXECUTABLE} doxyfile_py
  COMMAND ${DOXYGEN_EXECUTABLE} doxyfile
  COMMAND ${CMAKE_COMMAND} -E remove_directory tmp   
  VERBATIM
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

INSTALL(CODE "EXECUTE_PROCESS(COMMAND \"${CMAKE_COMMAND}\" --build ${PROJECT_BINARY_DIR} --target usr_docs)")
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/GMSHPLUGIN DESTINATION ${SALOME_INSTALL_DOC}/gui)
INSTALL(FILES images/head.png DESTINATION ${SALOME_INSTALL_DOC}/gui/GMSHPLUGIN)
INSTALL(FILES images/head.png DESTINATION ${SALOME_INSTALL_DOC}/gui/GMSHPLUGIN/gmshpluginpy_doc)

SET(make_clean_files GMSHPLUGIN gmshpluginpy_doc.tag)
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${make_clean_files}")
