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

INCLUDE(tests.set)

IF(SALOME_BUILD_GUI)
  SET(TEST_SCRIPTS ${TEST_NAMES})
ELSE()
  SET(TEST_SCRIPTS ${TUI_NAMES})
ENDIF()

SET(COMPONENT_NAME GMSHPLUGIN)
SET(SALOME_TEST_DRIVER "$ENV{KERNEL_ROOT_DIR}/bin/salome/appliskel/salome_test_driver.py")
SET(TIMEOUT        300)

SET(TEST_INSTALL_DIRECTORY ${SALOME_GMSHPLUGIN_INSTALL_TESTS})

# make test
SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

FOREACH(tfile ${TEST_SCRIPTS})
  SET(TEST_NAME ${COMPONENT_NAME}_${tfile})
  ADD_TEST(NAME ${TEST_NAME}
           COMMAND ${PYTHON_EXECUTABLE} -B ${SALOME_TEST_DRIVER} ${TIMEOUT} ${CMAKE_CURRENT_SOURCE_DIR}/${tfile}.py)
  SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES ENVIRONMENT "${tests_env}")
  SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES LABELS "${COMPONENT_NAME}")
ENDFOREACH()

# salome test
FOREACH(tfile ${TEST_NAMES})
  INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${tfile}.py
          DESTINATION ${TEST_INSTALL_DIRECTORY})
ENDFOREACH()

INSTALL(FILES CTestTestfileInstall.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)

INSTALL(FILES tests.set DESTINATION ${TEST_INSTALL_DIRECTORY})
