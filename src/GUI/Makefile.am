# Copyright (C) 2012-2013  ALNEOS
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
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
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

# header files 
salomeinclude_HEADERS =

# Libraries targets
lib_LTLIBRARIES = libGMSHPluginGUI.la

dist_libGMSHPluginGUI_la_SOURCES = \
    GMSHPluginGUI.cxx \
    GMSHPluginGUI.h \
    GMSHPluginGUI_HypothesisCreator.h \
    GMSHPluginGUI_HypothesisCreator.cxx

MOC_FILES = \
    GMSHPluginGUI_HypothesisCreator_moc.cxx

nodist_libGMSHPluginGUI_la_SOURCES= \
    $(MOC_FILES)

# additionnal information to compil and link file

libGMSHPluginGUI_la_CPPFLAGS = \
    $(CAS_CPPFLAGS) \
    $(GMSH_INCLUDES) \
    $(QT_INCLUDES) \
    $(PYTHON_INCLUDES) \
    $(KERNEL_CXXFLAGS) \
    $(GUI_CXXFLAGS) \
    $(GEOM_CXXFLAGS) \
    $(MED_CXXFLAGS) \
    $(SMESH_CXXFLAGS) \
    $(BOOST_CPPFLAGS) \
    $(CORBA_CXXFLAGS) \
    $(CORBA_INCLUDES) \
    -I$(srcdir)/../GMSHPlugin \
    -I$(top_builddir)/idl

libGMSHPluginGUI_la_LDFLAGS  = \
    ../../idl/libSalomeIDLGMSHPLUGIN.la \
    ../GMSHPlugin/libGMSHEngine.la \
    ${QT_MT_LIBS} \
    ${GUI_LDFLAGS} -lSalomeApp -lqtx -lsuit -lSalomeObject -lLightApp \
    ${SMESH_LDFLAGS} -lSMESH -lGeomSelectionTools \
    $(CAS_KERNEL)

# resources files
nodist_salomeres_DATA= \
    GMSHPlugin_images.qm \
    GMSHPlugin_msg_en.qm \
    GMSHPlugin_msg_fr.qm

