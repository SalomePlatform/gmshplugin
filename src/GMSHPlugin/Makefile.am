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
salomeinclude_HEADERS = \
        GMSHPlugin_GMSH_2D.hxx \
        GMSHPlugin_GMSH_2D_i.hxx \
        GMSHPlugin_GMSH.hxx \
        GMSHPlugin_GMSH_i.hxx \
        GMSHPlugin_Hypothesis.hxx \
        GMSHPlugin_Hypothesis_i.hxx \
        GMSHPlugin_Hypothesis_2D.hxx \
        GMSHPlugin_Hypothesis_2D_i.hxx \
        GMSHPlugin_Mesher.hxx \
        GMSHPlugin_Defs.hxx

# Libraries targets
lib_LTLIBRARIES = libGMSHEngine.la

dist_libGMSHEngine_la_SOURCES = \
    GMSHPlugin_GMSH_2D.cxx \
    GMSHPlugin_GMSH_2D_i.cxx \
    GMSHPlugin_GMSH.cxx \
    GMSHPlugin_GMSH_i.cxx \
    GMSHPlugin_Hypothesis.cxx \
    GMSHPlugin_Hypothesis_i.cxx \
    GMSHPlugin_Hypothesis_2D.cxx \
    GMSHPlugin_Hypothesis_2D_i.cxx \
    GMSHPlugin_Mesher.cxx \
    GMSHPlugin_i.cxx

libGMSHEngine_la_CPPFLAGS = \
    $(KERNEL_CXXFLAGS) \
    $(GUI_CXXFLAGS) \
    $(MED_CXXFLAGS) \
    $(GEOM_CXXFLAGS) \
    $(CAS_CPPFLAGS) \
    $(VTK_INCLUDES) \
    -I$(GMSHHOME)/include \
    $(SMESH_CXXFLAGS) \
    $(CORBA_CXXFLAGS) \
    $(CORBA_INCLUDES) \
    $(BOOST_CPPFLAGS) \
    -I$(top_builddir)/idl
    

libGMSHEngine_la_LDFLAGS = -L$(GMSHHOME)/lib -lGmsh

libGMSHEngine_la_LDFLAGS  += \
    ../../idl/libSalomeIDLGMSHPLUGIN.la \
    $(CAS_LDPATH) -lTKernel -lTKBRep -lTKShHealing -lTKSTEP -lTKXSBase -lTKIGES -lTKMesh -lTKSTL -lTKG3d -lTKTopAlgo -lTKG2d -lTKBool -lTKGeomAlgo -lTKOffset -lTKGeomBase -lTKBO \
                  -lTKMath -lTKFillet -lTKMeshVS -lTKPrim -lTKSTEPBase -lTKSTEPAttr -lTKSTEP209 -lTKXDESTEP -lTKXDEIGES -lTKXCAF -lTKLCAF -lFWOSPlugin \
    $(GEOM_LDFLAGS) -lGEOMbasic \
    $(MED_LDFLAGS) -lSalomeIDLMED \
    $(SMESH_LDFLAGS) -lSMESHimpl -lSMESHEngine -lSMESHUtils -lStdMeshersEngine -lStdMeshers -lSMESHDS -lSMDS -lSMESHControls \
    $(KERNEL_LDFLAGS) -lSalomeGenericObj -lSalomeNS -lSALOMELocalTrace -lOpUtil

# Scripts to be installed.
mypkgpythondir = $(salomepythondir)/salome/GMSHPlugin
mypkgpython_PYTHON = \
	__init__.py \
	GMSHPluginBuilder.py
