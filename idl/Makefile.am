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

BASEIDL_FILES= GMSHPlugin_Algorithm.idl

BASEIDL_FILES_PY=$(BASEIDL_FILES:%.idl=%_idl.py)

# This variable defines the files to be installed
dist_salomeidl_DATA = $(BASEIDL_FILES)

# GUI idl common library
lib_LTLIBRARIES = libSalomeIDLGMSHPLUGIN.la

# Sources built from idl files
nodist_libSalomeIDLGMSHPLUGIN_la_SOURCES = GMSHPlugin_AlgorithmSK.cc

# header files must be exported: other modules have to use this library
nodist_salomeinclude_HEADERS = $(BASEIDL_FILES:%.idl=%.hh)

libSalomeIDLGMSHPLUGIN_la_CPPFLAGS = \
    $(KERNEL_CXXFLAGS) \
    $(GEOM_CXXFLAGS) \
    $(SMESH_CXXFLAGS) \
    @CORBA_CXXFLAGS@ \
    @CORBA_INCLUDES@ \
    -I$(top_builddir)/idl

libSalomeIDLGMSHPLUGIN_la_LDFLAGS = -no-undefined -version-info=0:0:0
libSalomeIDLGMSHPLUGIN_la_LIBADD  = \
    $(KERNEL_LDFLAGS) -lSalomeIDLKernel \
    $(SMESH_LDFLAGS) -lSalomeIDLSMESH \
    $(GEOM_LDFLAGS) -lSalomeIDLGEOM \
    @CORBA_LIBS@

# These variables defines the building process of CORBA files
OMNIORB_IDL         = @OMNIORB_IDL@
OMNIORB_IDLCXXFLAGS = @OMNIORB_IDLCXXFLAGS@
OMNIORB_IDLPYFLAGS  = @OMNIORB_IDLPYFLAGS@ \
    -I$(top_builddir)/idl/salome \
    -I$(KERNEL_ROOT_DIR)/idl/salome \
    -I$(GEOM_ROOT_DIR)/idl/salome \
    -I$(SMESH_ROOT_DIR)/idl/salome

IDLCXXFLAGS = \
    -bcxx \
    @IDLCXXFLAGS@ \
    -I$(top_builddir)/idl/salome \
    -I$(KERNEL_ROOT_DIR)/idl/salome \
    -I$(GEOM_ROOT_DIR)/idl/salome \
    -I$(SMESH_ROOT_DIR)/idl/salome
IDLPYFLAGS  = \
    @IDLPYFLAGS@ \
    -I$(KERNEL_ROOT_DIR)/idl/salome \
    -I$(GEOM_ROOT_DIR)/idl/salome \
    -I$(SMESH_ROOT_DIR)/idl/salome

# potential problem on parallel make on the following - multiple outputs
SUFFIXES = .idl .hh SK.cc
.idlSK.cc:
	$(OMNIORB_IDL) $(IDLCXXFLAGS) $(OMNIORB_IDLCXXFLAGS) $<
.idl.hh:
	$(OMNIORB_IDL) $(IDLCXXFLAGS) $(OMNIORB_IDLCXXFLAGS) $<

install-exec-local: $(BASEIDL_FILES:%=$(top_srcdir)/idl/%)
	$(INSTALL) -d  $(DESTDIR)$(salomepythondir)
	ls $^ | while read file; do \
	  $(OMNIORB_IDL) $(IDLPYFLAGS) -C$(DESTDIR)$(salomepythondir) $$file ; \
	done

# we want to remove only staff generated for IDL files and nothing more
uninstall-local:
	@for modulen in GMSHPlugin ; do \
	  test -d $(DESTDIR)$(salomepythondir)/$${modulen} && echo "Removing $(DESTDIR)$(salomepythondir)/$${modulen}" && rm -rf $(DESTDIR)$(salomepythondir)/$${modulen} ; \
	  test -d $(DESTDIR)$(salomepythondir)/$${modulen}__POA && echo "Removing $(DESTDIR)$(salomepythondir)/$${modulen}__POA" && rm -rf $(DESTDIR)$(salomepythondir)/$${modulen}__POA ; \
	done ; \
	for filen in $(BASEIDL_FILES_PY) ; do \
	  echo "Removing $(DESTDIR)$(salomepythondir)/$${filen}" && rm -f $(DESTDIR)$(salomepythondir)/$${filen}* ; \
	done

mostlyclean-local:
	-rm -f *.hh *.cc .depidl

# we use cpp to generate dependencies between idl files.
# option x c tells the preprocessor to consider idl as a c file.
# if an idl is modified, all idl dependencies are rebuilt

.depidl: $(BASEIDL_FILES)
	@echo "" > $@
	@for dep in $^ dummy; do \
	  if [ $$dep != "dummy" ]; then \
	    echo Building dependencies for $$dep; \
	    $(CPP) $(C_DEPEND_FLAG) -x c -I$(srcdir) -I$(KERNEL_ROOT_DIR)/idl/salome -I$(GEOM_ROOT_DIR)/idl/salome -I$(SMESH_ROOT_DIR)/idl/salome $$dep 2>/dev/null | \
	    sed 's/\.o/\SK.cc/' >>$@; \
	  fi; \
	done ;

-include .depidl
