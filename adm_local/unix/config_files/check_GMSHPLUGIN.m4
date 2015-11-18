dnl Copyright (C) 2012-2013  ALNEOS
dnl
dnl This library is free software; you can redistribute it and/or
dnl modify it under the terms of the GNU Lesser General Public
dnl License as published by the Free Software Foundation; either
dnl version 2.1 of the License.
dnl
dnl This library is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public
dnl License along with this library; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
dnl
dnl See http://www.alneos.com/ or email : contact@alneos.fr
dnl

#------------------------------------------------------------
#  Check availability of Salome GMSH mesh plugin module
#   distribution
#------------------------------------------------------------

AC_DEFUN([CHECK_GMSHPLUGIN],[

AC_CHECKING(for GMSH mesh plugin)

NGplugin_ok=no

GMSHPLUGIN_LDLAGS=""
GMSHPLUGIN_CXXFLAGS=""

AC_ARG_WITH(gmshplugin,
	    [  --with-gmshplugin=DIR root directory path of GMSH mesh plugin installation ],
	    GMSHPLUGIN_DIR="$withval",GMSHPLUGIN_DIR="")

if test "x$GMSHPLUGIN_DIR" == "x" ; then

# no --with-gmshplugin-dir option used

   if test "x$GMSHPLUGIN_ROOT_DIR" != "x" ; then

    # GMSHPLUGIN_ROOT_DIR environment variable defined
      GMSHPLUGIN_DIR=$GMSHPLUGIN_ROOT_DIR

   fi
# 
fi

if test -f ${GMSHPLUGIN_DIR}/lib${LIB_LOCATION_SUFFIX}/salome/libGMSHEngine.so ; then
   NGplugin_ok=yes
   AC_MSG_RESULT(Using GMSH mesh plugin distribution in ${GMSHPLUGIN_DIR})

   if test "x$GMSHPLUGIN_ROOT_DIR" == "x" ; then
      GMSHPLUGIN_ROOT_DIR=${GMSHPLUGIN_DIR}
   fi
   AC_SUBST(GMSHPLUGIN_ROOT_DIR)

   GMSHPLUGIN_LDFLAGS=-L${GMSHPLUGIN_DIR}/lib${LIB_LOCATION_SUFFIX}/salome
   GMSHPLUGIN_CXXFLAGS=-I${GMSHPLUGIN_DIR}/include/salome

   AC_SUBST(GMSHPLUGIN_LDFLAGS)
   AC_SUBST(GMSHPLUGIN_CXXFLAGS)

else
   AC_MSG_WARN("Cannot find compiled GMSH mesh plugin distribution")
fi

AC_MSG_RESULT(for GMSH mesh plugin: $NGplugin_ok)
 
])dnl
