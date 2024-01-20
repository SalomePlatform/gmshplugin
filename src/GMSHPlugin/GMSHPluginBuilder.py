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

##
# @package GMSHPluginBuilder
# Python API for the GMSH meshing plug-in module.

from salome.smesh.smesh_algorithm import Mesh_Algorithm

# import GMSHPlugin module if possible
noGMSHPlugin = 0
try:
    import GMSHPlugin
except ImportError:
    noGMSHPlugin = 1
    pass

# Types of algorithms
GMSH = "GMSH"
GMSH_3D       = "GMSH_3D"
GMSH_2D       = "GMSH_2D"
GMSH_3D_Remote = "GMSH_3D_Remote"

## Base of all GMSH algorithms.
#
class GMSH_Algorithm(Mesh_Algorithm):
    
    meshMethod = "Tetrahedron"
    algoType   = GMSH
    
    def __init__(self, mesh, geom=0):
        Mesh_Algorithm.__init__(self)
        if noGMSHPlugin: print("Warning: GMSHPlugin module unavailable")
        self.Create(mesh, geom, self.algoType, "libGMSHEngine.so")
        self.params = None

    ## Defines hypothesis having several parameters
    #
    def Parameters(self):
        if self.algoType == GMSH_2D:
            hypType = "GMSH_Parameters_2D"
        elif self.algoType == GMSH:
            hypType = "GMSH_Parameters"
        elif self.algoType == GMSH_3D:
            hypType = "GMSH_Parameters_3D"
        elif self.algoType == GMSH_3D_Remote:
            hypType = "GMSH_Parameters_3D"

        if self.params and self.params.GetName() != hypType:
            self.mesh.RemoveHypothesis( self.params, self.geom )
            self.params = None
        if not self.params:
            self.params = self.Hypothesis(hypType, [],"libGMSHEngine.so",UseExisting=0)
        return self.params

class GMSH_2D_Algorithm(GMSH_Algorithm):
    
    meshMethod = "Triangle"
    algoType   = GMSH_2D
    
    ## Private constructor.
    def __init__(self, mesh, geom=0):
        GMSH_Algorithm.__init__(self, mesh, geom)

class GMSH_3D_Algorithm(GMSH_Algorithm):
    
    meshMethod = "Tetrahedron"
    algoType   = GMSH_3D

    ## Private constructor.
    def __init__(self, mesh, geom=0):
        GMSH_Algorithm.__init__(self, mesh, geom)


class GMSH_3D_Remote_Algorithm(GMSH_Algorithm):
    meshMethod = "Tetrahedron"
    algoType   = GMSH_3D_Remote

    ## Private constructor.
    def __init__(self, mesh, geom=0):
        GMSH_Algorithm.__init__(self, mesh, geom)