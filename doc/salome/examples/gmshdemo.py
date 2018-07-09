# -*- coding: iso-8859-1 -*-

import salome

salome.salome_init()

from salome.geom import geomBuilder
geompy = geomBuilder.New()

from salome.smesh import smeshBuilder
smesh = smeshBuilder.New()

from salome.GMSHPlugin import GMSHPluginBuilder

# create a box
Box = geompy.MakeBoxDXDYDZ(10, 10, 10)
geompy.addToStudy( Box, 'Box' )

# 1. Create a 2D mesh on the box with GMSH_2D algorithm
Mesh_2D = smesh.Mesh(Box, "Box : 2D mesh by GMSH_2D")
# create a Gmsh 2D algorithm for solids
Algo_2D = Mesh_2D.Triangle(algo=smeshBuilder.GMSH_2D)
# define hypotheses
Param_2D = Algo_2D.Parameters()
# define algorithm
Param_2D.Set2DAlgo( 0 )
# define min element
Param_2D.SetMinSize( 0 )
# define max element
Param_2D.SetMaxSize( 2 )

# 2. Create a 3D mesh on the box with GMSH_3D algorithm
Mesh_3D = smesh.Mesh(Box, "Box : 3D mesh by GMSH_3D")
# create a Gmsh 3D algorithm for solids
Algo_3D = Mesh_3D.Tetrahedron(algo=smeshBuilder.GMSH)
# define hypotheses
Param_3D = Algo_3D.Parameters()
# define algorithms
Param_3D.Set2DAlgo( 0 )
Param_3D.SetIs2d( 0 )
Param_3D.Set3DAlgo( 0 )
# define min element size
Param_3D.SetMinSize( 0 )
# define max element size
Param_3D.SetMaxSize( 2 )

# compute the meshes
Mesh_2D.Compute()
Mesh_3D.Compute()

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
