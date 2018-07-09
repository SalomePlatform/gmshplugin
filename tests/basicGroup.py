# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v8.5.0 with dump python functionality
###

import sys
import salome

salome.salome_init()

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Box_1 = geompy.MakeBoxDXDYDZ(200, 200, 200)
Cylinder_1 = geompy.MakeCylinderRH(100, 300)
Sphere_1 = geompy.MakeSphereR(100)
Fuse_1 = geompy.MakeFuseList([Box_1, Cylinder_1, Sphere_1], True, True)
Group_1 = geompy.CreateGroup(Fuse_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Group_1, [46, 37])
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Box_1, 'Box_1' )
geompy.addToStudy( Cylinder_1, 'Cylinder_1' )
geompy.addToStudy( Sphere_1, 'Sphere_1' )
geompy.addToStudy( Fuse_1, 'Fuse_1' )
geompy.addToStudyInFather( Fuse_1, Group_1, 'Group_1' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(Fuse_1)
GMSH = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH)
Gmsh_Parameters = GMSH.Parameters()
Gmsh_Parameters.Set2DAlgo( 0 )
Gmsh_Parameters.SetMinSize( 1 )
Gmsh_Parameters.SetMaxSize( 20 )
Gmsh_Parameters.SetIs2d( 0 )
Gmsh_Parameters.SetCompoundOnShape(Group_1)
isDone = Mesh_1.Compute()


## Set names of Mesh objects
smesh.SetName(GMSH.GetAlgorithm(), 'GMSH')
smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
