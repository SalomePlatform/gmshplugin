#!/usr/bin/env python

'''
This file creates a box using GEOM then meshes the same using Gmsh's different algorithms

1. Delaunay               ----> Gmsh_Parameters.Set3DAlgo( 0 )
2. Frontal Delaunay       ----> Gmsh_Parameters.Set3DAlgo( 1 )
3. MMG3D                  ----> Gmsh_Parameters.Set3DAlgo( 2 )
4. R-Tree                 ----> Gmsh_Parameters.Set3DAlgo( 3 )
5. Parallel Delaunay(HXT) ----> Gmsh_Parameters.Set3DAlgo( 4 )

This file is solely for the propose of testing and we do overwrite the meshes.
'''

import salome
salome.salome_init()

#------------------------------------ 
# GEOM: Creating a box of size 10^3
#------------------------------------
import GEOM
from salome.geom import geomBuilder

geompy = geomBuilder.New()
Box_1 = geompy.MakeBoxDXDYDZ(10, 10, 10)

#------------------------------------ 
# SMESH: Using Gmsh algorithm with size (3,10) (min,max)
#------------------------------------
import  SMESH
from salome.smesh import smeshBuilder
smesh = smeshBuilder.New()

Mesh_1 = smesh.Mesh(Box_1)
GMSH = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH)
Gmsh_Parameters = GMSH.Parameters()
Gmsh_Parameters.Set2DAlgo( 0 )
Gmsh_Parameters.SetMinSize( 3 )
Gmsh_Parameters.SetMaxSize( 10 )
Gmsh_Parameters.SetIs2d( 0 )


errorMsg=''
okMsg=''

#-------------------------------------
# Test: Frontal Delaunay
#-------------------------------------
try:
  Gmsh_Parameters.Set3DAlgo( 0 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using Delaunay algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using Delaunay algorithm from Gmsh\n'
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'

#-------------------------------------
# Test: Frontal Hex
#-------------------------------------
try:
  Gmsh_Parameters.Set3DAlgo( 1 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using Frontal Delaunay algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using Frontal Delaunay algorithm from Gmsh\n'
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'

#-------------------------------------
# Test: MMG3D
#-------------------------------------
try:
  Gmsh_Parameters.Set3DAlgo( 2 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using MMG3D algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using MMG3D algorithm from Gmsh\n'
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'
      
#-------------------------------------
# Test: R-Tree Algorithm
#-------------------------------------
try:
  Gmsh_Parameters.Set3DAlgo( 3 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using R-Tree algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using R-Tree algorithm from Gmsh\n'
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'


#-------------------------------------
# Test: R-Tree Algorithm
#-------------------------------------
try:
  Gmsh_Parameters.Set3DAlgo( 4 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using Parallel Delaunay HXT algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using Parallel Delaunay HXT algorithm from Gmsh\n'
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'

#-------------------------------------
# Message that test are OK or not
#-------------------------------------

if okMsg!= '':
  print (okMsg)
      
if errorMsg!= '':
  raise RuntimeError (errorMsg + "\n Test is KO.")

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
