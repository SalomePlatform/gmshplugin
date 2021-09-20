#!/usr/bin/env python

'''
This file creates a square using GEOM then meshes the same using Gmsh's different algorithms

1. Automatic              ----> Gmsh_Parameters.Set2DAlgo( 0 )
2. mesh Adapt             ----> Gmsh_Parameters.Set2DAlgo( 1 )
3. Delaunay               ----> Gmsh_Parameters.Set2DAlgo( 2 )
4. Frontal                ----> Gmsh_Parameters.Set2DAlgo( 3 )
5. Delaunay For Quads     ----> Gmsh_Parameters.Set2DAlgo( 3 )
6. Packing Parallelograms ----> Gmsh_Parameters.Set2DAlgo( 3 )

This file is solely for the propose of testing and we do overwrite the meshes.
'''

import salome
salome.salome_init()

#------------------------------------ 
# GEOM: Creating a square of size 10^2
#------------------------------------
import GEOM
from salome.geom import geomBuilder

geompy = geomBuilder.New()
Face_1 = geompy.MakeFaceHW(10, 10, 1)

#------------------------------------ 
# SMESH: Using Gmsh algorithm with size (1,5) (min,max)
#------------------------------------
import  SMESH
from salome.smesh import smeshBuilder
smesh = smeshBuilder.New()

Mesh_1 = smesh.Mesh(Face_1)
GMSH_2D = Mesh_1.Triangle(algo=smeshBuilder.GMSH_2D)
Gmsh_Parameters = GMSH_2D.Parameters()
Gmsh_Parameters.SetMinSize( 1 )
Gmsh_Parameters.SetMaxSize( 5 )
Gmsh_Parameters.SetIs2d( 1 )

errorMsg=''
okMsg=''


#-------------------------------------
# Test: Automatic
#-------------------------------------
try:
  Gmsh_Parameters.Set2DAlgo( 0 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using Automatic algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using Automatic algorithm from Gmsh\n'
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'

#-------------------------------------
# Test:  Mesh Adapt
#-------------------------------------
try:
  Gmsh_Parameters.Set2DAlgo( 1 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using Mesh Adapt algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using Mesh Adapt algorithm from Gmsh\n'    
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'

#-------------------------------------
# Test: Delaunay
#-------------------------------------      
try:
  Gmsh_Parameters.Set2DAlgo( 2 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using Delaunay algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using Delaunay algorithm from Gmsh\n'    
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'
      
#-------------------------------------
# Test: Frontal Algorithm
#-------------------------------------           
try:
  Gmsh_Parameters.Set2DAlgo( 3 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using  Frontal algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using Frontal algorithm from Gmsh\n'    
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'

#-------------------------------------
# Test: Delaunay For Quads Algorithm
#-------------------------------------           
try:
  Gmsh_Parameters.Set2DAlgo( 4 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using  Delaunay For Quads algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using Delaunay For Quads algorithm from Gmsh\n'    
except:
      errorMsg+='\n ERROR: Exception raised in Mesh computation'
      

#-------------------------------------
# Test: Packing Parallelolagrams  Algorithm
#-------------------------------------           
try:
  Gmsh_Parameters.Set2DAlgo( 5 )
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the box using  Packing Parallelolagrams algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the box using Packing Parallelolagrams algorithm from Gmsh\n'    
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
