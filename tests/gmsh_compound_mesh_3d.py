#!/usr/bin/env python

'''
This file creates a compound solid of cube and a sphere. Then meshes them with 
Gmsh using compound mesh feature of Gmsh.
'''
import salome
salome.salome_init()

#-------------------------------------
### SHAPER component
#-------------------------------------

from salome.shaper import model
model.begin()
partSet = model.moduleDocument()

### Create Part
Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

### Create Box
Box_1 = model.addBox(Part_1_doc, 10, 10, 10)

### Create Sphere
Sphere_1 = model.addSphere(Part_1_doc, model.selection("VERTEX", "PartSet/Origin"), 5)

### Create Fuse
Fuse_1 = model.addFuse(Part_1_doc, [model.selection("COMPOUND", "all-in-Box_1"), model.selection("COMPOUND", "all-in-Sphere_1")], keepSubResults = True)

### Create Group
Group_1_objects = [model.selection("FACE", "Box_1_1/Front"),
                   model.selection("FACE", "Box_1_1/Top"),
                   model.selection("FACE", "Fuse_1_1/Modified_Face&Box_1_1/Left"),
                   model.selection("FACE", "Box_1_1/Right"),
                   model.selection("FACE", "Fuse_1_1/Modified_Face&Box_1_1/Back"),
                   model.selection("FACE", "Fuse_1_1/Modified_Face&Box_1_1/Bottom")]
Group_1 = model.addGroup(Part_1_doc, "Faces", Group_1_objects)

model.end()

#-------------------------------------
### SHAPERSTUDY component
#-------------------------------------

model.publishToShaperStudy()
import SHAPERSTUDY
Fuse_1_1, Group_1_1, = SHAPERSTUDY.shape(model.featureStringId(Fuse_1))


#-------------------------------------
### SMESH component
#-------------------------------------

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

Mesh_1 = smesh.Mesh(Fuse_1_1)
GMSH = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH)
Gmsh_Parameters = GMSH.Parameters()
Gmsh_Parameters.Set2DAlgo( 0 )
Gmsh_Parameters.SetMinSize( 0.5 )
Gmsh_Parameters.SetMaxSize( 1 )
Gmsh_Parameters.SetIs2d( 0 )
Gmsh_Parameters.SetCompoundOnShape(Group_1_1)
Group_1_2 = Mesh_1.GroupOnGeom(Group_1_1,'Group_1',SMESH.FACE)
#isDone = Mesh_1.Compute()
#[ Group_1_2 ] = Mesh_1.GetGroups()

errorMsg=''
okMsg=''


#-------------------------------------
# Test: Frontal Delaunay
#-------------------------------------
try:
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the compound cube-sphere surface using Delaunay algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the compound cube-sphere surface using Delaunay algorithm from Gmsh\n'
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
  smesh.SetName(GMSH.GetAlgorithm(), 'GMSH')
  smesh.SetName(Group_1_2, 'Group_1')
  smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
  smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
  salome.sg.updateObjBrowser()
