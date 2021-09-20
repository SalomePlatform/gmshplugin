#!/usr/bin/env python

'''
This file creates a compound surface of a square and a circle. Then meshes them with 
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

### Create Sketch
Sketch_1 = model.addSketch(Part_1_doc, model.defaultPlane("XOY"))

### Create SketchLine
SketchLine_1 = Sketch_1.addLine(25, 0, 0, 0)

### Create SketchProjection
SketchProjection_1 = Sketch_1.addProjection(model.selection("VERTEX", "PartSet/Origin"), False)
SketchPoint_1 = SketchProjection_1.createdFeature()
Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchPoint_1.result())

### Create SketchLine
SketchLine_2 = Sketch_1.addLine(0, 0, 0, 25)

### Create SketchLine
SketchLine_3 = Sketch_1.addLine(0, 25, 25, 25)

### Create SketchLine
SketchLine_4 = Sketch_1.addLine(25, 25, 25, 0)
Sketch_1.setCoincident(SketchLine_4.endPoint(), SketchLine_1.startPoint())
Sketch_1.setCoincident(SketchLine_1.endPoint(), SketchLine_2.startPoint())
Sketch_1.setCoincident(SketchLine_2.endPoint(), SketchLine_3.startPoint())
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchLine_4.startPoint())
Sketch_1.setHorizontal(SketchLine_1.result())
Sketch_1.setVertical(SketchLine_2.result())
Sketch_1.setHorizontal(SketchLine_3.result())
Sketch_1.setVertical(SketchLine_4.result())
Sketch_1.setEqual(SketchLine_1.result(), SketchLine_4.result())
Sketch_1.setLength(SketchLine_1.result(), 25)

### Create SketchCircle
SketchCircle_1 = Sketch_1.addCircle(25, 25, 12.5)
Sketch_1.setCoincident(SketchLine_3.endPoint(), SketchCircle_1.center())
Sketch_1.setRadius(SketchCircle_1.results()[1], 12.5)
model.do()

### Create Face
Face_1_objects = [model.selection("FACE", "Sketch_1/Face-SketchLine_3f-SketchLine_4f-SketchCircle_1_2f-SketchCircle_1_2f"),
                  model.selection("FACE", "Sketch_1/Face-SketchLine_4r-SketchCircle_1_2r-SketchLine_3r-SketchLine_2r-SketchLine_1r"),
                  model.selection("FACE", "Sketch_1/Face-SketchCircle_1_2f-SketchLine_4r-SketchLine_3r")]
Face_1 = model.addFace(Part_1_doc, Face_1_objects)

### Create Group
Group_1_objects = [model.selection("FACE", "Face_1_2"),
                   model.selection("FACE", "Face_1_3"),
                   model.selection("FACE", "Face_1_1")]
Group_1 = model.addGroup(Part_1_doc, "Faces", Group_1_objects)

### Create Shell
Shell_1_objects = [model.selection("FACE", "Face_1_1"),
                   model.selection("FACE", "Face_1_2"),
                   model.selection("FACE", "Face_1_3")]
Shell_1 = model.addShell(Part_1_doc, Shell_1_objects)

### Create Group
Group_2_objects = [model.selection("FACE", "Shell_1_1/Modified_Face&Face_1_2/Face_1_2"),
                   model.selection("FACE", "Shell_1_1/Modified_Face&Face_1_3/Face_1_3"),
                   model.selection("FACE", "Shell_1_1/Modified_Face&Face_1_1/Face_1_1")]
Group_2 = model.addGroup(Part_1_doc, "Faces", Group_2_objects)

model.end()

#-------------------------------------
### SHAPERSTUDY component
#-------------------------------------

model.publishToShaperStudy()
import SHAPERSTUDY
Shell_1_1, Group_2_1, = SHAPERSTUDY.shape(model.featureStringId(Shell_1))

#-------------------------------------
### SMESH component
#-------------------------------------

import  SMESH
from salome.smesh import smeshBuilder
smesh = smeshBuilder.New()

Mesh_1 = smesh.Mesh(Shell_1_1)
GMSH_2D = Mesh_1.Triangle(algo=smeshBuilder.GMSH_2D)
Gmsh_Parameters = GMSH_2D.Parameters()
Gmsh_Parameters.Set2DAlgo( 0 )
Gmsh_Parameters.SetMaxSize( 10 )
Gmsh_Parameters.SetMinSize( 5 )
Gmsh_Parameters.SetIs2d( 1 )
Gmsh_Parameters.SetCompoundOnShape(Group_2_1)
isDone = Mesh_1.Compute()

errorMsg=''
okMsg=''

#-------------------------------------
# Test: Frontal Delaunay
#-------------------------------------
try:
  isDone = Mesh_1.Compute()
  if not isDone:
    errorMsg+= '\n ERROR: failed to mesh the compound square-circle surface using Delaunay algorithm from Gmsh\n'
  else:
    okMsg+= '\n PASSED: Successfully meshed the compound square-circle surface using Delaunay algorithm from Gmsh\n'
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
  smesh.SetName(GMSH_2D.GetAlgorithm(), 'GMSH_2D')
  smesh.SetName(Mesh_1.GetMesh(), 'Mesh_compound_square-circle')
  smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
  salome.sg.updateObjBrowser()
