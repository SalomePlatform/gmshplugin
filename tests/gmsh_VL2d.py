#!/usr/bin/env python

import sys
from salome.kernel import salome

salome.salome_init()
from salome.kernel import salome_notebook
notebook = salome_notebook.NoteBook()

#-------------------------------------
### SHAPER component
#-------------------------------------
from salome.shaper import model

model.begin()
partSet = model.moduleDocument()

Part_1 = model.addPart(partSet)
Part_1_doc = Part_1.document()

Box_1 = model.addBox(Part_1_doc, 10, 10, 10)

Face_1 = model.addFace(Part_1_doc, [model.selection("FACE", "Box_1_1/Left"), model.selection("FACE", "Box_1_1/Top")])

Compound_1 = model.addCompound(Part_1_doc, [model.selection("FACE", "Face_1_1"), model.selection("FACE", "Face_1_2")])

model.end()

#-------------------------------------
### SHAPERSTUDY component
#-------------------------------------
model.publishToShaperStudy()
import SHAPERSTUDY
Compound_1_1, = SHAPERSTUDY.shape(model.featureStringId(Compound_1))

#-------------------------------------
### SMESH component
#-------------------------------------
from salome.kernel import SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Compound_1_1,'Mesh_1')
GMSH_2D = Mesh_1.Triangle(algo=smeshBuilder.GMSH_2D)
Gmsh_Parameters = GMSH_2D.Parameters()
Gmsh_Parameters.Set2DAlgo( 0 )
Gmsh_Parameters.SetMinSize( 0 )
Gmsh_Parameters.SetMaxSize( 1e+22 )
Gmsh_Parameters.SetIs2d( 1 )
Gmsh_Parameters.SetVerbosityLevel( 1 )
Viscous_Layers_2D_1 = GMSH_2D.ViscousLayers2D(1,1,1,[ 17, 19 ],0)
Viscous_Layers_2D_1.SetTotalThickness( 1 )
Viscous_Layers_2D_1.SetNumberLayers( 1 )
Viscous_Layers_2D_1.SetStretchFactor( 1 )
Viscous_Layers_2D_1.SetEdges( [ 17 ], 0 )
isDone = Mesh_1.Compute()
Mesh_1.CheckCompute()


## Set names of Mesh objects
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
smesh.SetName(GMSH_2D.GetAlgorithm(), 'GMSH_2D')
smesh.SetName(Viscous_Layers_2D_1, 'Viscous Layers 2D_1')

assert Mesh_1.NbQuadrangles() == 6, "different number of quadrangles ({} instead of 6)".format(Mesh_1.NbTriangles())
assert Mesh_1.NbEdges() == 50, "different number of edges ({} instead of 50)".format(Mesh_1.NbEdges())
assert Mesh_1.NbTriangles() > 190, "wrong number of triangles ({} less than 200)".format(Mesh_1.NbFaces())

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()