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

model.end()

#-------------------------------------
### SHAPERSTUDY component
#-------------------------------------

model.publishToShaperStudy()
import SHAPERSTUDY
Box_1_1, = SHAPERSTUDY.shape(model.featureStringId(Box_1))

#-------------------------------------
### SMESH component
#-------------------------------------

from salome.kernel import SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Mesh_1 = smesh.Mesh(Box_1_1,'Mesh_1')
Regular_1D = Mesh_1.Segment()
nb_seg = 4
Number_of_Segments_1 = Regular_1D.NumberOfSegments(nb_seg)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
GMSH_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH_3D)
Gmsh_3D_Parameters = GMSH_3D.Parameters()
Gmsh_3D_Parameters.Set3DAlgo( 0 )
Gmsh_3D_Parameters.SetMinSize( 0 )
Gmsh_3D_Parameters.SetIs2d( 0 )
isDone = Mesh_1.Compute()
Mesh_1.CheckCompute()


## Set names of Mesh objects
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Gmsh_3D_Parameters, 'Gmsh 3D Parameters')
smesh.SetName(GMSH_3D.GetAlgorithm(), 'GMSH_3D')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')

expected_nb_faces = nb_seg*nb_seg*6
expected_nb_pyra = expected_nb_faces
assert Mesh_1.NbFaces() == expected_nb_faces, "different number of faces ({} instead of {})".format(Mesh_1.NbFaces(), expected_nb_faces)
assert Mesh_1.NbPyramids() == expected_nb_pyra, "different number of pyramids ({} instead of {})".format(Mesh_1.NbPyramids(), expected_nb_pyra)

mesh_volume = Mesh_1.GetVolume()
expected_mesh_volume = 10**3
assert abs(mesh_volume-expected_mesh_volume)/expected_mesh_volume < 1e-12, "Wrong volume of mesh. Got {} instead of {}".format(mesh_volume, expected_mesh_volume)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
