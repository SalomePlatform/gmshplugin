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

Group_1 = model.addGroup(Part_1_doc, "Faces", [model.selection("FACE", "Box_1_1/Left")])

model.end()

#-------------------------------------
### SHAPERSTUDY component
#-------------------------------------
model.publishToShaperStudy()
import SHAPERSTUDY
Box_1_1, Group_1_1, = SHAPERSTUDY.shape(model.featureStringId(Box_1))

#-------------------------------------
### SMESH component
#-------------------------------------
from salome.kernel import SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()
#smesh.SetEnablePublish( False ) # Set to False to avoid publish in study if not needed or in some particular situations:
                                 # multiples meshes built in parallel, complex and numerous mesh edition (performance)

Gmsh_3D_Parameters = smesh.CreateHypothesis('GMSH_Parameters_3D', 'GMSHEngine')
Gmsh_3D_Parameters.Set3DAlgo( 0 )
Gmsh_3D_Parameters.SetMinSize( 0 )
Gmsh_3D_Parameters.SetIs2d( 0 )
Mesh_1 = smesh.Mesh(Box_1_1,'Mesh_1')
GMSH = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH)
Gmsh_Parameters = GMSH.Parameters()
Gmsh_Parameters.Set2DAlgo( 0 )
Gmsh_Parameters.SetMinSize( 0 )
Gmsh_Parameters.SetMaxSize( 1e+22 )
Gmsh_Parameters.SetIs2d( 0 )
Gmsh_Parameters.SetVerbosityLevel( 1 )
Viscous_Layers_1 = GMSH.ViscousLayers(1,1,1,[ 13 ],0,smeshBuilder.SURF_OFFSET_SMOOTH)
Group_1_2 = Mesh_1.GroupOnGeom(Group_1_1,'Group_1',SMESH.FACE)
Regular_1D = Mesh_1.Segment(geom=Group_1_1)
Number_of_Segments_1 = Regular_1D.NumberOfSegments(6)
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE,geom=Group_1_1)
isDone = Mesh_1.Compute()
Mesh_1.CheckCompute()
[ Group_1_2 ] = Mesh_1.GetGroups()
Sub_mesh_1 = Regular_1D.GetSubMesh()


## Set names of Mesh objects
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Group_1_2, 'Group_1')
smesh.SetName(Sub_mesh_1, 'Sub-mesh_1')
smesh.SetName(Viscous_Layers_1, 'Viscous Layers_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(Gmsh_3D_Parameters, 'Gmsh 3D Parameters')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(GMSH.GetAlgorithm(), 'GMSH')
smesh.SetName(Gmsh_Parameters, 'Gmsh Parameters')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')

assert Mesh_1.NbQuadrangles() == 60, "different number of quadrangles ({} instead of 60)".format(Mesh_1.NbTriangles())
assert Mesh_1.NbPyramids() == 36, "different number of pyramids ({} instead of 36)".format(Mesh_1.NbPyramids())

mesh_volume = Mesh_1.GetVolume()
expected_mesh_volume = 10**3
assert abs(mesh_volume-expected_mesh_volume)/expected_mesh_volume < 1e-12, "Wrong volume of mesh. Got {} instead of {}".format(mesh_volume, expected_mesh_volume)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()