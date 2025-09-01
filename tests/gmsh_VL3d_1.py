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
Box_2 = model.addBox(Part_1_doc, 10, 10, 10)
Translation_1 = model.addTranslation(Part_1_doc, [model.selection("SOLID", "Box_2_1")], vector = [15, 0, 0], keepSubResults = True)
Compound_1 = model.addCompound(Part_1_doc, [model.selection("SOLID", "Box_1_1"), model.selection("SOLID", "Translation_1_1")])
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
NETGEN_1D_2D = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( 2.87228 )
NETGEN_2D_Parameters_1.SetMinSize( 0.957427 )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 2 )
NETGEN_2D_Parameters_1.SetChordalError( -1 )
NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_1.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_1.SetCheckChartBoundary( 0 )
GMSH_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.GMSH_3D)
Gmsh_3D_Parameters = GMSH_3D.Parameters()
Gmsh_3D_Parameters.Set3DAlgo( 0 )
Gmsh_3D_Parameters.SetMinSize( 0 )
Gmsh_3D_Parameters.SetIs2d( 0 )
Gmsh_3D_Parameters.SetVerbosityLevel( 1 )
status = Mesh_1.RemoveHypothesis(NETGEN_1D_2D)
Regular_1D = Mesh_1.Segment()
nb_seg = 6
Number_of_Segments_1 = Regular_1D.NumberOfSegments(nb_seg,None,[])
Quadrangle_2D = Mesh_1.Quadrangle(algo=smeshBuilder.QUADRANGLE)
status = Mesh_1.RemoveHypothesis(NETGEN_2D_Parameters_1)
Viscous_Layers_1 = GMSH_3D.ViscousLayers(1,1,1,[ 14, 24, 28, 32, 34, 38, 4, 48, 58, 62, 66, 68 ],0,smeshBuilder.SURF_OFFSET_SMOOTH)
isDone = Mesh_1.Compute()
Mesh_1.CheckCompute()


## Set names of Mesh objects
smesh.SetName(Viscous_Layers_1, 'Viscous Layers_1')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Number_of_Segments_1, 'Number of Segments_1')
smesh.SetName(Mesh_1.GetMesh(), 'Mesh_1')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 2D Parameters_1')
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(NETGEN_1D_2D.GetAlgorithm(), 'NETGEN 1D-2D')
smesh.SetName(Gmsh_3D_Parameters, 'Gmsh 3D Parameters')
smesh.SetName(GMSH_3D.GetAlgorithm(), 'GMSH_3D')


expected_nb_faces = nb_seg*nb_seg*6*2  # 2 boxes, each with nb_seg*nb_seg*6 faces
expected_nb_pyra = expected_nb_faces
expected_nb_hexas = expected_nb_faces
assert Mesh_1.NbFaces() == expected_nb_faces, "different number of faces ({} instead of {})".format(Mesh_1.NbFaces(), expected_nb_faces)
assert Mesh_1.NbPyramids() == expected_nb_pyra, "different number of pyramids ({} instead of {})".format(Mesh_1.NbPyramids(), expected_nb_pyra)
assert Mesh_1.NbHexas() == expected_nb_hexas, "different number of hexas ({} instead of {})".format(Mesh_1.NbPyramids(), expected_nb_pyra)

mesh_volume = Mesh_1.GetVolume()
expected_mesh_volume = 10**3 * 2  # 2 boxes, each with volume 10**3
assert abs(mesh_volume-expected_mesh_volume)/expected_mesh_volume < 1e-12, "Wrong volume of mesh. Got {} instead of {}".format(mesh_volume, expected_mesh_volume)

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()