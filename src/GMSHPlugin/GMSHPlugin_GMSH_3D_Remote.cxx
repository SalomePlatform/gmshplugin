// Copyright (C) 2007-2024  CEA, EDF, OPEN CASCADE
//
// Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

//=============================================================================
// File      : GMSHPlugin_GMSH_3D_Remote.hxx
// Created   : 09 Septembre 2023
// Author    : Cesar Conopoima (OCC)
// Project   : SALOME
//=============================================================================
//
//

#include "GMSHPlugin_GMSH_3D_Remote.hxx"
#include "GMSHPlugin_Hypothesis.hxx"
#include "Utils_SALOME_Exception.hxx"

#include <SMESH_Gen.hxx>
#include <SMESH_Mesh.hxx>
#include <SMESH_ParallelMesh.hxx>
#include <SMESH_MesherHelper.hxx>
#include <SMESH_DriverShape.hxx>
#include <SMESH_DriverMesh.hxx>
#include <SMESHDS_Mesh.hxx>
#include <SMESH_MeshLocker.hxx>
#include <SMESH_ProxyMesh.hxx>

#include <TopoDS.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

#include <QString>
#include <QProcess>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

//=============================================================================
/*!
 * Constructor
 */
//=============================================================================

GMSHPlugin_GMSH_3D_Remote::GMSHPlugin_GMSH_3D_Remote(int hypId, SMESH_Gen * gen)
  : GMSHPlugin_GMSH_3D(hypId, gen)
{
  _name = "GMSH_3D_Remote";
}

//=============================================================================
/*!
 * Destructor
 */
//=============================================================================

GMSHPlugin_GMSH_3D_Remote::~GMSHPlugin_GMSH_3D_Remote()
{
}

/**
 * @brief Fill the structure netgen_param with the information from the hypothesis
 *
 * @param param_file name of the file to saven the gmsh parameter
 * @param hyp the hypothesis
 */
void GMSHPlugin_GMSH_3D_Remote::exportGmshParams( const std::string param_file, const SMESHDS_Hypothesis* hyp )
{
  std::ofstream myfile(param_file);
  if ( const GMSHPlugin_Hypothesis* gmshHypo = dynamic_cast<const GMSHPlugin_Hypothesis*>(hyp) )
  {
    myfile << 1 << std::endl; // Mark existence of correct hypothesis
    myfile << (int) gmshHypo->Get2DAlgo()       << std::endl;
    myfile << (int) gmshHypo->Get3DAlgo()       << std::endl;
    myfile << (int) gmshHypo->GetRecomb2DAlgo() << std::endl;
    myfile << (int) gmshHypo->GetRecombineAll() << std::endl;
    myfile << (int) gmshHypo->GetSubdivAlgo()   << std::endl;
    myfile << (int) gmshHypo->GetRemeshAlgo()   << std::endl;
    myfile << (double) gmshHypo->GetSmouthSteps()  << std::endl;
    myfile << (double) gmshHypo->GetSizeFactor()   << std::endl;
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
    myfile << (double) gmshHypo->GetMeshCurvatureSize()   << std::endl;
#elif
    myfile << -1.0   << std::endl; // dummy value writte for conformity
#endif
    myfile << (double) gmshHypo->GetMaxSize()     << std::endl;
    myfile << (double) gmshHypo->GetMinSize()     << std::endl;
    myfile << (int) gmshHypo->GetSecondOrder()    << std::endl;
    myfile << (int) gmshHypo->GetUseIncomplElem() << std::endl;
    if ( gmshHypo->GetCompoundOnEntries().size() > 0 )
    {
      TCompound defCompounds =  gmshHypo->GetCompoundOnEntries();
      for (TCompound::const_iterator it = defCompounds.begin();  it != defCompounds.end(); ++it )
      {
        std::string token = (it !=--defCompounds.end()) ? "," : "";
        myfile << *it << token;
      }
      myfile << std::endl;
    }
    else
    {
      myfile << 0 << std::endl; // mark the absence of compounds
    }
  }
  else
  {
    //Mark tha absence of parameters in the file
    myfile << 0 << std::endl;
  }
  myfile.close();
}

//

/**
 * @brief write in a binary file the orientation for each surface element of the mesh
 *
 * @param aMesh The mesh
 * @param aShape the shape associated to the mesh
 * @param output_file name of the binary file
 */
void GMSHPlugin_GMSH_3D_Remote::exportElementOrientation(SMESH_Mesh& aMesh,
                                                          const TopoDS_Shape& aShape,
                                                          const std::string output_file)
{
  SMESH_MesherHelper helper(aMesh);
  SMESH_ProxyMesh::Ptr proxyMesh( new SMESH_ProxyMesh( aMesh ));
  std::map<vtkIdType, bool> elemOrientation;

  for ( TopExp_Explorer exFa( aShape, TopAbs_FACE ); exFa.More(); exFa.Next())
  {
    const TopoDS_Shape& aShapeFace = exFa.Current();
    int faceID = aMesh.GetMeshDS()->ShapeToIndex( aShapeFace );
    bool isRev = false;
    if ( helper.NbAncestors(aShapeFace, aMesh, aShape.ShapeType()) > 1 )
      // IsReversedSubMesh() can work wrong on strongly curved faces,
      // so we use it as less as possible
      isRev = helper.IsReversedSubMesh( TopoDS::Face( aShapeFace ));

    const SMESHDS_SubMesh * aSubMeshDSFace = proxyMesh->GetSubMesh( aShapeFace );
    if ( !aSubMeshDSFace ) continue;

    SMDS_ElemIteratorPtr iteratorElem = aSubMeshDSFace->GetElements();

    while ( iteratorElem->more() ) // loop on elements on a geom face
    {
      // check mesh face
      const SMDS_MeshElement* elem = iteratorElem->next();
      if ( !elem )
        error( COMPERR_BAD_INPUT_MESH, "Null element encounters");
      if ( elem->NbCornerNodes() != 3 )
        error( COMPERR_BAD_INPUT_MESH, "Not triangle element encounters");
      elemOrientation[elem->GetID()] = isRev;
    } // loop on elements on a face
  } // loop on faces of a SOLID or SHELL

  {
    std::ofstream df(output_file, ios::out|ios::binary);
    int size = elemOrientation.size();
    df.write((char*)&size, sizeof(int));
    for(std::map<vtkIdType,bool>::iterator iter = elemOrientation.begin(); iter != elemOrientation.end(); ++iter)
    {
      df.write((char*)&(iter->first), sizeof(vtkIdType));
      df.write((char*)&(iter->second), sizeof(bool));
    }
    
    df.close();
  }
}

/**
 * @brief Compute mesh associate to shape
 *
 * @param aMesh The mesh
 * @param aShape The shape
 * @return true fi there are some error
 */
bool GMSHPlugin_GMSH_3D_Remote::Compute(SMESH_Mesh&         aMesh,
                                           const TopoDS_Shape& aShape)
{
  {
    SMESH_MeshLocker myLocker(&aMesh);
    SMESH_Hypothesis::Hypothesis_Status hypStatus;
    GMSHPlugin_GMSH_3D::CheckHypothesis(aMesh, aShape, hypStatus); //in this call the _hypothesis is defined!
  }

  SMESH_ParallelMesh& aParMesh = dynamic_cast<SMESH_ParallelMesh&>(aMesh);
  // Temporary folder for run
#ifdef WIN32
  // On windows mesh does not have GetTmpFolder
  fs::path tmp_folder = aParMesh.GetTmpFolder() / fs::path("Volume-%%%%-%%%%");
#else
  fs::path tmp_folder = aParMesh.GetTmpFolder() / fs::unique_path(fs::path("Volume-%%%%-%%%%"));
#endif
  fs::create_directories(tmp_folder);
  // Using MESH2D generated after all triangles where created.
  fs::path mesh_file= aParMesh.GetTmpFolder() / fs::path("Mesh2D.med");
  fs::path element_orientation_file=tmp_folder / fs::path("element_orientation.dat");
  fs::path new_element_file=tmp_folder / fs::path("new_elements.dat");
  //fs::path tmp_mesh_file=tmp_folder / fs::path("tmp_mesh.med");
  // Not used kept for debug
  // fs::path output_mesh_file=tmp_folder / fs::path("output_mesh.med");
  fs::path shape_file=tmp_folder / fs::path("shape.brep");
  fs::path param_file=tmp_folder / fs::path("gmsh_param.txt");
  fs::path log_file=tmp_folder / fs::path("run.log");
  fs::path cmd_file=tmp_folder / fs::path("cmd.txt");  
  
  std::string mesh_name = "MESH";

  {
    SMESH_MeshLocker myLocker(&aMesh);
    //Writing Shape
    SMESH_DriverShape::exportShape(shape_file.string(), aShape);

    //Writing hypothesis to file
    exportGmshParams(param_file.string(), _hypothesis);

    // Exporting element orientation
    exportElementOrientation(aMesh, aShape, element_orientation_file.string());
  }

  // Calling run_mesher
  // Path to mesher script
  fs::path mesher_launcher = fs::path(std::getenv("SMESH_ROOT_DIR"))/
       fs::path("bin")/
       fs::path("salome")/
       fs::path("mesher_launcher.py");

  std::string s_program="python3";
  std::list<std::string> params;
  params.push_back(mesher_launcher.string());
  params.push_back("GMSH3D");
  params.push_back(mesh_file.string());
  params.push_back(shape_file.string());
  params.push_back(param_file.string());
  params.push_back("--elem-orient-file=" + element_orientation_file.string());
  params.push_back("--new-element-file=" + new_element_file.string());
  // params.push_back("--output-mesh-file=" + output_mesh_file.string());

   // Parallelism method parameters
  int method = aParMesh.GetParallelismMethod();
  if(method == ParallelismMethod::MultiThread){
    params.push_back("--method=local");
  } else if (method == ParallelismMethod::MultiNode){
    params.push_back("--method=cluster");
    params.push_back("--resource="+aParMesh.GetResource());
    params.push_back("--wc-key="+aParMesh.GetWcKey());
    params.push_back("--nb-proc=1");
    params.push_back("--nb-proc-per-node="+std::to_string(aParMesh.GetNbProcPerNode()));
    params.push_back("--nb-node="+std::to_string(aParMesh.GetNbNode()));
    params.push_back("--walltime="+aParMesh.GetWalltime());
  } else {
    throw SALOME_Exception("Unknown parallelism method "+method);
  }
  
  std::string cmd = "";
  cmd += s_program;
  for(auto arg: params){
    cmd += " " + arg;
  }
  MESSAGE("Running command: ");
  MESSAGE(cmd);
  // Writing command in cmd.log
  {
    std::ofstream flog(cmd_file.string());
    flog << cmd << endl;
  }

   // Building arguments for QProcess
  QString program = QString::fromStdString(s_program);
  QStringList arguments;
  for(auto arg : params){
    arguments << arg.c_str();
  }

  QString out_file = log_file.string().c_str();
  QProcess myProcess;
  // myProcess.setProcessChannelMode(QProcess::MergedChannels);
  myProcess.setProcessChannelMode(QProcess::ForwardedChannels);
  myProcess.setStandardOutputFile(out_file);

  myProcess.start(program, arguments);
  // Waiting for process to finish (argument -1 make it wait until the end of
  // the process otherwise it just waits 30 seconds)
  bool finished = myProcess.waitForFinished(-1);
  int ret = myProcess.exitCode();
  if(ret != 0 || !finished){
    // Run crahed
    std::string msg = "Issue with mesh_launcher: \n";
    msg += "See log for more details: " + log_file.string() + "\n";
    msg += cmd + "\n";
    throw SALOME_Exception(msg);
  }                 
  {
    SMESH_MeshLocker myLocker(&aMesh);
    // Binary file written from SA version of the mesher
    std::ifstream df(new_element_file.string(), ios::binary);

    int numberOfNodes;
    int totalNumberOfNodes;
    int numberOfVolumes;
    double meshNodes[3];
    int    volNodes[4];
    int nodeID;

    SMESH_MesherHelper helper(aMesh);
    // This function is mandatory for setElementsOnShape to work
    helper.IsQuadraticSubMesh(aShape);
    helper.SetElementsOnShape( true );

    // Number of nodes in intial mesh
    df.read((char*) &numberOfNodes, sizeof(int));
    // Number of nodes added by netgen
    df.read((char*) &totalNumberOfNodes, sizeof(int));
    // Filling nodevec (correspondence netgen numbering mesh numbering)
    std::vector< const SMDS_MeshNode* > nodeVec ( totalNumberOfNodes + 1 );
    SMESHDS_Mesh * meshDS = helper.GetMeshDS();
    for (int nodeIndex = 1 ; nodeIndex <= numberOfNodes; ++nodeIndex )
    {
      //Id of the point
      df.read((char*) &nodeID, sizeof(int));
      nodeVec.at(nodeIndex) = meshDS->FindNode(nodeID);
    }

    // Add new points and update nodeVec
    for (int nodeIndex = numberOfNodes +1; nodeIndex <= totalNumberOfNodes; ++nodeIndex )
    {
      df.read((char *) &meshNodes, sizeof(double)*3);
      nodeVec.at(nodeIndex) = helper.AddNode(meshNodes[0], meshNodes[1], meshNodes[2]);
    }

    // Add tetrahedrons
    df.read((char*) &numberOfVolumes, sizeof(int));
    for ( int elemIndex = 1; elemIndex <= numberOfVolumes; ++elemIndex )
    {
      df.read((char*) &volNodes, sizeof(int)*4);
      auto n0 = meshDS->FindNode(nodeVec[volNodes[0]]->GetID());
      auto n1 = meshDS->FindNode(nodeVec[volNodes[1]]->GetID());
      auto n2 = meshDS->FindNode(nodeVec[volNodes[2]]->GetID());
      auto n3 = meshDS->FindNode(nodeVec[volNodes[3]]->GetID());
      if ( n0 && n1 && n2 && n3 )
        helper.AddVolume( n0, n2, n1, n3 );      
    }
  }

  return true;
}

/**
 * @brief Assign submeshes to compute
 *
 * @param aSubMesh submesh to add
 */
void GMSHPlugin_GMSH_3D_Remote::setSubMeshesToCompute(SMESH_subMesh * aSubMesh)
{
  SMESH_MeshLocker myLocker(aSubMesh->GetFather());
  SMESH_Algo::setSubMeshesToCompute(aSubMesh);
}
