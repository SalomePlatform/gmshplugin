// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016-2023  EDF
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
// See http://www.alneos.com/ or email : contact@alneos.fr
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#include "GMSHPlugin_GMSH_3D.hxx"
#include "GMSHPlugin_GMSH_3D_SA.hxx"
#include "GMSHPlugin_Hypothesis_2D.hxx"
#include "GMSHPlugin_Mesher.hxx"

#include <SMESH_Gen.hxx>
#include <SMESH_ControlsDef.hxx>
#include <SMESHDS_Mesh.hxx>
#include <utilities.h>

#include <list>

//=============================================================================
/*!
 *  
 */
//=============================================================================

GMSHPlugin_GMSH_3D_SA::GMSHPlugin_GMSH_3D_SA()
  : GMSHPlugin_GMSH_3D(0, new SMESH_Gen() )
{
  MESSAGE("GMSHPlugin_GMSH_3D_SA::GMSHPlugin_GMSH_3D_SA");
  _name = "GMSH_3D_SA";
}

//=============================================================================
/*!
 *  
 */
//=============================================================================

GMSHPlugin_GMSH_3D_SA::~GMSHPlugin_GMSH_3D_SA()
{
  MESSAGE("GMSHPlugin_GMSH_3D_SA::~GMSHPlugin_GMSH_3D_SA");
}

void GMSHPlugin_GMSH_3D_SA::importGMSHParameters( const std::string hypo_file )
{
  GMSHPlugin_Hypothesis * hypParameters = new GMSHPlugin_Hypothesis(0, GetGen());
  std::ifstream myfile(hypo_file);
  std::string line;

  std::getline(myfile, line);
  int hasParams = std::stoi(line);
  if ( hasParams )
  {
    std::getline(myfile, line);
    hypParameters->Set2DAlgo( (GMSHPlugin_Hypothesis::Algo2D)(std::stoi(line)) );
    std::getline(myfile, line);
    hypParameters->Set3DAlgo( (GMSHPlugin_Hypothesis::Algo3D)(std::stoi(line)) );
    std::getline(myfile, line);
    hypParameters->SetRecomb2DAlgo( (GMSHPlugin_Hypothesis::Recomb2DAlgo)(std::stoi(line)) );
    std::getline(myfile, line);
    hypParameters->SetRecombineAll( std::stoi(line) );  
    std::getline(myfile, line);
    hypParameters->SetSubdivAlgo( (GMSHPlugin_Hypothesis::SubdivAlgo)(std::stoi(line)) );  
    std::getline(myfile, line);
    hypParameters->SetRemeshAlgo( (GMSHPlugin_Hypothesis::RemeshAlgo)(std::stoi(line)) );
    std::getline(myfile, line);
    hypParameters->SetSmouthSteps( std::stod(line) );
    std::getline(myfile, line);
    hypParameters->SetSizeFactor( std::stod(line) );
    std::getline(myfile, line);
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
    hypParameters->SetMeshCurvatureSize( std::stod(line) );
#endif
    std::getline(myfile, line);
    hypParameters->SetMaxSize( std::stod(line) );
    std::getline(myfile, line);
    hypParameters->SetMinSize( std::stod(line) );
    std::getline(myfile, line);
    hypParameters->SetSecondOrder( std::stoi(line) );
    std::getline(myfile, line);
    hypParameters->SetUseIncomplElem( std::stoi(line) );
    std::getline(myfile, line); // Get Compound names comma separated
    if ( line != "0" ) // Mark no presence of groups
    {
      std::stringstream compoundLine (line);
      std::string compound;
      while ( std::getline (compoundLine, compound, ',') ) 
        hypParameters->SetCompoundOnEntry( compound );
    }
    _hypothesis = dynamic_cast< const GMSHPlugin_Hypothesis *> (hypParameters);
  } 
  else
  {
    // in case the parameters are not defined the _hypothesis should remain undefined
    // so default parameters are used to mesh!
    // _hypothesis = NULL; 
  }   
}

/**
 * @brief Write a binary file containing information on the elements/nodes
 *        created by the mesher
 *
 * @param nodeVec mapping between the mesh id and the netgen structure id
 * @param mesher gmhsplugin mesher
 * @param new_element_file Name of the output file
 * @return true if there are some error
 */
void GMSHPlugin_GMSH_3D_SA::fillNewElementFile( std::vector< const SMDS_MeshNode* > &nodeVec, GMSHPlugin_Mesher &mesher, std::string new_element_file )
{
  GModel* gModel = mesher.GetGModel();
  const int numberOfNodes = nodeVec.size() - 1;
  const int numOfVolumens = gModel->getNumMeshElements( 3 );
  int numberOfTotalNodes  = numberOfNodes;
  std::map< MVertex *,int> nodeMap;
  bool isOK = ( numOfVolumens > 0 );
  if ( isOK && !new_element_file.empty() )
  {
    MESSAGE("Writting new elements");
    
    std::ofstream df(new_element_file, ios::out|ios::binary);
    double points[3];
    int    volumens[4];

    df.write((char*) &numberOfNodes, sizeof(int));

    // To get 3D elements we need to iterate in regions
    std::map<int,MVertex*> vertexIdToCoordinate;
    // Index nodes of new volumetric elements in order from the vol region
    for ( GModel::riter it = gModel->firstRegion(); it != gModel->lastRegion(); ++it)
    {
      GRegion *gRegion = *it;
      for( size_t i = 0; i < gRegion->mesh_vertices.size(); i++)
      {
        MVertex *v = gRegion->mesh_vertices[i];
        const SMDS_MeshNode * preMeshedNode =  mesher.PremeshedNode( v );
        if ( !preMeshedNode )
        {
          numberOfTotalNodes++;
          vertexIdToCoordinate[ numberOfTotalNodes ] = v;
          nodeMap.insert({ v, numberOfTotalNodes });
        }           
      }
    }

    df.write((char*) &numberOfTotalNodes, sizeof(int));

    for (int nodeIndex = 1 ; nodeIndex <= numberOfNodes; ++nodeIndex )
    {
      //Id of the point
      int id = nodeVec.at(nodeIndex)->GetID();
      df.write((char*) &id, sizeof(int));
    }
        
    // Writing info on new points
    for (int nodeIndex = numberOfNodes+1; nodeIndex <= numberOfTotalNodes; ++nodeIndex )
    {
      if ( vertexIdToCoordinate[nodeIndex] )
      {
        df.write((char *) &vertexIdToCoordinate[nodeIndex]->x(), sizeof(double) );
        df.write((char *) &vertexIdToCoordinate[nodeIndex]->y(), sizeof(double) );
        df.write((char *) &vertexIdToCoordinate[nodeIndex]->z(), sizeof(double) );
      }      
    }

    df.write((char*) &numOfVolumens, sizeof(int));
    for ( GModel::riter it = gModel->firstRegion(); it != gModel->lastRegion(); ++it)
    {
      GRegion *gRegion = *it;
      std::vector<MVertex*> verts;

      for( size_t i = 0; i < gRegion->getNumMeshElements(); i++)
      {
        MElement *element = gRegion->getMeshElement(i);
        verts.clear();
        element->getVertices(verts);
        for( MVertex* v : verts )
        {
          const SMDS_MeshNode * node = mesher.PremeshedNode( v );
          auto it = find(nodeVec.begin(), nodeVec.end(), node ); 
          int nodeId = node ? (it - nodeVec.begin()) : nodeMap[ v ];
          df.write((char*) &nodeId, sizeof(int) );      
        }
      }               
    }
    df.close();
  }
}

/**
 * @brief Fill the list of elements in order and associate the read orientation (if any) to then. 
 *        Index the elements and orientation to the ordered map so we are sure the write orientation 
 *        match the element from the load mesh. IF no orientation file was written, then all the faces
 *        of the read mesh match the face orientation of his original associated face
 * @param aMesh the loaded mesh
 * @param element_orientation_file name of the file where to read the orientation
 * @param listElements the ordered map of elements with his orientation associated
 * @return the listElements filled.
 */
void GMSHPlugin_GMSH_3D_SA::fillListOfElementsOriented( SMESH_Mesh& aMesh, std::string element_orientation_file, 
                                                        std::map<const SMDS_MeshElement*, bool, TIDCompare>& listElements )
{
   
  // fill the elements found in the surface to the listOfElements so it can be used by the mesher!
  SMESHDS_Mesh* meshDS = aMesh.GetMeshDS();
  std::map<vtkIdType, bool> elemOrientation;
  {
    // Setting all element orientation to false if there no element orientation file
    if(element_orientation_file.empty())
    {
      SMDS_ElemIteratorPtr iteratorElem = meshDS->elementsIterator(SMDSAbs_Face);
      while ( iteratorElem->more() ) // loop on elements on a geom face
      {
        // check mesh face
        const SMDS_MeshElement* elem = iteratorElem->next();
        listElements[elem] = false;
      }
    } 
    else 
    {
      std::ifstream df(element_orientation_file, ios::binary|ios::in);
      int nbElement;
      bool orient;

      vtkIdType id;
      df.read((char*)&nbElement, sizeof(int));

      for(int ielem=0;ielem<nbElement;++ielem){
        df.read((char*) &id, sizeof(vtkIdType));
        df.read((char*) &orient, sizeof(bool));
        elemOrientation[id] = orient;
      }
      df.close();
      //
      SMDS_ElemIteratorPtr iteratorElem = meshDS->elementsIterator(SMDSAbs_Face);
      while ( iteratorElem->more() ) // loop on elements on a geom face
      {
        // check mesh face
        const SMDS_MeshElement* elem = iteratorElem->next();
        // only index registered elements
        bool isIn = elemOrientation.count(elem->GetID())==1;
        if(!isIn) continue;
        listElements[elem] = elemOrientation[elem->GetID()];
      }
    }
  }  
}


/**
 * @brief Compute the mesh by first filling premeshed faces to gmsh and then calling the mesher for the upper dimension
 * @param aShape the loaded shape 
 * @param aMesh the read Mesh (contain 2D elements covering the entire geometry)
 * @param new_element_file output file containing info the elements created by the mesher
 * @param element_orientation_file Binary file containing the orientation of surface elemnts
 * @param output_mesh whether or not write the created elements into the mesh
 * @return negation of mesh fail: true, false
 * */
bool GMSHPlugin_GMSH_3D_SA::Compute( TopoDS_Shape &aShape, SMESH_Mesh& aMesh, std::string new_element_file,
                                      std::string element_orientation_file, bool output_mesh )
{
  GMSHPlugin_Mesher mesher(&aMesh, aShape,/*2d=*/false, true);
  std::vector< const SMDS_MeshNode* > nodeVec; // to save premeshed elements
  mesher.SetParameters(dynamic_cast<const GMSHPlugin_Hypothesis*>(_hypothesis));

  std::map<const SMDS_MeshElement*, bool, TIDCompare> listElements;
  fillListOfElementsOriented( aMesh, element_orientation_file, listElements );
  bool Compute = mesher.Compute3D( nodeVec, listElements, output_mesh );
  
  fillNewElementFile( nodeVec, mesher, new_element_file );
  mesher.finalizeGModel();
  return Compute;
}

/**
 * @brief Running the mesher on the given files
 *
 * @param input_mesh_file Mesh file (containing 2D elements)
 * @param shape_file Shape file (BREP or STEP format)
 * @param hypo_file Ascii file containing the gmsh parameters
 * @param element_orientation_file Binary file containing the orientation of surface elements
 * @param new_element_file output file containing info the elements created by the mesher
 * @param output_mesh_file output mesh file (if empty it will not be created)
 * @return int
 */
int GMSHPlugin_GMSH_3D_SA::run(const std::string input_mesh_file,
                                const std::string shape_file,
                                const std::string hypo_file,
                                const std::string element_orientation_file,
                                const std::string new_element_file,
                                const std::string output_mesh_file)
{
  std::unique_ptr<SMESH_Mesh> myMesh(_gen->CreateMesh(false));

  SMESH_DriverMesh::importMesh(input_mesh_file, *myMesh);

  // Importing shape
  TopoDS_Shape myShape;
  SMESH_DriverShape::importShape(shape_file, myShape);

  // Define _hypothesis to then be able to call SetParameters( hypothesis )
  importGMSHParameters(hypo_file);

  MESSAGE("Meshing with gmsh3d");
  int ret = Compute(myShape, *myMesh, new_element_file, element_orientation_file, !output_mesh_file.empty());

  if(ret){
    std::cerr << "Meshing failed" << std::endl;
    return ret;
  }

  if(!output_mesh_file.empty()){
    std::string meshName = "MESH";
    SMESH_DriverMesh::exportMesh(output_mesh_file, *myMesh, meshName);
  }

  return ret;
}