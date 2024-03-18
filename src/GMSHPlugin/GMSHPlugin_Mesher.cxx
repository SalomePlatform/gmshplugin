// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016-2024  EDF
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
#include "GMSHPlugin_Mesher.hxx"
#include "GMSHPlugin_Hypothesis_2D.hxx"

#include <SMDS_MeshElement.hxx>
#include <SMDS_MeshNode.hxx>
#include <SMESHDS_Mesh.hxx>
#include <SMESH_Comment.hxx>
#include <SMESH_ComputeError.hxx>
#include <SMESH_Gen_i.hxx>
#include <SMESH_Mesh.hxx>
#include <SMESH_MesherHelper.hxx>
#include <SMESH_subMesh.hxx>
#include <StdMeshers_FaceSide.hxx>
#include <utilities.h>
#include <StdMeshers_QuadToTriaAdaptor.hxx>

//CAS
#include <BRep_Tool.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <gp_Pnt.hxx>

#include <gmsh.h>
#include <vector>
#include <limits>

#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include <MLine.h>
#include <MTriangle.h>
#include <MQuadrangle.h>
#if GMSH_MAJOR_VERSION >=4
#include <GmshGlobal.h>
#include <gmsh/Context.h>
#endif

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
#include <omp.h>
#endif

using namespace std;

namespace
{
  struct ShapeBounds
  {
    SBoundingBox3d _bounds;
    TopoDS_Shape   _shape;
  };

  //================================================================================
  /*!
   * \brief Retrieve ShapeBounds from a compound GEdge
   */
  //================================================================================

  bool getBoundsOfShapes( GEdge*                       gEdge,
                          std::vector< ShapeBounds > & topoEdges )
  {
    topoEdges.clear();
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    for ( size_t i = 0; i < gEdge->compound.size(); ++i )
    {
      GEdge* gE = static_cast< GEdge* >( gEdge->compound[ i ]);
      topoEdges.push_back( ShapeBounds{ gE->bounds(), *((TopoDS_Edge*)gE->getNativePtr()) });
    }
#else
    for ( size_t i = 0; i < gEdge->_compound.size(); ++i )
    {
      GEdge* gE = static_cast< GEdge* >( gEdge->_compound[ i ]);
      topoEdges.push_back( ShapeBounds{ gE->bounds(), *((TopoDS_Edge*)gE->getNativePtr()) });
    }
#endif
    return topoEdges.size();
  }

  //================================================================================
  /*!
   * \brief Retrieve ShapeBounds from a compound GFace
   */
  //================================================================================

  bool getBoundsOfShapes( GFace*                       gFace,
                          std::vector< ShapeBounds > & topoFaces )
  {
    topoFaces.clear();
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    for ( size_t i = 0; i < gFace->compound.size(); ++i )
    {
      GFace* gF = static_cast< GFace* >( gFace->compound[ i ]);
      topoFaces.push_back( ShapeBounds{ gF->bounds(), *((TopoDS_Face*)gF->getNativePtr()) });
    }
#else
    for ( size_t i = 0; i < gFace->_compound.size(); ++i )
    {
      GFace* gF = static_cast< GFace* >( gFace->_compound[ i ]);
      topoFaces.push_back( ShapeBounds{ gF->bounds(), *((TopoDS_Face*)gF->getNativePtr()) });
    }
#endif
    return topoFaces.size();
  }
  //================================================================================
  /*!
   * \brief Find a shape whose bounding box includes a given point
   */
  //================================================================================

  TopoDS_Shape getShapeAtPoint( const SPoint3& point, const std::vector< ShapeBounds > & shapes )
  {
    TopoDS_Shape shape;
    float distmin = std::numeric_limits<float>::max();
    for ( size_t i = 0; i < shapes.size(); ++i )
    {
      float dist = GMSHPlugin_Mesher::DistBoundingBox( shapes[i]._bounds, point );
      if (dist < distmin)
      {
        shape = shapes[i]._shape;
        distmin = dist;
        if ( distmin == 0. )
          break;
      }
    }
    return shape;
  }

  double segmentSize( const UVPtStructVec& nodeParam, size_t i )
  {
    double l1 = SMESH_NodeXYZ( nodeParam[i].node ).Distance( nodeParam[i-1].node );
    double l2 = SMESH_NodeXYZ( nodeParam[i].node ).Distance( nodeParam[i+1].node );
    return 0.5 * ( l1 + l2 );
  }
}

//=============================================================================
/*!
 *
 */
//=============================================================================

GMSHPlugin_Mesher::GMSHPlugin_Mesher (SMESH_Mesh*         mesh,
                                      const TopoDS_Shape& aShape,
                                      bool                is2D,
                                      bool                is3D)
  : _mesh    (mesh),
    _shape   (aShape),
    _is2d    (is2D),
    _is3d    (is3D)
{
  // il faudra peut être mettre un truc par defaut si l'utilisateur ne rentre rien en para
  //defaultParameters();
}

//void GMSHPlugin_Mesher::defaultParameters(){}

void GMSHPlugin_Mesher::SetParameters(const GMSHPlugin_Hypothesis* hyp)
{
  if (hyp != NULL)
  {
    _algo2d          = hyp->Get2DAlgo();
    _algo3d          = hyp->Get3DAlgo();
    _recomb2DAlgo    = hyp->GetRecomb2DAlgo();
    _recombineAll    = hyp->GetRecombineAll();
    _subdivAlgo      = hyp->GetSubdivAlgo();
    _remeshAlgo      = hyp->GetRemeshAlgo();
    _remeshPara      = hyp->GetRemeshPara();
    _smouthSteps     = hyp->GetSmouthSteps();
    _sizeFactor      = hyp->GetSizeFactor();
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
    _meshCurvatureSize = hyp->GetMeshCurvatureSize();
#endif
    _minSize         = hyp->GetMinSize();
    _maxSize         = hyp->GetMaxSize();
    _secondOrder     = hyp->GetSecondOrder();
    _useIncomplElem  = hyp->GetUseIncomplElem();
    _verbLvl         = hyp->GetVerbosityLevel();
    _compounds       = hyp->GetCompoundOnEntries();
    // 6 in the enum corresponds to 99 in gmsh
    if(_verbLvl == 6)
      _verbLvl = 99;
  }
  else
  {
    _algo2d          = 0;
    _algo3d          = 0;
    _recomb2DAlgo    = 0;
    _recombineAll    = false;
    _subdivAlgo      = 0;
    _remeshAlgo      = 0;
    _remeshPara      = 0;
    _smouthSteps     = 1;
    _sizeFactor      = 1;
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
    _meshCurvatureSize = 0;
#endif
    _minSize         = 0;
    _maxSize         = 1e22;
    _secondOrder     = false;
    _useIncomplElem  = true;
    _compounds.clear();
  }
}


//================================================================================
/*!
 * \brief Set maximum number of threads to be used by Gmsh
 */
//================================================================================

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
void GMSHPlugin_Mesher::SetMaxThreadsGmsh()
{
  MESSAGE("GMSHPlugin_Mesher::SetMaxThreadsGmsh");
  // compound meshing (_compounds.size() > 0) and quad meshing (_algo2d >= 5) will
  // not be multi-threaded
  if (_compounds.size() > 0 || _algo2d >= 5){
    _maxThreads = 1;
  }
  else
    _maxThreads = omp_get_max_threads();
}
#endif


//================================================================================
/*!
 * \brief Check that the passed nodes are all IN the face
 * \param element the element
 * \param F the geom Face
 * \param uvValues a vector of size elem->NbCornerNodes() to save the uv coordinate points on the face
 * \return true if all the nodes are IN the face
 */
 //================================================================================
bool GMSHPlugin_Mesher::IsAllNodesInSameFace( const SMDS_MeshElement* element, const TopoDS_Face& F, 
                                                std::vector<gp_XY>& uvValues )
{
  Handle(ShapeAnalysis_Surface) sprojector = new ShapeAnalysis_Surface( BRep_Tool::Surface( F ));
  double tol = BRep_Tool::MaxTolerance( F, TopAbs_FACE );
  int nbN = element->NbCornerNodes();
  gp_Pnt surfPnt(0,0,0);
  for ( int i = 0; i < nbN; ++i )
  {
    SMESH_NodeXYZ nXYZ( element->GetNode( i ) );
    gp_XY uv = sprojector->ValueOfUV( nXYZ, tol ).XY();
    surfPnt = sprojector->Value( uv );
    double dist = surfPnt.Distance( nXYZ );
    if ( dist > tol )
      return false;
    else
      uvValues[ i ] = uv;
  }
  return true;
}

//================================================================================
/*!
 * \brief Associate mesh elements to geometrical faces.
 * \param the list of elements
 * \return a map between faces (incremental order) and mesh elements found to be placed on the face
 */
 //================================================================================
std::map<int,std::vector<std::tuple<smIdType,bool,std::vector<gp_XY>>>> GMSHPlugin_Mesher::AssociateElementsToFaces( std::map<const SMDS_MeshElement*, bool, TIDCompare>& listElements )
{
  // Map faces to elementId and uv of nodes.
  // Index by face id
  // Index vector with element smIdType, [ gp_XY_0, gp_XY_1, gp_XY_2 ]
  std::map<int,std::vector<std::tuple<smIdType,bool,std::vector<gp_XY>>>> elementToFaceMap;

  for(std::map<const SMDS_MeshElement*, bool>::iterator iter = listElements.begin(); iter != listElements.end(); ++iter)
  { 
    const SMDS_MeshElement* elem = iter->first;
    bool IsReverse = iter->second;
    int nbN = elem->NbCornerNodes();
    std::vector<gp_XY> uvValues(nbN);
    if ( nbN > 4 /*this restriction might be eliminated. Have to adapt FillGeomMapMeshUsing2DMeshIterator function too */)
      throw std::string("Polygon sub-meshes not supported");

    int faceId = 1;
    for( GModel::fiter it = _gModel->firstFace(); it != _gModel->lastFace(); ++it )
    {
      GFace *gFace = *it;
      TopoDS_Face topoFace = *((TopoDS_Face*)gFace->getNativePtr());
      if ( IsAllNodesInSameFace( elem, topoFace, uvValues ) )
      {
        elementToFaceMap[ faceId ].push_back( std::make_tuple( elem->GetID(), IsReverse, uvValues ) );     
        break;
      }
      faceId++;
    }
  }
  return elementToFaceMap;
}

//================================================================================
/*!
 * \brief Add the elements found associated to the face as gmsh elements
 */
 //================================================================================
void GMSHPlugin_Mesher::Set2DMeshes( std::vector< const SMDS_MeshNode* >& nodeVec, std::map<const SMDS_MeshElement*, bool, TIDCompare>& listElements )
{
  std::map< const SMDS_MeshNode* , const MVertex * > nodes2mvertMap;
  SMESHDS_Mesh* meshDS = _mesh->GetMeshDS();
  auto elementToFaceMap = AssociateElementsToFaces( listElements );
  int faceId = 1;
  std::vector<MVertex *> mVertices;

  for(GModel::fiter it = _gModel->firstFace(); it != _gModel->lastFace(); ++it)
  {
    GFace *gFace = *it;
    gFace->deleteMesh();

    auto element2uv = elementToFaceMap.find( faceId )->second;

    const int numberOfEntries = element2uv.size();
    
    for (int el = 0; el < numberOfEntries; el++)
    {
      const smIdType elementId      = std::get<0>( element2uv[ el ] ); // smesh element id
      bool isReverse                = std::get<1>( element2uv[ el ] );
      const SMDS_MeshElement* elem = meshDS->FindElement( elementId );

      int nbN = elem->NbCornerNodes();
      mVertices.resize( nbN );

      for ( int i = 0; i < nbN; ++i )
      {
        const SMDS_MeshNode* n = elem->GetNode( i );
        MVertex *           mv = nullptr;
        auto n2v = nodes2mvertMap.find( n );
        if ( n2v != nodes2mvertMap.end() )
        {
          mv = const_cast< MVertex*>( n2v->second );
        }
        else
        {
          if ( n->GetPosition()->GetDim() < 2 )
            throw std::string("Wrong mapping of edge nodes to GMSH nodes");
          SMESH_NodeXYZ xyz = n;
          gp_XY uv = std::get<2>(element2uv[ el ])[ i ];
          mv = new MFaceVertex( xyz.X(), xyz.Y(), xyz.Z(), gFace, uv.X(), uv.Y() );
          gFace->mesh_vertices.push_back( mv );
          nodes2mvertMap.insert({ n, mv });
          _nodeMap.insert      ({ mv, n });
          _premeshednodeMap.insert({ mv, n });
        }
        mVertices[ i ] = mv;
      }
       // create GMSH mesh faces
      switch ( nbN ) {
      case 3:
        if ( isReverse )
          gFace->triangles.push_back (new MTriangle(mVertices[0], mVertices[2], mVertices[1]));
        else
          gFace->triangles.push_back (new MTriangle(mVertices[0], mVertices[1], mVertices[2]));
        break;
      case 4:
        if ( isReverse )
          gFace->quadrangles.push_back (new MQuadrangle(mVertices[0], mVertices[3],
                                                        mVertices[2], mVertices[1]));
        else
          gFace->quadrangles.push_back (new MQuadrangle(mVertices[0], mVertices[1],
                                                        mVertices[2], mVertices[3]));
        break;
      default:;
      }
    }        
    faceId++; // face counter
  } // iterator in the face

  // Fill the node 
  nodeVec.resize( nodes2mvertMap.size() + 1, 0 );
  int count = 1;
  for (auto k : nodes2mvertMap )
  {
    nodeVec[ count ] = k.first; // Index the node id to the smesh node itself        
    count++;
  }
}


//================================================================================
/*!
 * \brief Initialize GMSH model with mesh elements as geometry objects. 
 *          Nodes are vertexes and element connections are geom lines
 */
 //================================================================================
void GMSHPlugin_Mesher::FillGeomMapMeshUsing2DMeshIterator( std::map<const SMDS_MeshElement*, bool, TIDCompare>& listElements )
{
  gmsh::initialize();
  gmsh::model::add("mesh");
    // typedef for maps
  typedef map< const SMDS_MeshNode*, int, TIDCompare > TNodeToIDMap;
  typedef TNodeToIDMap::value_type                     TN2ID;
  typedef map<std::pair<int, int>, int> TLineToIDMap;
  TNodeToIDMap aNodeToID;
  TLineToIDMap aLineToID;
  
  int aNbOfNodes = 0;
  int aNbOfLines = 0;

  const int invalid_ID = -1;
  std::vector<int> aTrinagle( 3, 0 );
  
  for(std::map<const SMDS_MeshElement*, bool>::iterator iter = listElements.begin(); iter != listElements.end(); ++iter)
  { 
    const SMDS_MeshElement* elem = iter->first;
    bool IsReverse = iter->second;
    if ( elem->NbCornerNodes() != 3 )
      return;

    for (int iN = 0; iN < 3; ++iN)
    {
      const SMDS_MeshNode* aNode = elem->GetNode(iN);
      
      int& ngID = aNodeToID.insert(TN2ID(aNode, invalid_ID)).first->second;
      if (ngID == invalid_ID)
      {
        ngID = ++aNbOfNodes;
        gmsh::model::occ::addPoint(aNode->X(), aNode->Y(), aNode->Z(), 1.e-2, ngID);
      }

      aTrinagle[ IsReverse ? 2 - iN : iN ] = ngID;
    }
    // add triangle
    if ((aTrinagle[0] == aTrinagle[1] ||
        aTrinagle[0] == aTrinagle[2] ||
        aTrinagle[2] == aTrinagle[1]))
      continue;
    
    std::vector<int> LinesID(3, 0);
    for (int anIndex = 0; anIndex < 3; ++anIndex)
    {
      int aNextIndex = (anIndex + 1) % 3;
      if (aLineToID.find({ aTrinagle[anIndex], aTrinagle[aNextIndex] }) == aLineToID.end()
        && aLineToID.find({ aTrinagle[aNextIndex], aTrinagle[anIndex] }) == aLineToID.end())
      {
        LinesID[anIndex] = aLineToID.insert({ { aTrinagle[aNextIndex], aTrinagle[anIndex] }, ++aNbOfLines }).first->second;
        gmsh::model::occ::addLine(aTrinagle[anIndex], aTrinagle[aNextIndex], LinesID[anIndex]);
      }
      else
      {
        LinesID[anIndex] = aLineToID.find({ aTrinagle[anIndex], aTrinagle[aNextIndex] })->second;
        if (LinesID[anIndex] == 0)
          LinesID[anIndex] = aLineToID.find({ aTrinagle[aNextIndex], aTrinagle[anIndex] })->second;

      }
    }
    // if (!aProxyMesh->IsTemporary(ls.first))
    //   swap(aTrinagle[1], aTrinagle[2]);
    gmsh::model::occ::addCurveLoop(LinesID);
  }
}

//================================================================================
/*!
 * \brief Initialize GMSH model
 */
 //================================================================================
void GMSHPlugin_Mesher::FillGMSHMesh()
{
  gmsh::initialize();
  gmsh::model::add("mesh");

  SMESHDS_Mesh* meshDS = _mesh->GetMeshDS();

  int aNbOfNodes = 0;
  int aNbOfLines = 0;
  std::vector<int> aTrinagle(3, 0);

  const int invalid_ID = -1;

  // typedef for maps
  typedef map< const SMDS_MeshNode*, int, TIDCompare > TNodeToIDMap;
  typedef TNodeToIDMap::value_type                     TN2ID;
  typedef map<std::pair<int, int>, int> TLineToIDMap;
  TNodeToIDMap aNodeToID;
  TLineToIDMap aLineToID;

  TopAbs_ShapeEnum aMainType = _mesh->GetShapeToMesh().ShapeType();
  bool aCheckReverse = (aMainType == TopAbs_COMPOUND || aMainType == TopAbs_COMPSOLID);

  SMESH_MesherHelper aHelper(*_mesh);
  SMESH_ProxyMesh::Ptr aProxyMesh(new SMESH_ProxyMesh(*_mesh));
  if (_mesh->NbQuadrangles() > 0)
  {
    StdMeshers_QuadToTriaAdaptor* Adaptor = new StdMeshers_QuadToTriaAdaptor;
    Adaptor->Compute(*_mesh, _shape, aProxyMesh.get());
    aProxyMesh.reset(Adaptor);
  }

  std::map<const SMDS_MeshElement*, bool, TIDCompare> listElements;
  for (TopExp_Explorer exFa(_shape, TopAbs_FACE); exFa.More(); exFa.Next())
  {
    const TopoDS_Shape& aShapeFace = exFa.Current();
    int faceID = meshDS->ShapeToIndex(aShapeFace);
    bool isRev = false;
    if (aCheckReverse && aHelper.NbAncestors(aShapeFace, *_mesh, _shape.ShapeType()) > 1)
      // IsReversedSubMesh() can work wrong on strongly curved faces,
      // so we use it as less as possible
      isRev = aHelper.IsReversedSubMesh(TopoDS::Face(aShapeFace));

    const SMESHDS_SubMesh* aSubMeshDSFace = aProxyMesh->GetSubMesh(aShapeFace);
    if (!aSubMeshDSFace) 
      continue;
    SMDS_ElemIteratorPtr iteratorElem = aSubMeshDSFace->GetElements();
    if (aHelper.IsQuadraticSubMesh(_shape) &&
      dynamic_cast<const SMESH_ProxyMesh::SubMesh*>(aSubMeshDSFace))
    {
      // add medium nodes of proxy triangles to helper
      while (iteratorElem->more())
        aHelper.AddTLinks(static_cast<const SMDS_MeshFace*>(iteratorElem->next()));

      iteratorElem = aSubMeshDSFace->GetElements();
    }
    while (iteratorElem->more()) // loop on elements on a geom face
    {
      // check mesh face
      const SMDS_MeshElement* elem = iteratorElem->next();
      if (!elem)
        return;
      if (elem->NbCornerNodes() != 3)
        return;
      listElements[elem] = isRev;
    }
  }

  for (auto const& ls : listElements)
  {
    // Add nodes of triangles and triangles them-selves to netgen mesh
    // add three nodes of
    bool hasDegen = false;
    for (int iN = 0; iN < 3; ++iN)
    {
      const SMDS_MeshNode* aNode = ls.first->GetNode(iN);
      const int shapeID = aNode->getshapeId();
      if (aNode->GetPosition()->GetTypeOfPosition() == SMDS_TOP_EDGE &&
        aHelper.IsDegenShape(shapeID))
      {
        // ignore all nodes on degeneraged edge and use node on its vertex instead
        TopoDS_Shape vertex = TopoDS_Iterator(meshDS->IndexToShape(shapeID)).Value();
        aNode = SMESH_Algo::VertexNode(TopoDS::Vertex(vertex), meshDS);
        hasDegen = true;
      }
      int& ngID = aNodeToID.insert(TN2ID(aNode, invalid_ID)).first->second;
      if (ngID == invalid_ID)
      {
        ngID = ++aNbOfNodes;
        gmsh::model::occ::addPoint(aNode->X(), aNode->Y(), aNode->Z(), 1.e-2, ngID);
      }
      aTrinagle[ls.second ? 2 - iN : iN] = ngID;
    }
    // add triangle
    if (hasDegen && (aTrinagle[0] == aTrinagle[1] ||
      aTrinagle[0] == aTrinagle[2] ||
      aTrinagle[2] == aTrinagle[1]))      
        continue;    


    std::vector<int> LinesID(3, 0);
    for (int anIndex = 0; anIndex < 3; ++anIndex)
    {
      int aNextIndex = (anIndex + 1) % 3;
      if (aLineToID.find({ aTrinagle[anIndex], aTrinagle[aNextIndex] }) == aLineToID.end()
        && aLineToID.find({ aTrinagle[aNextIndex], aTrinagle[anIndex] }) == aLineToID.end())
      {
        LinesID[anIndex] = aLineToID.insert({ { aTrinagle[aNextIndex], aTrinagle[anIndex] }, ++aNbOfLines }).first->second;
        gmsh::model::occ::addLine(aTrinagle[anIndex], aTrinagle[aNextIndex], LinesID[anIndex]);
      }
      else
      {
        LinesID[anIndex] = aLineToID.find({ aTrinagle[anIndex], aTrinagle[aNextIndex] })->second;
        if (LinesID[anIndex] == 0)
          LinesID[anIndex] = aLineToID.find({ aTrinagle[aNextIndex], aTrinagle[anIndex] })->second;

      }
    }
    if (!aProxyMesh->IsTemporary(ls.first))
      swap(aTrinagle[1], aTrinagle[2]);

    gmsh::model::occ::addCurveLoop(LinesID);
  }

  // Generate 1D and 2D mesh
  _gModel->mesh( /*dim=*/ 1);
  Set1DSubMeshes(_gModel);
  _gModel->mesh( /*dim=*/ 2);
}

//================================================================================
/*!
 * \brief Set Gmsh Options
 */
 //================================================================================
void GMSHPlugin_Mesher::SetGmshOptions()
{
  MESSAGE("GMSHPlugin_Mesher::SetGmshOptions");
  /*
  printf("We chose _algo2d         %d \n", _algo2d        );
  printf("We chose _algo3d         %d \n", _algo3d        );
  printf("We chose _recomb2DAlgo   %d \n", _recomb2DAlgo  );
  printf("We chose _recombineAll   %d \n", (_recombineAll)?1:0);
  printf("We chose _subdivAlgo     %d \n", _subdivAlgo    );
  printf("We chose _remeshAlgo     %d \n", _remeshAlgo    );
  printf("We chose _remeshPara     %d \n", _remeshPara    );
  printf("We chose _smouthSteps    %e \n", _smouthSteps   );
  printf("We chose _sizeFactor     %e \n", _sizeFactor    );
  printf("We chose _minSize        %e \n", _minSize       );
  printf("We chose _maxSize        %e \n", _maxSize       );
  printf("We chose _secondOrder    %d \n", (_secondOrder)?1:0);
  printf("We chose _useIncomplElem %d \n", (_useIncomplElem)?1:0);
  printf("We are in dimension      %d \n", (_is2d)?2:3);
  //*/

  std::map <int,double> mapAlgo2d;
  mapAlgo2d[0]=2; // Automatic
  mapAlgo2d[1]=1; // MeshAdapt
  mapAlgo2d[2]=5; // Delaunay
  mapAlgo2d[3]=6; // Frontal-Delaunay
  mapAlgo2d[4]=8; // DelQuad (Frontal-Delaunay for Quads)
  mapAlgo2d[5]=9; // Packing of parallelograms
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
  mapAlgo2d[6]=11;// Quasistructured quads with cross-fields
#endif

  std::map <int,double> mapAlgo3d;
  mapAlgo3d[0]=1; // Delaunay
  mapAlgo3d[1]=4; // Frontal
  mapAlgo3d[2]=7; // MMG3D
  mapAlgo3d[3]=9; // R-tree
  mapAlgo3d[4]=10;// HXT

  int ok;
  ok = GmshSetOption("Mesh", "Algorithm"                , mapAlgo2d[_algo2d])    ;
  ASSERT(ok);
  if ( !_is2d)
  {
    ok = GmshSetOption("Mesh", "Algorithm3D"            , mapAlgo3d[_algo3d])    ;
    ASSERT(ok);
  }
  ok = GmshSetOption("Mesh", "RecombinationAlgorithm"   , (double)_recomb2DAlgo) ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "RecombineAll"             , (_recombineAll)?1.:0.) ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "SubdivisionAlgorithm"     , (double)_subdivAlgo)   ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "RemeshAlgorithm"          , (double)_remeshAlgo)   ;
  //ASSERT(ok);
  ok = GmshSetOption("Mesh", "RemeshParametrization"    , (double)_remeshPara)   ;
  //ASSERT(ok);
  ok = GmshSetOption("Mesh", "Smoothing"                , (double)_smouthSteps)  ;
  //ASSERT(ok);
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
  ok = GmshSetOption("Mesh", "MeshSizeFromCurvature"       , _meshCurvatureSize) ;
  ASSERT(ok);
#endif
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
  ok = GmshSetOption("Mesh", "MeshSizeFactor"              , _sizeFactor)     ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "MeshSizeMin"                 , _minSize)        ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "MeshSizeMax"                 , _maxSize)        ;
  ASSERT(ok);
#else
  ok = GmshSetOption("Mesh", "CharacteristicLengthFactor"  , _sizeFactor)     ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "CharacteristicLengthMin"     , _minSize)        ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "CharacteristicLengthMax"     , _maxSize)        ;
  ASSERT(ok);
#endif
  ok = GmshSetOption("Mesh", "ElementOrder"             , (_secondOrder)?2.:1.)  ;
  ASSERT(ok);
  if (_secondOrder)
  {
    ok = GmshSetOption("Mesh", "SecondOrderIncomplete"  ,(_useIncomplElem)?1.:0.);
    ASSERT(ok);
  }

  ok = GmshSetOption("General", "Verbosity"            , (double) _verbLvl )  ; // Verbosity level
  ASSERT(ok);

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
/*ok = GmshSetOption("Mesh", "MaxNumThreads1D"          , 0. )  ; // Coarse-grain algo threads
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "MaxNumThreads2D"          , 0. )  ; // Coarse-grain algo threads
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "MaxNumThreads3D"          , 0. )  ; // Fine-grain algo threads HXT
  ASSERT(ok);
**/
  ok = GmshSetOption("General", "NumThreads"            , _maxThreads )  ; // system default i.e. OMP_NUM_THREADS
  ASSERT(ok);
#ifdef WIN32
  if ( GMSHPlugin_Hypothesis::Algo3D::hxt == _algo3d ){
    MESSAGE("GMSHPlugin_Mesher::SetGmshOptions: HXT algorithm is being used. Setting number of threads to 1.");
    ok = GmshSetOption("Mesh", "MaxNumThreads3D"       , 1. );
    ASSERT(ok);
  } // hxt
#endif
#endif
}

//================================================================================
/*!
 * \brief Create and add Compounds into GModel _gModel.
 */
//================================================================================

void GMSHPlugin_Mesher::CreateGmshCompounds()
{
  MESSAGE("GMSHPlugin_Mesher::CreateGmshCompounds");

  SMESH_Gen_i* smeshGen_i = SMESH_Gen_i::GetSMESHGen();

  OCC_Internals* occgeo = _gModel->getOCCInternals();
  bool toSynchronize = false;

  for(std::set<std::string>::const_iterator its = _compounds.begin();its != _compounds.end(); ++its )
  {
    GEOM::GEOM_Object_var aGeomObj;
    TopoDS_Shape geomShape = TopoDS_Shape();
    SALOMEDS::SObject_var aSObj = SMESH_Gen_i::GetSMESHGen()->getStudyServant()->FindObjectID( (*its).c_str() );
    SALOMEDS::GenericAttribute_var anAttr;
    if (!aSObj->_is_nil() && aSObj->FindAttribute(anAttr, "AttributeIOR"))
    {
      SALOMEDS::AttributeIOR_var anIOR = SALOMEDS::AttributeIOR::_narrow(anAttr);
      CORBA::String_var aVal = anIOR->Value();
      CORBA::Object_var obj = SMESH_Gen_i::GetSMESHGen()->getStudyServant()->ConvertIORToObject(aVal);
      aGeomObj = GEOM::GEOM_Object::_narrow(obj);
    }
    geomShape = smeshGen_i->GeomObjectToShape( aGeomObj.in() );
    if ( geomShape.IsNull() )
      continue;

    TopAbs_ShapeEnum geomType = geomShape.ShapeType();
    if ( geomType == TopAbs_COMPOUND)// voir s'il ne faut pas mettre une erreur dans le cas contraire
    {
      MESSAGE("shapeType == TopAbs_COMPOUND");
      TopoDS_Iterator it(geomShape);
      if ( !it.More() )
        continue;
      TopAbs_ShapeEnum shapeType = it.Value().ShapeType();
      std::vector< std::pair< int, int > > dimTags;
      for ( ; it.More(); it.Next())
      {
        const TopoDS_Shape& topoShape = it.Value();
        ASSERT(topoShape.ShapeType() == shapeType);
        if ( _mesh->GetMeshDS()->ShapeToIndex( topoShape ) > 0 )
          occgeo->importShapes( &topoShape, false, dimTags );
        else
        {
          TopoDS_Shape face = TopExp_Explorer( _shape, shapeType ).Current();
          SMESH_subMesh* sm = _mesh->GetSubMesh( face );
          sm->GetComputeError() =
            SMESH_ComputeError::New
            ( COMPERR_WARNING, "Compound shape does not belong to the main geometry. Ignored");
        }
      }
      std::vector<int> tags;
      int dim = ( shapeType == TopAbs_EDGE ) ? 1 : 2;
      for ( size_t i = 0; i < dimTags.size(); ++i )
      {
        if ( dimTags[i].first == dim )
          tags.push_back( dimTags[i].second );
      }
      if ( !tags.empty() )
      {
        _gModel->getGEOInternals()->setCompoundMesh( dim, tags );
        toSynchronize = true;
      }
      if ( toSynchronize )
        _gModel->getGEOInternals()->synchronize(_gModel);
    }
  }
}

//================================================================================
/*!
 * \brief For a compound mesh set the mesh components to be transmitted to SMESH
 */
//================================================================================

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
void GMSHPlugin_Mesher::SetCompoundMeshVisibility()
{

  // Loop over all faces, if the face belongs to a compound entry then
  // for all (boundary) edges whithin the face visibility is set to 0,
  // if the face doesn't belong to a compound entry then visibility is
  // set to 1 for all its (boundary) edges. Later, in FillSMesh() func
  // getVisibility() (returns either 1 or 0) is used to decide weather
  // the mesh of edge should be transmitted  to SMESH or not.

  for ( GModel::fiter itF = _gModel->firstFace(); itF != _gModel->lastFace(); ++itF )
  {
    std::vector< GEdge *> faceEdges = (*itF)->edges();

    for ( auto itE = faceEdges.begin(); itE != faceEdges.end(); ++itE )
    {
      if ( ((*itF)->compound.size()) )
        (*itE)->setVisibility(0);
      else
        (*itE)->setVisibility(1);
    }
  }


  // Loop over all edges, if the  edge belongs to a compound entry then
  // for all (boundary) vertices whithin the  edge visibility is set to
  // 0, if the edge doesn't belong to a  compound entry then visibility
  // is set to 1 for all its (boundary) vertices. Later, in FillSMesh()
  // func getVisibility() (returns either 1 or 0) is used to decide we-
  // ather the mesh of vertices should be transmitted  to SMESH or not.

  for ( GModel::eiter itE = _gModel->firstEdge(); itE != _gModel->lastEdge(); ++itE )
  {
    std::vector<GVertex *> bndVerticies = (*itE)->vertices();

    for( auto itV = bndVerticies.begin(); itV != bndVerticies.end(); ++itV )
    {
      if(((*itE)->compound.size()))
        (*itV)->setVisibility(0);
      else
        (*itV)->setVisibility(1);
    }
  }

}
#endif

//================================================================================
/*!
 * \brief Get a node by a GMSH mesh vertex
 */
//================================================================================

const SMDS_MeshNode* GMSHPlugin_Mesher::Node( const MVertex* v )
{
  std::map< const MVertex *, const SMDS_MeshNode* >::iterator v2n = _nodeMap.find( v );
  if ( v2n != _nodeMap.end() )
    return v2n->second;

  return nullptr;
}

//================================================================================
/*!
 * \brief Get a node by a GMSH mesh vertex
 */
//================================================================================

const SMDS_MeshNode* GMSHPlugin_Mesher::PremeshedNode( const MVertex* v )
{
  std::map< const MVertex *, const SMDS_MeshNode* >::iterator v2n = _premeshednodeMap.find( v );
  if ( v2n != _premeshednodeMap.end() )
    return v2n->second;

  return nullptr;
}


//================================================================================
/*!
 * \brief Return a corresponding sub-mesh if a shape is meshed
 */
//================================================================================

SMESHDS_SubMesh* GMSHPlugin_Mesher::HasSubMesh( const TopoDS_Shape& s )
{
  if ( SMESHDS_SubMesh*  sm = _mesh->GetMeshDS()->MeshElements( s ))
  {
    if ( s.ShapeType() == TopAbs_VERTEX )
      return ( sm->NbNodes() > 0 ) ? sm : nullptr;
    else
      return ( sm->NbElements() > 0 ) ? sm : nullptr;
  }
  return nullptr;
}

//================================================================================
/*!
 * \brief Write mesh from GModel instance to SMESH instance
 */
//================================================================================

void GMSHPlugin_Mesher::FillSMesh()
{
  SMESHDS_Mesh* meshDS = _mesh->GetMeshDS();

  // ADD 0D ELEMENTS
  for ( GModel::viter it = _gModel->firstVertex(); it != _gModel->lastVertex(); ++it)
  {
    GVertex *gVertex = *it;

    // GET topoVertex CORRESPONDING TO gVertex
    TopoDS_Vertex topoVertex = *((TopoDS_Vertex*)gVertex->getNativePtr());

    if (gVertex->getVisibility() == 0) // belongs to a compound
    {
      SMESH_subMesh* sm = _mesh->GetSubMesh(topoVertex);
      sm->SetIsAlwaysComputed(true); // prevent from displaying errors
      continue;
    }

    // FILL SMESH FOR topoVertex
    //nodes
    for( size_t i = 0; i < gVertex->mesh_vertices.size(); i++)
    {
      MVertex *v = gVertex->mesh_vertices[i];
      if(v->getIndex() >= 0)
      {
        if ( SMESHDS_SubMesh* sm = HasSubMesh( topoVertex ))
        {
          const SMDS_MeshNode *node = sm->GetNodes()->next();
          _nodeMap.insert({ v, node });
        }
        else
        {
          SMDS_MeshNode *node = meshDS->AddNode( v->x(),v->y(),v->z() );
          meshDS->SetNodeOnVertex( node, topoVertex );
          _nodeMap.insert({ v, node });
        }
      }
    }
    // WE DON'T ADD 0D ELEMENTS because it does not follow the salome meshers philosophy
    //elements
    // for(unsigned int i = 0; i < gVertex->getNumMeshElements(); i++)
    // {
    //   MElement *e = gVertex->getMeshElement(i);
    //   std::vector<MVertex*> verts;
    //   e->getVertices(verts);
    //   ASSERT(verts.size()==1);
    //   SMDS_Mesh0DElement* zeroDElement = 0;
    //   zeroDElement = meshDS->Add0DElementWithID(verts[0]->getNum(),e->getNum());
    //   meshDS->SetMeshElementOnShape(zeroDElement, topoVertex);
    // }
  }

  // ADD 1D ELEMENTS
  for(GModel::eiter it = _gModel->firstEdge(); it != _gModel->lastEdge(); ++it)
  {
    GEdge *gEdge = *it;

    // GET topoEdge CORRESPONDING TO gEdge
    TopoDS_Edge topoEdge;
    std::vector< ShapeBounds > topoEdges;
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if(gEdge->haveParametrization())
#else
    if ( gEdge->geomType() != GEntity::CompoundCurve )
#endif
    {
      topoEdge = *((TopoDS_Edge*)gEdge->getNativePtr());
      if (gEdge->getVisibility() == 0) // belongs to a compound
      {
        SMESH_subMesh* sm = _mesh->GetSubMesh(topoEdge);
        sm->SetIsAlwaysComputed(true); // prevent from displaying errors
        continue;
      }
      if ( HasSubMesh( topoEdge ))
        continue; // a meshed sub-mesh
    }
    bool isCompound = getBoundsOfShapes( gEdge, topoEdges );

    // FILL SMESH FOR topoEdge
    //nodes
    for ( size_t i = 0; i < gEdge->mesh_vertices.size(); i++ )
    {
      MVertex *v = gEdge->mesh_vertices[i];
      if ( v->getIndex() >= 0 )
      {
        SMDS_MeshNode *node = meshDS->AddNode( v->x(),v->y(),v->z() );

        if ( isCompound )
          topoEdge = TopoDS::Edge( getShapeAtPoint( v->point(), topoEdges ));

        // Based on BLSURFPlugin_BLSURF
        gp_Pnt point3D( v->x(),v->y(),v->z() );
        Standard_Real p0 = 0.0;
        Standard_Real p1 = 1.0;
        TopLoc_Location loc;
        Handle(Geom_Curve) curve = BRep_Tool::Curve(topoEdge, loc, p0, p1);

        if ( !curve.IsNull() )
        {
          if ( !loc.IsIdentity() )
            point3D.Transform( loc.Transformation().Inverted() );

          GeomAPI_ProjectPointOnCurve proj(point3D, curve, p0, p1);

          double pa = 0.;
          if ( proj.NbPoints() > 0 )
            pa = (double)proj.LowerDistanceParameter();

          meshDS->SetNodeOnEdge( node, topoEdge, pa );
        }
        else
        {
          meshDS->SetNodeOnEdge( node, topoEdge );
        }
        //END on BLSURFPlugin_BLSURF


        _nodeMap.insert({ v, node });
      }
    }
  }

  for ( GModel::eiter it = _gModel->firstEdge(); it != _gModel->lastEdge(); ++it )
  {
    GEdge *gEdge = *it;
    if ( gEdge->getVisibility() == 0) // belongs to a compound
      continue;

    TopoDS_Edge topoEdge;
    std::vector< ShapeBounds > topoEdges;
    bool isCompound = getBoundsOfShapes( gEdge, topoEdges );
    if ( !isCompound )
      topoEdge = *((TopoDS_Edge*)gEdge->getNativePtr());

    if ( HasSubMesh( topoEdge ))
      continue; // a meshed sub-mesh

    //elements
    std::vector<MVertex*> verts(3);
    for ( size_t i = 0; i < gEdge->getNumMeshElements(); i++ )
    {
      MElement *e = gEdge->getMeshElement(i);
      verts.clear();
      e->getVertices(verts);

      // if a node wasn't set, it is assigned here
      for ( size_t j = 0; j < verts.size(); j++ )
      {
        if ( verts[j]->onWhat()->getVisibility() == 0 )
        {
          SMDS_MeshNode *node = meshDS->AddNode(verts[j]->x(),verts[j]->y(),verts[j]->z() );

          gp_Pnt point3D( verts[j]->x(),verts[j]->y(),verts[j]->z() );
          Standard_Real p0 = 0.0;
          Standard_Real p1 = 1.0;
          TopLoc_Location loc;
          Handle(Geom_Curve) curve = BRep_Tool::Curve(topoEdge, loc, p0, p1);

          if ( !curve.IsNull() )
          {
            if ( !loc.IsIdentity() )
              point3D.Transform( loc.Transformation().Inverted() );

            GeomAPI_ProjectPointOnCurve proj(point3D, curve, p0, p1);

            double pa = 0.;
            if ( proj.NbPoints() > 0 )
              pa = (double)proj.LowerDistanceParameter();

            meshDS->SetNodeOnEdge( node, topoEdge, pa );
          }
          else
          {
            meshDS->SetNodeOnEdge( node, topoEdge );
          }

          verts[j]->setEntity(gEdge);
          _nodeMap.insert({ verts[j], node });
        }
      }

      SMDS_MeshEdge* edge = 0;
      switch (verts.size())
      {
        case 2:
          edge = meshDS->AddEdge(Node( verts[0]),
                                 Node( verts[1]));
          break;
        case 3:
          edge = meshDS->AddEdge(Node( verts[0]),
                                 Node( verts[1]),
                                 Node( verts[2]));
          break;
        default:
          ASSERT(false);
          continue;
      }
      if ( isCompound )
        topoEdge = TopoDS::Edge( getShapeAtPoint( e->barycenter(), topoEdges ));

      meshDS->SetMeshElementOnShape( edge, topoEdge );
    }
  }

  // ADD 2D ELEMENTS
  for ( GModel::fiter it = _gModel->firstFace(); it != _gModel->lastFace(); ++it)
  {
    GFace *gFace = *it;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    // Gmsh since version 4.3 is now producing extra surface and mesh when
    // compounds are involved. Since in Gmsh meshing procedure needs acess
    // to each of the original topology and the meshed topology. Hence  we
    // bypass the additional mesh in case of compounds. Note, similar cri-
    // teria also occurs in the following 'for' loop.
    if ( _compounds.size() && gFace->geomType() == GEntity::DiscreteSurface )
      continue;
#endif

    // GET topoFace CORRESPONDING TO gFace
    TopoDS_Face topoFace;
    std::vector< ShapeBounds > topoFaces;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if(gFace->haveParametrization())
#else
    if ( gFace->geomType() != GEntity::CompoundSurface )
#endif
    {
      topoFace = *((TopoDS_Face*)gFace->getNativePtr());
      if (gFace->getVisibility() == 0) // belongs to a compound
      {
        SMESH_subMesh* sm = _mesh->GetSubMesh(topoFace);
        sm->SetIsAlwaysComputed(true); // prevent from displaying errors
        continue;
      }
      if ( HasSubMesh( topoFace ))
        continue; // a meshed sub-mesh
    }
    bool isCompound = getBoundsOfShapes( gFace, topoFaces );

    // FILL SMESH FOR topoFace
    //nodes
    for ( size_t i = 0; i < gFace->mesh_vertices.size(); i++ )
    {
      MVertex *v = gFace->mesh_vertices[i];
      if ( v->getIndex() >= 0 )
      {
        double U,V;
        gFace->XYZtoUV( v->x(),v->y(),v->z(), U, V, 1.0 );

        SMDS_MeshNode *node = meshDS->AddNode( v->x(),v->y(),v->z() );

        if ( isCompound )
          topoFace = TopoDS::Face( getShapeAtPoint( v->point(), topoFaces ));

        meshDS->SetNodeOnFace( node, topoFace, U, V );
        _nodeMap.insert({ v, node });
      }
    }
  }

  for ( GModel::fiter it = _gModel->firstFace(); it != _gModel->lastFace(); ++it)
  {
    GFace *gFace = *it;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if ( _compounds.size() && gFace->geomType() == GEntity::DiscreteSurface )
      continue;

    bool isCompound = (!gFace->haveParametrization());
#else
    bool isCompound = ( gFace->geomType() == GEntity::CompoundSurface );
#endif

    if ( !isCompound && gFace->getVisibility() == 0 )
      continue;  // belongs to a compound

    TopoDS_Face topoFace;
    std::vector< ShapeBounds > topoFaces;
    if ( isCompound )
      getBoundsOfShapes( gFace, topoFaces );
    else
      topoFace = *((TopoDS_Face*)gFace->getNativePtr());

    if ( HasSubMesh( topoFace ))
      continue; // a meshed sub-mesh

    //elements
    std::vector<MVertex*> verts;
    for ( size_t i = 0; i < gFace->getNumMeshElements(); i++ )
    {
      MElement *e = gFace->getMeshElement(i);
      verts.clear();
      e->getVertices(verts);
      SMDS_MeshFace* face = 0;

      // if a node wasn't set, it is assigned here
      for ( size_t j = 0; j < verts.size(); j++)
      {
        if(verts[j]->onWhat()->getVisibility() == 0)
        {
          SMDS_MeshNode *node = meshDS->AddNode(verts[j]->x(),verts[j]->y(),verts[j]->z());
          double U,V;
          gFace->XYZtoUV( verts[j]->x(),verts[j]->y(),verts[j]->z(), U, V, 1.0 );
          meshDS->SetNodeOnFace( node, topoFace, U, V );
          _nodeMap.insert({ verts[j], node });
          verts[j]->setEntity(gFace);
        }
      }
      switch (verts.size())
      {
        case 3:
          face = meshDS->AddFace(Node( verts[0]),
                                 Node( verts[1]),
                                 Node( verts[2]));
          break;
        case 4:
          face = meshDS->AddFace(Node( verts[0]),
                                 Node( verts[1]),
                                 Node( verts[2]),
                                 Node( verts[3]));
          break;
        case 6:
          face = meshDS->AddFace(Node( verts[0]),
                                 Node( verts[1]),
                                 Node( verts[2]),
                                 Node( verts[3]),
                                 Node( verts[4]),
                                 Node( verts[5]));
          break;
        case 8:
          face = meshDS->AddFace(Node( verts[0]),
                                 Node( verts[1]),
                                 Node( verts[2]),
                                 Node( verts[3]),
                                 Node( verts[4]),
                                 Node( verts[5]),
                                 Node( verts[6]),
                                 Node( verts[7]));
          break;
        case 9:
          face = meshDS->AddFace(Node( verts[0]),
                                 Node( verts[1]),
                                 Node( verts[2]),
                                 Node( verts[3]),
                                 Node( verts[4]),
                                 Node( verts[5]),
                                 Node( verts[6]),
                                 Node( verts[7]),
                                 Node( verts[8]));
          break;
        default:
          ASSERT(false);
          continue;
      }

      if ( isCompound )
        topoFace = TopoDS::Face( getShapeAtPoint( e->barycenter(), topoFaces ));

      meshDS->SetMeshElementOnShape(face, topoFace);
    }
  }

  // ADD 3D ELEMENTS
  for ( GModel::riter it = _gModel->firstRegion(); it != _gModel->lastRegion(); ++it)
  {
    GRegion *gRegion = *it;
    if (gRegion->getVisibility() == 0)
      continue;

    // GET topoSolid CORRESPONDING TO gRegion
    TopoDS_Solid topoSolid = *((TopoDS_Solid*)gRegion->getNativePtr());

    // FILL SMESH FOR topoSolid

    //nodes
    for( size_t i = 0; i < gRegion->mesh_vertices.size(); i++)
    {
      MVertex *v = gRegion->mesh_vertices[i];
      if(v->getIndex() >= 0)
      {
        SMDS_MeshNode *node = meshDS->AddNode( v->x(),v->y(),v->z() );
        meshDS->SetNodeInVolume( node, topoSolid );
        _nodeMap.insert({ v, node });
      }
    }

    //elements
    std::vector<MVertex*> verts;
    for( size_t i = 0; i < gRegion->getNumMeshElements(); i++)
    {
      MElement *e = gRegion->getMeshElement(i);
      verts.clear();
      e->getVertices(verts);
      SMDS_MeshVolume* volume = 0;
      switch (verts.size()){
      case 4:
        volume = meshDS->AddVolume(Node( verts[0]),
                                   Node( verts[2]),
                                   Node( verts[1]),
                                   Node( verts[3]));
        break;
      case 5:
        volume = meshDS->AddVolume(Node( verts[0]),
                                   Node( verts[3]),
                                   Node( verts[2]),
                                   Node( verts[1]),
                                   Node( verts[4]));
        break;
      case 6:
        volume = meshDS->AddVolume(Node( verts[0]),
                                   Node( verts[2]),
                                   Node( verts[1]),
                                   Node( verts[3]),
                                   Node( verts[5]),
                                   Node( verts[4]));
        break;
      case 8:
        volume = meshDS->AddVolume(Node( verts[0]),
                                   Node( verts[3]),
                                   Node( verts[2]),
                                   Node( verts[1]),
                                   Node( verts[4]),
                                   Node( verts[7]),
                                   Node( verts[6]),
                                   Node( verts[5]));
        break;
      case 10:
        volume = meshDS->AddVolume(Node( verts[0]),
                                   Node( verts[2]),
                                   Node( verts[1]),
                                   Node( verts[3]),
                                   Node( verts[6]),
                                   Node( verts[5]),
                                   Node( verts[4]),
                                   Node( verts[7]),
                                   Node( verts[8]),
                                   Node( verts[9]));
        break;
      case 13:
        volume = meshDS->AddVolume(Node( verts[0] ),
                                   Node( verts[3] ),
                                   Node( verts[2] ),
                                   Node( verts[1] ),
                                   Node( verts[4] ),
                                   Node( verts[6] ),
                                   Node( verts[10] ),
                                   Node( verts[8] ),
                                   Node( verts[5] ),
                                   Node( verts[7] ),
                                   Node( verts[12] ),
                                   Node( verts[11] ),
                                   Node( verts[9]));
        break;
      case 14: // same as case 13, because no pyra14 in smesh
        volume = meshDS->AddVolume(Node( verts[0] ),
                                   Node( verts[3] ),
                                   Node( verts[2] ),
                                   Node( verts[1] ),
                                   Node( verts[4] ),
                                   Node( verts[6] ),
                                   Node( verts[10] ),
                                   Node( verts[8] ),
                                   Node( verts[5] ),
                                   Node( verts[7] ),
                                   Node( verts[12] ),
                                   Node( verts[11] ),
                                   Node( verts[9]));
        break;
      case 15:
        volume = meshDS->AddVolume(Node( verts[0] ),
                                   Node( verts[2] ),
                                   Node( verts[1] ),
                                   Node( verts[3] ),
                                   Node( verts[5] ),
                                   Node( verts[4] ),
                                   Node( verts[7] ),
                                   Node( verts[9] ),
                                   Node( verts[6] ),
                                   Node( verts[13] ),
                                   Node( verts[14] ),
                                   Node( verts[12] ),
                                   Node( verts[8] ),
                                   Node( verts[11] ),
                                   Node( verts[10]));
        break;
      case 18: // same as case 15, because no penta18 in smesh
        volume = meshDS->AddVolume(Node( verts[0] ),
                                   Node( verts[2] ),
                                   Node( verts[1] ),
                                   Node( verts[3] ),
                                   Node( verts[5] ),
                                   Node( verts[4] ),
                                   Node( verts[7] ),
                                   Node( verts[9] ),
                                   Node( verts[6] ),
                                   Node( verts[13] ),
                                   Node( verts[14] ),
                                   Node( verts[12] ),
                                   Node( verts[8] ),
                                   Node( verts[11] ),
                                   Node( verts[10]));
        break;
      case 20:
        volume = meshDS->AddVolume(Node( verts[0] ),
                                   Node( verts[3] ),
                                   Node( verts[2] ),
                                   Node( verts[1] ),
                                   Node( verts[4] ),
                                   Node( verts[7] ),
                                   Node( verts[6] ),
                                   Node( verts[5] ),
                                   Node( verts[9] ),
                                   Node( verts[13] ),
                                   Node( verts[11] ),
                                   Node( verts[8] ),
                                   Node( verts[17] ),
                                   Node( verts[19] ),
                                   Node( verts[18] ),
                                   Node( verts[16] ),
                                   Node( verts[10] ),
                                   Node( verts[15] ),
                                   Node( verts[14] ),
                                   Node( verts[12]));
        break;
      case 27:
        volume = meshDS->AddVolume(Node( verts[0] ),
                                   Node( verts[3] ),
                                   Node( verts[2] ),
                                   Node( verts[1] ),
                                   Node( verts[4] ),
                                   Node( verts[7] ),
                                   Node( verts[6] ),
                                   Node( verts[5] ),
                                   Node( verts[9] ),
                                   Node( verts[13] ),
                                   Node( verts[11] ),
                                   Node( verts[8] ),
                                   Node( verts[17] ),
                                   Node( verts[19] ),
                                   Node( verts[18] ),
                                   Node( verts[16] ),
                                   Node( verts[10] ),
                                   Node( verts[15] ),
                                   Node( verts[14] ),
                                   Node( verts[12] ),
                                   Node( verts[20] ),
                                   Node( verts[22] ),
                                   Node( verts[24] ),
                                   Node( verts[23] ),
                                   Node( verts[21] ),
                                   Node( verts[25] ),
                                   Node( verts[26] ));
        break;
      default:
        ASSERT(false);
        continue;
      }
      meshDS->SetMeshElementOnShape(volume, topoSolid);
    }
  }

  //return 0;
}

//================================================================================
/*!
 * \brief Find if SPoint point is in SBoundingBox3d bounds
 */
//================================================================================

float GMSHPlugin_Mesher::DistBoundingBox(const SBoundingBox3d& bounds, const SPoint3& point)
{
  SPoint3 min = bounds.min();
  SPoint3 max = bounds.max();

  float x,y,z;

  if (point.x() < min.x())
    x = min.x()-point.x();
  else if (point.x() > max.x())
    x = point.x()-max.x();
  else
    x = 0.;

  if (point.y() < min.y())
    y = min.y()-point.y();
  else if (point.y() > max.y())
    y = point.y()-max.y();
  else
    y = 0.;

  if (point.z() < min.z())
    z = min.z()-point.z();
  else if (point.z() > max.z())
    z = point.z()-max.z();
  else
    z = 0.;

  return x*x+y*y+z*z;
}
//================================================================================
/*!
 * \brief Reimplemented GmshMessage call. Actions done if errors occurs
 *        during gmsh meshing. We define here what to display and what to do.
 */
//================================================================================
void  GMSHPlugin_Mesher::mymsg::operator()(std::string level, std::string msg)
{
  //MESSAGE("level="<< level.c_str() << ", msg=" << msg.c_str()<< "\n");
  printf("level=%s msg=%s\n", level.c_str(), msg.c_str());

  if(level == "Fatal" || level == "Error")
  {
    std::ostringstream oss;
    if (level == "Fatal")
      oss << "Fatal error during Generation of Gmsh Mesh\n";
    else
      oss << "Error during Generation of Gmsh Mesh\n";
    oss << "  " << msg.c_str() << "\n";
    GEntity *e = _gModel->getCurrentMeshEntity();
    if(e)
    {
      oss << "  error occurred while meshing entity:\n" <<
             "    tag=" << e->tag() << "\n" <<
             "    dimension=" << e->dim() << "\n" <<
             "    native pointer=" << e->getNativePtr();
      //if(e->geomType() != GEntity::CompoundCurve and e->geomType() != GEntity::CompoundSurface)
      //{
        //SMESH_subMesh *sm = _mesh->GetSubMesh(*((TopoDS_Shape*)e->getNativePtr()));
        //SMESH_ComputeErrorPtr& smError = sm->GetComputeError();
        //SMESH_Comment comment;
        //comment << SMESH_Comment(oss.str);
        //std::string str = oss.str();
        //smError.reset( new SMESH_ComputeError( str ));

        // plutot que de faire de la merde ici, on pourait simplement
        // remplir une liste pour dire sur quelles entités gmsh se plante
        // (puis faire le fillsmesh)
        // puis faire une nouvelle routine qui réécrit les messages d'erreur
        // probleme : gmsh peut planté en Fatal, dans ce cas pas de fillsmesh
      //}
    }
    if (level == "Fatal")
    {
        CTX::instance()->lock = 0;
        throw oss.str();
    }
    else
        printf("%s\n", oss.str().c_str());
  }
}

bool GMSHPlugin_Mesher::Compute3D( std::vector< const SMDS_MeshNode* >& nodeVec, std::map<const SMDS_MeshElement*, bool, TIDCompare>& listElements, bool addElements )
{
  MESSAGE("GMSHPlugin_Mesher::Compute3D");
  int err = 0;
  _maxThreads = 1;
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=3
  _maxThreads = 1;
#endif

  char* argv[] = {"-noenv"};
  GmshInitialize(1,argv);
  SetGmshOptions();
  _gModel = new GModel();
  mymsg msg(_gModel);
  GmshSetMessageHandler(&msg);

  _gModel->importOCCShape((void*)&_shape);
  try
  {
    HideComputedEntities( _gModel, true );   
    // fill geometry with elements as geom objects
    FillGeomMapMeshUsing2DMeshIterator( listElements );
    Set2DMeshes( nodeVec, listElements );
    _gModel->mesh( /*dim=*/ 3);
  }
  catch (std::string& str)
  {
    err = 1;
    MESSAGE(str);
  }
  catch (...)
  {
    err = 1;
    MESSAGE("Unrecoverable error during Generation of Gmsh Mesh");
  }

  if (!err)
  {
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if (_compounds.size() > 0)
      SetCompoundMeshVisibility();
#endif
  }
  
  if ( addElements )
    FillSMesh();
  MESSAGE("GMSHPlugin_Mesher::Compute3D:End");
  return err;
}

void GMSHPlugin_Mesher::finalizeGModel()
{
  if ( _gModel )
  {
    GmshSetMessageHandler(nullptr);
    delete _gModel;
    GmshFinalize();
  }
}

//=============================================================================
/*!
 * Here we are going to use the GMSH mesher
 */
//=============================================================================

bool GMSHPlugin_Mesher::Compute()
{
  MESSAGE("GMSHPlugin_Mesher::Compute");

  int err = 0;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
  SetMaxThreadsGmsh();
#endif
  //RNV: to avoid modification of PATH and PYTHONPATH
  char* argv[] = {"-noenv"};
  GmshInitialize(1,argv);
  SetGmshOptions();
  _gModel = new GModel();
  mymsg msg(_gModel);
  GmshSetMessageHandler(&msg);
  _gModel->importOCCShape((void*)&_shape);
  if (_compounds.size() > 0) CreateGmshCompounds();
  try
  {

    HideComputedEntities(_gModel);
    if (_is3d)
    {
      FillGMSHMesh();
      Set2DSubMeshes(_gModel);
      _gModel->mesh( /*dim=*/ 3);
    }
    else
    {
      //CTX::instance()->mesh.maxNumThreads1D=1;
      _gModel->mesh( /*dim=*/ 1);

      Set1DSubMeshes(_gModel);

      //_gModel->writeUNV("/tmp/1D.unv", 1,0,0);
      //CTX::instance()->mesh.maxNumThreads2D=1;

      _gModel->mesh( /*dim=*/ 2);

      if (!_is2d)
      {
        Set2DSubMeshes(_gModel);

        //CTX::instance()->mesh.maxNumThreads3D=1;

        _gModel->mesh( /*dim=*/ 3);
      }
      RestoreVisibility(_gModel);
    }
#ifdef WITH_SMESH_CANCEL_COMPUTE

#endif
  }
  catch (std::string& str)
  {
    err = 1;
    std::cerr << "GMSH: exception caught: " << str << std::endl;
    MESSAGE(str);
  }
  catch (...)
  {
    err = 1;
    std::cerr << "GMSH: Unknown exception caught: " << std::endl;
    MESSAGE("Unrecoverable error during Generation of Gmsh Mesh");
  }

  if (!err)
  {
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if (_compounds.size() > 0)
      SetCompoundMeshVisibility();
#endif
    FillSMesh();
  }
  GmshSetMessageHandler(nullptr);
  delete _gModel;
  GmshFinalize();
  MESSAGE("GMSHPlugin_Mesher::Compute:End");
  return !err;
}

//================================================================================
/*!
 * \brief Set 1D sub-meshes to GModel. GMSH 1D mesh is made by now.
 *  \param [inout] _gModel - GMSH model
 */
//================================================================================

void GMSHPlugin_Mesher::Set1DSubMeshes( GModel* gModel )
{
  SMESHDS_Mesh* meshDS = _mesh->GetMeshDS();

  for(GModel::eiter it = gModel->firstEdge(); it != gModel->lastEdge(); ++it)
  {
    GEdge *gEdge = *it;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if ( !gEdge->haveParametrization())
#else
    if ( gEdge->geomType() == GEntity::CompoundCurve )
#endif
      continue;

    TopoDS_Edge topoEdge = *((TopoDS_Edge*)gEdge->getNativePtr());
    if ( !HasSubMesh( topoEdge ))
      continue; // empty sub-mesh

    gEdge->deleteMesh();

    // get node parameters on topoEdge
    StdMeshers_FaceSide side( TopoDS_Face(), topoEdge, _mesh, /*fwd=*/true, /*skpMedium=*/true);
    const UVPtStructVec& nodeParam = side.GetUVPtStruct();
    if ( nodeParam.empty() )
      throw std::string("Pb with StdMeshers_FaceSide::GetUVPtStruct()");

    // get GMSH mesh vertices on VERTEX'es
    std::vector<MVertex *> mVertices( nodeParam.size(), nullptr );
    GVertex * gV0 = gEdge->getBeginVertex(), *gV1 = gEdge->getEndVertex();
    mVertices[0]     = gV0->mesh_vertices[ 0 ];
    mVertices.back() = gV1->mesh_vertices[ 0 ];
    TopoDS_Vertex v01 = *((TopoDS_Vertex*) gV0->getNativePtr());
    TopoDS_Shape  v02 = SMESH_MesherHelper::GetSubShapeByNode( nodeParam[0].node, meshDS );
    bool      reverse = !v01.IsSame( v02 );
    if ( mVertices[0] == mVertices.back() )
      reverse = ( nodeParam[0].param > nodeParam.back().param );
    const SMDS_MeshNode* n0 = reverse ? nodeParam.back().node : nodeParam[0].node;
    const SMDS_MeshNode* n1 = reverse ? nodeParam[0].node : nodeParam.back().node;
    _nodeMap.insert({ mVertices[ 0 ],   n0 });
    _nodeMap.insert({ mVertices.back(), n1 });

    // create GMSH mesh vertices on gEdge
    for ( size_t i = 1; i < nodeParam.size() - 1; ++i )
    {
      size_t iN = reverse ? ( nodeParam.size() - 1 - i ) : i;
      SMESH_NodeXYZ xyz = nodeParam[ iN ].node;
      double lc = segmentSize( nodeParam, iN );
      // SVector3 der = gEdge->firstDer(nodeParam[ iN ].param);
      // double lc = norm(der) / segmentSize( nodeParam, i );

      mVertices[ i ] = new MEdgeVertex( xyz.X(), xyz.Y(), xyz.Z(),
                                        gEdge, nodeParam[ iN ].param, 0, lc);
      gEdge->mesh_vertices.push_back( mVertices[ i ]);
      _nodeMap.insert({ mVertices[ i ], nodeParam[ iN ].node });
    }
    // create GMSH mesh edges
    for ( size_t i = 1; i < mVertices.size(); ++i )
    {
      gEdge->lines.push_back( new MLine( mVertices[ i - 1 ],
                                         mVertices[ i ]));
    }
    /*{
      cout << endl << "EDGE " << gEdge->tag() <<
        ( topoEdge.Orientation() == TopAbs_FORWARD ? " F" : " R") << endl;
      MVertex* mv = gV0->mesh_vertices[ 0 ];
      cout << "V0: " << mv->x() << ", " << mv->y() << ", " << mv->z() << ", "<<endl;
      for ( size_t i = 0; i < gEdge->mesh_vertices.size(); ++i )
      {
        MEdgeVertex* mv = (MEdgeVertex*) gEdge->mesh_vertices[i];
        cout << i << ": " << mv->x() << ", " << mv->y() << ", " << mv->z() << ", ";
        double t;
        mv->getParameter(0, t );
        cout << ":\t t = "<< t << " lc = " << mv->getLc() << endl;
      }
      mv = gV1->mesh_vertices[ 0 ];
      cout << "V1: " << mv->x() << ", " << mv->y() << ", " << mv->z() << ", "<<endl;
      }*/
  }
  return;
}

//================================================================================
/*!
 * \brief Set 2D sub-meshes to GModel. GMSH 2D mesh is made by now.
 *  \param [inout] _gModel - GMSH model
 */
//================================================================================

void GMSHPlugin_Mesher::Set2DSubMeshes( GModel* gModel )
{
  if ( _nodeMap.empty() )
    return; // no sub-meshes

  SMESH_MesherHelper helper( *_mesh );

  std::map< const SMDS_MeshNode* , const MVertex * > nodes2mvertMap;
  for ( auto & v2n : _nodeMap )
    nodes2mvertMap.insert({ v2n.second, v2n.first });

  std::vector<MVertex *> mVertices;

  for(GModel::fiter it = gModel->firstFace(); it != gModel->lastFace(); ++it)
  {
    GFace *gFace = *it;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if ( !gFace->haveParametrization())
#else
      if ( gFace->geomType() == GEntity::CompoundSurface )
#endif
        continue;

    TopoDS_Face topoFace = *((TopoDS_Face*)gFace->getNativePtr());
    SMESHDS_SubMesh*  sm = HasSubMesh( topoFace );
    if ( !sm )
      continue;
    //_gModel->writeUNV("/tmp/befDEl.unv", 1,0,0);

    gFace->deleteMesh();

    bool reverse = false;
    if ( gFace->getRegion(0) )
    {
      //GRegion * gRegion = gFace->getRegion(0);
      TopoDS_Shape topoSolid = *((TopoDS_Shape*)gFace->getNativePtr());
      TopAbs_Orientation faceOriInSolid = helper.GetSubShapeOri( topoSolid, topoFace );
      if ( faceOriInSolid >= 0 )
        reverse =
          helper.IsReversedSubMesh( TopoDS::Face( topoFace.Oriented( faceOriInSolid )));
    }

    for ( SMDS_ElemIteratorPtr fIt = sm->GetElements(); fIt->more(); )
    {
      const SMDS_MeshElement* f = fIt->next();

      int nbN = f->NbCornerNodes();
      if ( nbN > 4 )
        throw std::string("Polygon sub-meshes not supported");

      mVertices.resize( nbN );
      for ( int i = 0; i < nbN; ++i )
      {
        const SMDS_MeshNode* n = f->GetNode( i );
        MVertex *           mv = nullptr;
        auto n2v = nodes2mvertMap.find( n );
        if ( n2v != nodes2mvertMap.end() )
        {
          mv = const_cast< MVertex*>( n2v->second );
        }
        else
        {
          if ( n->GetPosition()->GetDim() < 2 )
            throw std::string("Wrong mapping of edge nodes to GMSH nodes");
          SMESH_NodeXYZ xyz = n;
          bool ok = true;
          gp_XY uv = helper.GetNodeUV( topoFace, n, nullptr, &ok );
          mv = new MFaceVertex( xyz.X(), xyz.Y(), xyz.Z(), gFace, uv.X(), uv.Y() );
          gFace->mesh_vertices.push_back( mv );
          nodes2mvertMap.insert({ n, mv });
          _nodeMap.insert      ({ mv, n });
        }
        mVertices[ i ] = mv;
      }
      // create GMSH mesh faces
      switch ( nbN ) {
      case 3:
        if ( reverse )
          gFace->triangles.push_back (new MTriangle(mVertices[0], mVertices[2], mVertices[1]));
        else
          gFace->triangles.push_back (new MTriangle(mVertices[0], mVertices[1], mVertices[2]));
        break;
      case 4:
        if ( reverse )
          gFace->quadrangles.push_back (new MQuadrangle(mVertices[0], mVertices[3],
                                                        mVertices[2], mVertices[1]));
        else
          gFace->quadrangles.push_back (new MQuadrangle(mVertices[0], mVertices[1],
                                                        mVertices[2], mVertices[3]));
        break;
      default:;
      }
    }
  } // loop on GMSH faces

  return;
}

//================================================================================
/*!
 * \brief Set visibility 0 to already computed geom entities
 *        to prevent their meshing
 */
//================================================================================

void GMSHPlugin_Mesher::HideComputedEntities( GModel* gModel, bool hideAnyway )
{
  CTX::instance()->mesh.meshOnlyVisible = true;

  for(GModel::eiter it = gModel->firstEdge(); it != gModel->lastEdge(); ++it)
  {
    GEdge *gEdge = *it;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if ( !gEdge->haveParametrization())
#else
      if ( gEdge->geomType() == GEntity::CompoundCurve )
#endif
        continue;

    TopoDS_Edge topoEdge = *((TopoDS_Edge*)gEdge->getNativePtr());
    if ( HasSubMesh( topoEdge ) || hideAnyway )
      gEdge->setVisibility(0);
  }


  for(GModel::fiter it = gModel->firstFace(); it != gModel->lastFace(); ++it)
  {
    GFace *gFace = *it;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if ( !gFace->haveParametrization())
#else
      if ( gFace->geomType() == GEntity::CompoundSurface )
#endif
        continue;

    TopoDS_Face topoFace = *((TopoDS_Face*)gFace->getNativePtr());
    if ( HasSubMesh( topoFace ) || hideAnyway )
      gFace->setVisibility(0);
  }
}

//================================================================================
/*!
 * \brief Restore visibility of all geom entities
 */
//================================================================================

void GMSHPlugin_Mesher::RestoreVisibility( GModel* gModel )
{
  for(GModel::eiter it = gModel->firstEdge(); it != gModel->lastEdge(); ++it)
  {
    GEdge *gEdge = *it;
    gEdge->setVisibility(1);
  }
  for(GModel::fiter it = gModel->firstFace(); it != gModel->lastFace(); ++it)
  {
    GFace *gFace = *it;
    gFace->setVisibility(1);
  }
}

/*
  void GMSHPlugin_Mesher::toPython( GModel* )
  {
  const char*  pyFile = "/tmp/gMesh.py";
  ofstream outfile( pyFile, ios::out );
  if ( !outfile ) return;

  outfile << "import salome, SMESH" << std::endl
          << "from salome.smesh import smeshBuilder" << std::endl
          << "smesh = smeshBuilder.New()" << std::endl
          << "mesh = smesh.Mesh()" << std::endl << std::endl;

  outfile << "## VERTICES" << endl;
  for ( GModel::viter it = _gModel->firstVertex(); it != _gModel->lastVertex(); ++it)
  {
    GVertex *gVertex = *it;

    for(unsigned int i = 0; i < gVertex->mesh_vertices.size(); i++)
    {
      MVertex *v = gVertex->mesh_vertices[i];
      if ( v->getIndex() >= 0)
      {
        outfile << "n" << v->getNum() << " = mesh.AddNode("
                << v->x() << ", " << v->y() << ", " << v->z()<< ")"
                << " ## tag = " << gVertex->tag() << endl;
      }
    }
  }

  for(GModel::eiter it = _gModel->firstEdge(); it != _gModel->lastEdge(); ++it)
  {
    GEdge *gEdge = *it;
    outfile << "## GEdge " << gEdge->tag() << endl;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if(gEdge->haveParametrization())
#else
      if ( gEdge->geomType() != GEntity::CompoundCurve )
#endif
        for ( size_t i = 0; i < gEdge->mesh_vertices.size(); i++ )
        {
          MVertex *v = gEdge->mesh_vertices[i];
          if ( v->getIndex() >= 0 )
          {
            outfile << "n" << v->getNum() << " = mesh.AddNode("
                    << v->x() << ", " << v->y() << ", " << v->z()<< ")"
                    << " ## tag = " << gEdge->tag() << endl;
          }
        }
  }

  for ( GModel::fiter it = _gModel->firstFace(); it != _gModel->lastFace(); ++it)
  {
    GFace *gFace = *it;
    if ( _compounds.size() && gFace->geomType() == GEntity::DiscreteSurface )
      continue;
    outfile << "## GFace " << gFace->tag() << endl;

    for ( size_t i = 0; i < gFace->mesh_vertices.size(); i++ )
    {
      MVertex *v = gFace->mesh_vertices[i];
      if ( v->getIndex() >= 0 )
      {
        outfile << "n" << v->getNum() << " = mesh.AddNode("
                << v->x() << ", " << v->y() << ", " << v->z()<< ")"
                << " ## tag = " << gFace->tag() << endl;
      }
    }
  }

  std::vector<MVertex*> verts(3);
  for(GModel::eiter it = _gModel->firstEdge(); it != _gModel->lastEdge(); ++it)
  {
    GEdge *gEdge = *it;
    outfile << "## GEdge " << gEdge->tag() << endl;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    if(gEdge->haveParametrization())
#else
      if ( gEdge->geomType() != GEntity::CompoundCurve )
#endif
    for ( size_t i = 0; i < gEdge->getNumMeshElements(); i++ )
    {
      MElement *e = gEdge->getMeshElement(i);
      verts.clear();
      e->getVertices(verts);

      outfile << "e" << e->getNum() << " = mesh.AddEdge(["
              << "n" << verts[0]->getNum() << ","
              << "n" << verts[1]->getNum();
      if ( verts.size() == 3 )
        outfile << "n" << verts[2]->getNum();
      outfile << "])"<< endl;
    }
  }

  for ( GModel::fiter it = _gModel->firstFace(); it != _gModel->lastFace(); ++it)
  {
    GFace *gFace = *it;
    if ( _compounds.size() && gFace->geomType() == GEntity::DiscreteSurface )
      continue;
    outfile << "## GFace " << gFace->tag() << endl;

    for ( size_t i = 0; i < gFace->getNumMeshElements(); i++ )
    {
      MElement *e = gFace->getMeshElement(i);
      verts.clear();
      e->getVertices(verts);

      outfile << "f" << e->getNum() << " = mesh.AddFace([";
      for ( size_t j = 0; j < verts.size(); j++)
      {
        outfile << "n" << verts[j]->getNum();
        if ( j < verts.size()-1 )
          outfile << ", ";
      }
      outfile << "])" << endl;
    }
  }
  std::cout << "Write " << pyFile << std::endl;
}
*/
