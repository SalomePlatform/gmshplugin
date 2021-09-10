// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016-2021  EDF R&D
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

#include <SMDS_FaceOfNodes.hxx>
#include <SMDS_MeshElement.hxx>
#include <SMDS_MeshNode.hxx>
#include <SMESHDS_Mesh.hxx>
#include <SMESH_Block.hxx>
#include <SMESH_Comment.hxx>
#include <SMESH_ComputeError.hxx>
#include <SMESH_File.hxx>
#include <SMESH_Gen_i.hxx>
#include <SMESH_Mesh.hxx>
#include <SMESH_MesherHelper.hxx>
#include <SMESH_subMesh.hxx>
#include <utilities.h>

#include <vector>
#include <limits>

#include <BRep_Tool.hxx>
#include <Bnd_B3d.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <NCollection_Map.hxx>
#include <OSD_File.hxx>
#include <OSD_Path.hxx>
#include <Standard_ErrorHandler.hxx>
#include <Standard_ProgramError.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_DataMapIteratorOfDataMapOfShapeInteger.hxx>
#include <TopTools_DataMapIteratorOfDataMapOfShapeShape.hxx>
#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_DataMapOfShapeShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>

#if GMSH_MAJOR_VERSION >=4
#include <GmshGlobal.h>
#include <Context.h>
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
}

//=============================================================================
/*!
 *
 */
//=============================================================================

GMSHPlugin_Mesher::GMSHPlugin_Mesher (SMESH_Mesh* mesh,
                                      const TopoDS_Shape& aShape)
  : _mesh    (mesh),
    _shape   (aShape)
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
    _minSize         = hyp->GetMinSize();
    _maxSize         = hyp->GetMaxSize();
    _secondOrder     = hyp->GetSecondOrder();
    _useIncomplElem  = hyp->GetUseIncomplElem();
    _is2d            = hyp->GetIs2d();
    _compounds       = hyp->GetCompoundOnEntries();
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
    _minSize         = 0;
    _maxSize         = 1e22;
    _secondOrder     = false;
    _useIncomplElem  = true;
    _is2d            = false;
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
  if (_compounds.size() > 0)
    _maxThreads = 1;
  else
    _maxThreads = omp_get_max_threads();
}
#endif

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

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
/*ok = GmshSetOption("Mesh", "MaxNumThreads1D"          , 0. )  ; // Coarse-grain algo threads
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "MaxNumThreads2D"          , 0. )  ; // Coarse-grain algo threads
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "MaxNumThreads3D"          , 0. )  ; // Fine-grain algo threads HXT
  ASSERT(ok);
/**/
  ok = GmshSetOption("General", "NumThreads"            , _maxThreads )  ; // system default i.e. OMP_NUM_THREADS
  ASSERT(ok);
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
    for(unsigned int i = 0; i < gVertex->mesh_vertices.size(); i++)
    {
      MVertex *v = gVertex->mesh_vertices[i];
      if(v->getIndex() >= 0)
      {
        SMDS_MeshNode *node = meshDS->AddNodeWithID(v->x(),v->y(),v->z(),v->getNum());
        meshDS->SetNodeOnVertex( node, topoVertex );
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
    }
    bool isCompound = getBoundsOfShapes( gEdge, topoEdges );

    // FILL SMESH FOR topoEdge
    //nodes
    for ( size_t i = 0; i < gEdge->mesh_vertices.size(); i++ )
    {
      MVertex *v = gEdge->mesh_vertices[i];
      if ( v->getIndex() >= 0 )
      {
        SMDS_MeshNode *node = meshDS->AddNodeWithID(v->x(),v->y(),v->z(),v->getNum());

        if ( isCompound )
          topoEdge = TopoDS::Edge( getShapeAtPoint( v->point(), topoEdges ));

        meshDS->SetNodeOnEdge( node, topoEdge );
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
          SMDS_MeshNode *node = meshDS->AddNodeWithID(verts[i]->x(),verts[j]->y(),verts[j]->z(),verts[j]->getNum());
          meshDS->SetNodeOnEdge( node, topoEdge );
          verts[j]->setEntity(gEdge);
        }
      }

      SMDS_MeshEdge* edge = 0;
      switch (verts.size())
      {
        case 2:
          edge = meshDS->AddEdgeWithID(verts[0]->getNum(),
                                       verts[1]->getNum(),e->getNum());
          break;
        case 3:
          edge = meshDS->AddEdgeWithID(verts[0]->getNum(),
                                       verts[1]->getNum(),
                                       verts[2]->getNum(),e->getNum());
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
    // teria also occus in the following 'for' loop.
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
    }
    bool isCompound = getBoundsOfShapes( gFace, topoFaces );

    // FILL SMESH FOR topoFace
    //nodes
    for ( size_t i = 0; i < gFace->mesh_vertices.size(); i++ )
    {
      MVertex *v = gFace->mesh_vertices[i];
      if ( v->getIndex() >= 0 )
      {
        SMDS_MeshNode *node = meshDS->AddNodeWithID(v->x(),v->y(),v->z(),v->getNum());

        if ( isCompound )
          topoFace = TopoDS::Face( getShapeAtPoint( v->point(), topoFaces ));

        meshDS->SetNodeOnFace( node, topoFace );
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
          SMDS_MeshNode *node = meshDS->AddNodeWithID(verts[j]->x(),verts[j]->y(),verts[j]->z(),verts[j]->getNum());
          meshDS->SetNodeOnFace( node, topoFace );
          verts[i]->setEntity(gFace);
        }
      }
      switch (verts.size())
      {
        case 3:
          face = meshDS->AddFaceWithID(verts[0]->getNum(),
                                       verts[1]->getNum(),
                                       verts[2]->getNum(),e->getNum());
          break;
        case 4:
          face = meshDS->AddFaceWithID(verts[0]->getNum(),
                                       verts[1]->getNum(),
                                       verts[2]->getNum(),
                                       verts[3]->getNum(),e->getNum());
          break;
        case 6:
          face = meshDS->AddFaceWithID(verts[0]->getNum(),
                                       verts[1]->getNum(),
                                       verts[2]->getNum(),
                                       verts[3]->getNum(),
                                       verts[4]->getNum(),
                                       verts[5]->getNum(),e->getNum());
          break;
        case 8:
          face = meshDS->AddFaceWithID(verts[0]->getNum(),
                                       verts[1]->getNum(),
                                       verts[2]->getNum(),
                                       verts[3]->getNum(),
                                       verts[4]->getNum(),
                                       verts[5]->getNum(),
                                       verts[6]->getNum(),
                                       verts[7]->getNum(),e->getNum());
          break;
        case 9:
          face = meshDS->AddFaceWithID(verts[0]->getNum(),
                                       verts[1]->getNum(),
                                       verts[2]->getNum(),
                                       verts[3]->getNum(),
                                       verts[4]->getNum(),
                                       verts[5]->getNum(),
                                       verts[6]->getNum(),
                                       verts[7]->getNum(),
                                       verts[8]->getNum(),e->getNum());
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
    for(unsigned int i = 0; i < gRegion->mesh_vertices.size(); i++)
    {
      MVertex *v = gRegion->mesh_vertices[i];
      if(v->getIndex() >= 0)
      {
        SMDS_MeshNode *node = meshDS->AddNodeWithID(v->x(),v->y(),v->z(),v->getNum());
        meshDS->SetNodeInVolume( node, topoSolid );
      }
    }

    //elements
    std::vector<MVertex*> verts;
    for(unsigned int i = 0; i < gRegion->getNumMeshElements(); i++)
    {
      MElement *e = gRegion->getMeshElement(i);
      verts.clear();
      e->getVertices(verts);
      SMDS_MeshVolume* volume = 0;
      switch (verts.size()){
        case 4:
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[3]->getNum(),e->getNum());
          break;
        case 5:
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[3]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[4]->getNum(),e->getNum());
          break;
        case 6:
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[3]->getNum(),
                                           verts[5]->getNum(),
                                           verts[4]->getNum(),e->getNum());
          break;
        case 8:
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[3]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[4]->getNum(),
                                           verts[7]->getNum(),
                                           verts[6]->getNum(),
                                           verts[5]->getNum(),e->getNum());
          break;
        case 10:
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[3]->getNum(),
                                           verts[6]->getNum(),
                                           verts[5]->getNum(),
                                           verts[4]->getNum(),
                                           verts[7]->getNum(),
                                           verts[8]->getNum(),
                                           verts[9]->getNum(),e->getNum());
          break;
        case 13:
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[3]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[4]->getNum(),
                                           verts[6]->getNum(),
                                           verts[10]->getNum(),
                                           verts[8]->getNum(),
                                           verts[5]->getNum(),
                                           verts[7]->getNum(),
                                           verts[12]->getNum(),
                                           verts[11]->getNum(),
                                           verts[9]->getNum(),e->getNum());
          break;
        case 14: // same as case 13, because no pyra14 in smesh
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[3]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[4]->getNum(),
                                           verts[6]->getNum(),
                                           verts[10]->getNum(),
                                           verts[8]->getNum(),
                                           verts[5]->getNum(),
                                           verts[7]->getNum(),
                                           verts[12]->getNum(),
                                           verts[11]->getNum(),
                                           verts[9]->getNum(),e->getNum());
          break;
        case 15:
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[3]->getNum(),
                                           verts[5]->getNum(),
                                           verts[4]->getNum(),
                                           verts[7]->getNum(),
                                           verts[9]->getNum(),
                                           verts[6]->getNum(),
                                           verts[13]->getNum(),
                                           verts[14]->getNum(),
                                           verts[12]->getNum(),
                                           verts[8]->getNum(),
                                           verts[11]->getNum(),
                                           verts[10]->getNum(),e->getNum());
          break;
        case 18: // same as case 15, because no penta18 in smesh
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[3]->getNum(),
                                           verts[5]->getNum(),
                                           verts[4]->getNum(),
                                           verts[7]->getNum(),
                                           verts[9]->getNum(),
                                           verts[6]->getNum(),
                                           verts[13]->getNum(),
                                           verts[14]->getNum(),
                                           verts[12]->getNum(),
                                           verts[8]->getNum(),
                                           verts[11]->getNum(),
                                           verts[10]->getNum(),e->getNum());
          break;
        case 20:
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[3]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[4]->getNum(),
                                           verts[7]->getNum(),
                                           verts[6]->getNum(),
                                           verts[5]->getNum(),
                                           verts[9]->getNum(),
                                           verts[13]->getNum(),
                                           verts[11]->getNum(),
                                           verts[8]->getNum(),
                                           verts[17]->getNum(),
                                           verts[19]->getNum(),
                                           verts[18]->getNum(),
                                           verts[16]->getNum(),
                                           verts[10]->getNum(),
                                           verts[15]->getNum(),
                                           verts[14]->getNum(),
                                           verts[12]->getNum(),e->getNum());
          break;
        case 27:
          volume = meshDS->AddVolumeWithID(verts[0]->getNum(),
                                           verts[3]->getNum(),
                                           verts[2]->getNum(),
                                           verts[1]->getNum(),
                                           verts[4]->getNum(),
                                           verts[7]->getNum(),
                                           verts[6]->getNum(),
                                           verts[5]->getNum(),
                                           verts[9]->getNum(),
                                           verts[13]->getNum(),
                                           verts[11]->getNum(),
                                           verts[8]->getNum(),
                                           verts[17]->getNum(),
                                           verts[19]->getNum(),
                                           verts[18]->getNum(),
                                           verts[16]->getNum(),
                                           verts[10]->getNum(),
                                           verts[15]->getNum(),
                                           verts[14]->getNum(),
                                           verts[12]->getNum(),
                                           verts[20]->getNum(),
                                           verts[22]->getNum(),
                                           verts[24]->getNum(),
                                           verts[23]->getNum(),
                                           verts[21]->getNum(),
                                           verts[25]->getNum(),
                                           verts[26]->getNum(),
                                           e->getNum());
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
  GmshInitialize();
  SetGmshOptions();
  _gModel = new GModel();
  mymsg msg(_gModel);
  GmshSetMessageHandler(&msg);
  _gModel->importOCCShape((void*)&_shape);
  if (_compounds.size() > 0) CreateGmshCompounds();
  MESSAGE("GModel::Mesh");
  try
  {
    _gModel->mesh((_is2d)?2:3);
#ifdef WITH_SMESH_CANCEL_COMPUTE

#endif
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
    FillSMesh();
  }
  delete _gModel;
  GmshFinalize();
  MESSAGE("GMSHPlugin_Mesher::Compute:End");
  return !err;
}
