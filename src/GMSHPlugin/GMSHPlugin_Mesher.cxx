// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016-2019  EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

using namespace std;

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
  }
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
  mapAlgo2d[0]=2; mapAlgo2d[1]=1; mapAlgo2d[2]=5; mapAlgo2d[3]=6; mapAlgo2d[4]=8; mapAlgo2d[5]=9;
  std::map <int,double> mapAlgo3d;
  mapAlgo3d[0]=1; mapAlgo3d[1]=4; mapAlgo3d[2]=5; mapAlgo3d[3]=6; mapAlgo3d[4]=7; mapAlgo3d[5]=9;

  int ok;
  ok = GmshSetOption("Mesh", "Algorithm"                , mapAlgo2d[_algo2d])    ;
  ASSERT(ok);
  if ( !_is2d)
    {
    ok = GmshSetOption("Mesh", "Algorithm3D"            , mapAlgo2d[_algo3d])    ;
    ASSERT(ok);
    }
  ok = GmshSetOption("Mesh", "RecombinationAlgorithm"   , (double)_recomb2DAlgo) ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "RecombineAll"             , (_recombineAll)?1.:0.) ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "SubdivisionAlgorithm"     , (double)_subdivAlgo)   ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "RemeshAlgorithm"          , (double)_remeshAlgo)   ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "RemeshParametrization"    , (double)_remeshPara)   ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "Smoothing"                , (double)_smouthSteps)  ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "CharacteristicLengthFactor", _sizeFactor)          ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "CharacteristicLengthMin"   , _minSize)        ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "CharacteristicLengthMax"   , _maxSize)        ;
  ASSERT(ok);
  ok = GmshSetOption("Mesh", "ElementOrder"             , (_secondOrder)?2.:1.)  ;
  ASSERT(ok);
  if (_secondOrder)
    {
    ok = GmshSetOption("Mesh", "SecondOrderIncomplete"  ,(_useIncomplElem)?1.:0.);
    ASSERT(ok);
    }
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
    SALOMEDS::SObject_var aSObj = SMESH_Gen_i::getStudyServant()->FindObjectID( (*its).c_str() );
    SALOMEDS::GenericAttribute_var anAttr;
    if (!aSObj->_is_nil() && aSObj->FindAttribute(anAttr, "AttributeIOR"))
    {
      SALOMEDS::AttributeIOR_var anIOR = SALOMEDS::AttributeIOR::_narrow(anAttr);
      CORBA::String_var aVal = anIOR->Value();
      CORBA::Object_var obj = SMESH_Gen_i::getStudyServant()->ConvertIORToObject(aVal);
      aGeomObj = GEOM::GEOM_Object::_narrow(obj);
    }
    if ( !aGeomObj->_is_nil() )
      geomShape = smeshGen_i->GeomObjectToShape( aGeomObj.in() );

    TopAbs_ShapeEnum geomType = geomShape.ShapeType();
    if ( geomType == TopAbs_COMPOUND)// voir s'il ne faut pas mettre une erreur dans le cas contraire
    {
      MESSAGE("shapeType == TopAbs_COMPOUND");
      TopoDS_Iterator it(geomShape);
      TopAbs_ShapeEnum shapeType = it.Value().ShapeType();
#if GMSH_MAJOR_VERSION >=4
      std::vector< std::pair< int, int > > dimTags;
      for ( ; it.More(); it.Next())
      {
        const TopoDS_Shape& topoShape = it.Value();
        ASSERT(topoShape.ShapeType() == shapeType);
        occgeo->importShapes( &topoShape, false, dimTags );
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
#else
      // compound of edges
      if (shapeType == TopAbs_EDGE)
      {
        MESSAGE("    shapeType == TopAbs_EDGE :");
        int num = _gModel->getNumEdges()+1;
        Curve *curve = CreateCurve(num, MSH_SEGM_COMPOUND, 1, NULL, NULL, -1, -1, 0., 1.);
        for ( ; it.More(); it.Next())
        {
          TopoDS_Shape topoShape = it.Value();
          ASSERT(topoShape.ShapeType() == shapeType);
          curve->compound.push_back(occgeo->addEdgeToModel(_gModel, (TopoDS_Edge&)topoShape)->tag());
        }
        toSynchronize = true;
        Tree_Add(_gModel->getGEOInternals()->Curves, &curve);
        //_gModel->importGEOInternals();
      }
      // compound of faces
      else if (shapeType == TopAbs_FACE)
      {
        MESSAGE("    shapeType == TopAbs_FACE :");
        int num = _gModel->getNumFaces()+1;
        Surface *surface = CreateSurface(num, MSH_SURF_COMPOUND);
        for ( ; it.More(); it.Next())
        {
          TopoDS_Shape topoShape = it.Value();
          ASSERT(topoShape.ShapeType() == shapeType);
          surface->compound.push_back(occgeo->addFaceToModel(_gModel, (TopoDS_Face&)topoShape)->tag());
        }
        toSynchronize = true;
        Tree_Add(_gModel->getGEOInternals()->Surfaces, &surface);
      }
#endif
      if ( toSynchronize )
        _gModel->getGEOInternals()->synchronize(_gModel);
    }
  }
}

//================================================================================
/*!
 * \brief Write mesh from GModel instance to SMESH instance
 */
//================================================================================

void GMSHPlugin_Mesher::FillSMesh()
{
  MESSAGE("GMSHPlugin_Mesher::FillSMesh");
  // GET MESH FROM SMESH
  SMESHDS_Mesh* meshDS = _mesh->GetMeshDS();
  
  // ADD 0D ELEMENTS
  for(GModel::viter it = _gModel->firstVertex(); it != _gModel->lastVertex(); ++it){
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
    //elements
    for(unsigned int i = 0; i < gVertex->getNumMeshElements(); i++)
    {
      MElement *e = gVertex->getMeshElement(i);
      std::vector<MVertex*> verts;
      e->getVertices(verts);
      ASSERT(verts.size()==1);
      //SMDS_Mesh0DElement* zeroDElement = 0;
      // WE DON'T ADD 0D ELEMENTS because it does not follow the salome meshers philosophy
      //zeroDElement = meshDS->Add0DElementWithID(verts[0]->getNum(),e->getNum());
      //meshDS->SetMeshElementOnShape(zeroDElement, topoVertex);
    }
  }
  
  // ADD 1D ELEMENTS
  for(GModel::eiter it = _gModel->firstEdge(); it != _gModel->lastEdge(); ++it)
  {
    GEdge *gEdge = *it;
    
    // GET topoEdge CORRESPONDING TO gEdge (map if compound)
    TopoDS_Edge topoEdge;
    std::map<GEdge*,TopoDS_Edge> topoEdges;
    
    if(gEdge->geomType() != GEntity::CompoundCurve)
    {
      topoEdge = *((TopoDS_Edge*)gEdge->getNativePtr());
      if (gEdge->getVisibility() == 0) // belongs to a compound
      {
        SMESH_subMesh* sm = _mesh->GetSubMesh(topoEdge);
        sm->SetIsAlwaysComputed(true); // prevent from displaying errors
        continue;
      }
    }
#if GMSH_MAJOR_VERSION >=4
    for ( size_t i = 0; i < gEdge->_compound.size(); ++i )
    {
      GEdge* gE = static_cast< GEdge* >( gEdge->_compound[ i ]);
      topoEdges.insert( std::make_pair( gE, *((TopoDS_Edge*)gE->getNativePtr())) );
    }
#else
    else 
    {
      // compound case, map construction GEdge/TopoDS_Edge
      std::vector<GEdge*> gEdges = ((GEdgeCompound*)gEdge)->getCompounds();
      for(std::vector<GEdge*>::const_iterator itv = gEdges.begin();itv != gEdges.end(); ++itv)
      {
        topoEdges.insert( pair<GEdge*,TopoDS_Edge>(*itv,*((TopoDS_Edge*)(*itv)->getNativePtr())) );
      }
    }
#endif

    // FILL SMESH FOR topoEdge
    //nodes
    for(unsigned int i = 0; i < gEdge->mesh_vertices.size(); i++)
    {
      MVertex *v = gEdge->mesh_vertices[i];
      if(v->getIndex() >= 0)
      {
        SMDS_MeshNode *node = meshDS->AddNodeWithID(v->x(),v->y(),v->z(),v->getNum());
        
        if (topoEdges.size() != 0)
        {
          // get back corresponding topoEdge in compound case
          topoEdge.Nullify();
          SPoint3 point = v->point();
          float distmin = std::numeric_limits<float>::max();
          for(map<GEdge*,TopoDS_Edge>::const_iterator itm=topoEdges.begin();itm != topoEdges.end(); itm++)
          { 
            SBoundingBox3d bounds = (*itm).first->bounds();
            float dist = DistBoundingBox(bounds,point);
            if (dist < distmin)
            {
              topoEdge = (*itm).second;
              distmin = dist;
              if (distmin == 0.)
                break;
            }
          }
        }
        
        meshDS->SetNodeOnEdge( node, topoEdge );
      }
    }
    //elements
    for(unsigned int i = 0; i < gEdge->getNumMeshElements(); i++)
    {
      MElement *e = gEdge->getMeshElement(i);
      std::vector<MVertex*> verts;
      e->getVertices(verts);
      SMDS_MeshEdge* edge = 0;
      
      if (topoEdges.size() !=0)
      {
        // get back corresponding topoEdge in compound case
        topoEdge.Nullify();
        SPoint3 point = e->barycenter();
        float distmin = std::numeric_limits<float>::max();
        for(map<GEdge*,TopoDS_Edge>::const_iterator itm=topoEdges.begin();itm != topoEdges.end(); itm++)
        {
          SBoundingBox3d bounds = (*itm).first->bounds();
          float dist = DistBoundingBox(bounds,point);
          if (dist < distmin)
          {
            topoEdge = (*itm).second;
            distmin = dist;
            if (distmin == 0.)
              break;
          }
        }
      }
      
      // if a node wasn't set, it is assigned here
      for (unsigned j = 0; j < verts.size(); j++)
      {
        if(verts[j]->onWhat()->getVisibility() == 0)
        {
          SMDS_MeshNode *node = meshDS->AddNodeWithID(verts[i]->x(),verts[j]->y(),verts[j]->z(),verts[j]->getNum());
          meshDS->SetNodeOnEdge( node, topoEdge );
          verts[j]->setEntity(gEdge);
        }
      }
      
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
      meshDS->SetMeshElementOnShape(edge, topoEdge);
    }
  }

  // ADD 2D ELEMENTS
  for(GModel::fiter it = _gModel->firstFace(); it != _gModel->lastFace(); ++it){
    GFace *gFace = *it;
    
    // GET topoFace CORRESPONDING TO gFace (map if compound)
    TopoDS_Face topoFace;
    std::map<GFace*,TopoDS_Face> topoFaces;
    
    if(gFace->geomType() != GEntity::CompoundSurface)
    {
      topoFace = *((TopoDS_Face*)gFace->getNativePtr());
      if (gFace->getVisibility() == 0) // belongs to a compound
      {
        SMESH_subMesh* sm = _mesh->GetSubMesh(topoFace);
        sm->SetIsAlwaysComputed(true); // prevent from displaying errors
        continue;
      }
    }
#if GMSH_MAJOR_VERSION >=4
    for ( size_t i = 0; i < gFace->_compound.size(); ++i )
    {
      GFace* gF = static_cast< GFace* >( gFace->_compound[ i ]);
      topoFaces.insert( std::make_pair( gF, *((TopoDS_Face*)gF->getNativePtr())) );
    }
#else
    else 
    {
      // compound case, map construction GFace/TopoDS_Face
      std::list<GFace*> gFaces = ((GFaceCompound*)gFace)->getCompounds();
      for(std::list<GFace*>::const_iterator itl = gFaces.begin();itl != gFaces.end(); ++itl)
      {
        topoFaces.insert( std::make_pair( *itl, *((TopoDS_Face*)(*itl)->getNativePtr())) );
      }
    }
#endif
    
    // FILL SMESH FOR topoFace
    //nodes
    for(unsigned int i = 0; i < gFace->mesh_vertices.size(); i++)
    {
      MVertex *v = gFace->mesh_vertices[i];
      if(v->getIndex() >= 0)
      {
        SMDS_MeshNode *node = meshDS->AddNodeWithID(v->x(),v->y(),v->z(),v->getNum()); 
        
        if (topoFaces.size() != 0)
        {
          // get back corresponding topoFace in compound case
          topoFace.Nullify();
          SPoint3 point = v->point();
          float distmin = std::numeric_limits<float>::max();
          for(map<GFace*,TopoDS_Face>::const_iterator itm=topoFaces.begin();itm != topoFaces.end(); itm++)
          { 
            SBoundingBox3d bounds = (*itm).first->bounds();
            float dist = DistBoundingBox(bounds,point);
            if (dist < distmin)
            {
              topoFace = (*itm).second;
              distmin = dist;
              if (distmin == 0.)
                break;
            }
          }
        }
        
        meshDS->SetNodeOnFace( node, topoFace );
      }
    }
    //elements
    for(unsigned int i = 0; i < gFace->getNumMeshElements(); i++)
    {
      MElement *e = gFace->getMeshElement(i);
      std::vector<MVertex*> verts;
      e->getVertices(verts);
      SMDS_MeshFace* face = 0;
      
      if (topoFaces.size() !=0)
      {
        // get back corresponding topoFace in compound case
        topoFace.Nullify();
        SPoint3 point = e->barycenter();
        float distmin = std::numeric_limits<float>::max();
        for(map<GFace*,TopoDS_Face>::const_iterator itm=topoFaces.begin();itm != topoFaces.end(); itm++)
        {
          SBoundingBox3d bounds = (*itm).first->bounds();
          float dist = DistBoundingBox(bounds,point);
          if (dist < distmin)
          {
            topoFace = (*itm).second;
            distmin = dist;
            if (distmin == 0.)
              break;
          }
        }
      }
      
      // if a node wasn't set, it is assigned here
      for (unsigned j = 0; j < verts.size(); j++)
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
      meshDS->SetMeshElementOnShape(face, topoFace);
    }
  }
  
  // ADD 3D ELEMENTS
  for(GModel::riter it = _gModel->firstRegion(); it != _gModel->lastRegion(); ++it){
    GRegion *gRegion = *it;
    if (gRegion->getVisibility() == 0)
      continue;
    
    // GET topoSolid CORRESPONDING TO gRegion
    TopoDS_Solid topoSolid = *((TopoDS_Solid*)gRegion->getNativePtr());
    
    //printf("volume %d (%d) contains %d elements\n", meshDS->ShapeToIndex(topoSolid), gRegion->tag(), gRegion->getNumMeshElements());
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
    for(unsigned int i = 0; i < gRegion->getNumMeshElements(); i++)
    {
      MElement *e = gRegion->getMeshElement(i);
      std::vector<MVertex*> verts;
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

float GMSHPlugin_Mesher::DistBoundingBox(SBoundingBox3d& bounds,SPoint3& point)
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
  
  return sqrt(x*x+y*y+z*z);
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
        printf(oss.str().c_str());
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
#if GMSH_MAJOR_VERSION < 4
    if (_compounds.size() > 0) _gModel->setCompoundVisibility();
#endif
    FillSMesh();
  }
  delete _gModel;
  GmshFinalize();
  MESSAGE("GMSHPlugin_Mesher::Compute:End");
  return !err;
}
