// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016-2026  EDF
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

#ifndef _GMSHPlugin_Mesher_HXX_
#define _GMSHPlugin_Mesher_HXX_

#include <stdio.h>
#include "GmshVersion.h"
#if GMSH_MAJOR_VERSION >=4
#include "gmsh.h"
#else
#include "Gmsh.h"
#endif
#include "GmshConfig.h"
#include "GModelIO_OCC.h"
#include "GModelIO_GEO.h"
#include "Geo.h"
#if GMSH_MAJOR_VERSION >=4
#include "GEdge.h"
#include "GFace.h"
#else
#include "GEdgeCompound.h"
#include "GFaceCompound.h"
#endif
#include "MElement.h"

#include <SMDS_MeshElement.hxx>
#include "GMSHPlugin_Defs.hxx"
#include "SMESH_Algo.hxx"

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <unordered_map>

class SMESH_Mesh;
class SMESH_Comment;
class SMESHDS_Mesh;
class TopoDS_Shape;
class GMSHPlugin_Hypothesis;
class StdMeshers_ViscousLayers;

//=============================================================================
/*!
 * \brief This class calls the GMSH mesher of OCC geometry
 */
//=============================================================================

class GMSHPLUGIN_EXPORT GMSHPlugin_Mesher
{
 public:
  // ---------- PUBLIC METHODS ----------

  GMSHPlugin_Mesher (SMESH_Mesh* mesh, const TopoDS_Shape& aShape, bool is2D, bool is3D);

  void SetParameters(const GMSHPlugin_Hypothesis*    hyp);
  void SetParameters(const StdMeshers_ViscousLayers* hyp );
  void SetViscousLayers2DAssigned(bool isAssigned) { _isViscousLayers2D = isAssigned; }

  bool Compute3D( std::vector< const SMDS_MeshNode* >& nodeVec, 
                  std::map<const SMDS_MeshElement*, bool, TIDCompare>& listElements, 
                  bool addElements );
  bool Compute();

  bool Evaluate(MapShapeNbElems& aResMap);

  static float DistBoundingBox(const SBoundingBox3d& bounds, const SPoint3& point);

  void FillGeomMapMeshUsing2DMeshIterator( std::map<const SMDS_MeshElement*, bool, TIDCompare>& listElements );
  GModel* GetGModel(){ return _gModel;};
  void finalizeGModel();
  const SMDS_MeshNode* Node( const MVertex* v );
  const SMDS_MeshNode* PremeshedNode( const MVertex* v );

 private:
  SMESH_Mesh*          _mesh;
  const TopoDS_Shape&  _shape;
  int                  _algo2d;
  int                  _algo3d;
  int                  _recomb2DAlgo;
  bool                 _recombineAll;
  int                  _subdivAlgo;
  int                  _remeshAlgo;
  int                  _remeshPara;
  double               _smouthSteps;
  double               _sizeFactor;
  double               _minSize, _maxSize;
  bool                 _secondOrder, _useIncomplElem;
  bool                 _is2d;
  bool                 _is3d;
  int                  _verbLvl;
  bool                 _isViscousLayers2D;
  GModel*              _gModel;
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=3
  double               _maxThreads;
#endif
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
  double               _meshCurvatureSize;
#endif

//================================================================================
// Structure for mesh geometry elements used to build a new mesh.
// - Triangle: three points - used for 3d mesh generation
//================================================================================
  struct GeomElement
  {
    gp_Pnt p1, p2, p3;

    //triangle
    GeomElement(const gp_Pnt& point1, const gp_Pnt& point2, const gp_Pnt& point3)
      : p1(point1), p2(point2), p3(point3) {}
  };

//================================================================================
// GeomPtr: smart pointer to GeomElement
// GeomVector: list of elements for a shape
// ShapeGeomMap: maps TopoDS_Shape to its elements
//================================================================================
  using GeomPtr = std::unique_ptr<GeomElement>;
  using GeomVector = std::vector<GeomPtr>;
  using ShapeGeomMap = std::unordered_map<TopoDS_Shape, GeomVector>;

  std::set<std::string> _compounds;

  std::unordered_set<TopoDS_Shape> _meshedShapes;
  std::unordered_map<int, TopoDS_Shape > _shapeGMSHMap;
  std::unordered_map<TopoDS_Shape, TopoDS_Shape > _faceGMSHMap;
  std::map< const MVertex *, const SMDS_MeshNode* > _nodeMap;
  std::map< const MVertex *, const SMDS_MeshNode* > _premeshednodeMap; // used for the SA version
  const StdMeshers_ViscousLayers* _viscousLayersHyp;
  SMESHDS_SubMesh* hasSubMesh( const TopoDS_Shape& s );

  void setGmshOptions();
  void createGmshCompounds();
  void fillSMesh(GModel* gModel);
  void fillSMesh0D(GModel* gModel, SMESHDS_Mesh* meshDS);
  void fillSMesh1D(GModel* gModel, SMESHDS_Mesh* meshDS);
  void fillSMesh2D(GModel* gModel, SMESHDS_Mesh* meshDS, bool allowFaceEdit = false);
  void fillSMesh3D(GModel* gModel, SMESHDS_Mesh* meshDS);
  void hideComputedEntities( GModel* gModel, bool hideAnyway = false );
  void restoreVisibility( GModel* gModel );
  void set1DSubMeshes( GModel* );
  void set2DSubMeshes( GModel* );
  void toPython( GModel* );
  bool isAllNodesInSameFace( const SMDS_MeshElement* triangle, const TopoDS_Face& F, std::vector<gp_XY>& uvValues );
  std::map<int,std::vector<std::tuple<smIdType,bool,std::vector<gp_XY>>>> associateElementsToFaces( std::map<const SMDS_MeshElement*, bool, TIDCompare>& listElements );
  void set2DMeshes(std::vector< const SMDS_MeshNode* >& nodeVec,
                   std::map<const SMDS_MeshElement*, bool, TIDCompare>& listElements );

//================================================================================
  void adapt2DQuadMeshAndComputeGmsh3D();

  /**
   * @brief Generates a GMSH mesh.
   * Depending on the mesh dimension, either `theShape` or `theShapeGeomMap` is used:
   * - For 2D meshes (e.g. surface meshing), `theShape` must be provided and is used as the input geometry.
   * - For 3D meshes (e.g. volume meshing), `theShapeGeomMap` must be provided and maps each shape to its geometric definition.
   *
   * @param theUniqBoundNodes   Map of unique boundary nodes used for gmsh mesh transition to smesh mesh.
   * @param theMeshedShapes     Auxiliary vector to collect resulting TopoDS_Shapes that were meshed.
   * @param meshDimension       Mesh dimension (2 or 3).
   * @param theShapeGeomMap     Optional geometry map used for 3D meshing.
   * @param theShape            Optional shape used for 2D meshing.
   */
  void genGMSHMesh(std::unordered_map<gp_Pnt, const SMDS_MeshNode*> & theUniqBoundNodes,
                   std::vector<TopoDS_Shape>& theMeshedShapes,
                   int meshDimension,
                   ShapeGeomMap* theShapeGeomMap = nullptr,
                   const TopoDS_Shape& theShape = TopoDS_Shape());

  void prepareTriangles(const std::map<const SMDS_MeshElement*, bool, TIDCompare>& theListElements,
                        std::unordered_map<gp_Pnt, const SMDS_MeshNode*> & theUniqBoundNodes,
                        SMESH_MesherHelper& theMeshHelper,
                        GeomVector& theTriangles);

  bool fill2DViscousLayer();
  bool fill3DViscousLayerAndQuadAdapt();
//================================================================================

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=3
  void SetMaxThreadsGmsh();
  void SetCompoundMeshVisibility();
#endif
  class mymsg : public GmshMessage
  {
  private:
    GModel* _gModel;
  public:
    mymsg(GModel* _gModel) :  GmshMessage(), _gModel(_gModel) {}
    void operator()(std::string level, std::string msg);
  };
};

#endif
