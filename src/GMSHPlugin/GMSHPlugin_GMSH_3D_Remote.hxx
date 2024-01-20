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
#ifndef _GMSHPlugin_GMSH_3D_REMOTE_HXX_
#define _GMSHPlugin_GMSH_3D_REMOTE_HXX_

#include "GMSHPlugin_GMSH_3D.hxx"

#include <vector>
#include <map>

class SMDS_MeshNode;

class GMSHPLUGIN_EXPORT GMSHPlugin_GMSH_3D_Remote: public GMSHPlugin_GMSH_3D
{
 public:
  GMSHPlugin_GMSH_3D_Remote(int hypId, SMESH_Gen* gen);
  virtual ~GMSHPlugin_GMSH_3D_Remote();

  // Function whould not be used with remote Computing
  bool CheckHypothesis (SMESH_Mesh&         aMesh,
                        const TopoDS_Shape& aShape,
                        Hypothesis_Status&  aStatus) override {(void)aMesh;(void)aShape;aStatus = HYP_OK;return true;};

  bool Compute(SMESH_Mesh& aMesh, const TopoDS_Shape& aShape) override;

  void setSubMeshesToCompute(SMESH_subMesh * aSubMesh) override;


 protected:
  void exportElementOrientation(SMESH_Mesh& aMesh, const TopoDS_Shape& aShape, const std::string output_file);
  void exportGmshParams( const std::string param_file,  const SMESHDS_Hypothesis* hyp );
  typedef std::set<std::string> TCompound;
};

#endif
