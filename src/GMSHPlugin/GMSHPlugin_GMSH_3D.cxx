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
#include "GMSHPlugin_GMSH_3D.hxx"
#include "GMSHPlugin_Hypothesis_2D.hxx"

#include <SMESH_Gen.hxx>
#include <SMESH_ControlsDef.hxx>
#include <SMESHDS_Mesh.hxx>
#include <utilities.h>

#include <list>

using namespace std;

//=============================================================================
/*!
 *  
 */
//=============================================================================

GMSHPlugin_GMSH_3D::GMSHPlugin_GMSH_3D(int hypId, SMESH_Gen* gen)
  : SMESH_3D_Algo(hypId, gen)
{
  MESSAGE("GMSHPlugin_GMSH_3D::GMSHPlugin_GMSH_3D");
  _name = "GMSH_3D";
  _shapeType = (1 << TopAbs_SOLID); // 1 bit /shape type
  _compatibleHypothesis.push_back("GMSH_Parameters_3D");
  _onlyUnaryInput = false;
  _hypothesis = NULL;
  _supportSubmeshes = true;

  _requireShape = false; // can work without shape
}

//=============================================================================
/*!
 *  
 */
//=============================================================================

GMSHPlugin_GMSH_3D::~GMSHPlugin_GMSH_3D()
{
  MESSAGE("GMSHPlugin_GMSH_3D::~GMSHPlugin_GMSH_3D");
}

//=============================================================================
/*!
 *  
 */
//=============================================================================

bool GMSHPlugin_GMSH_3D::CheckHypothesis
                         (SMESH_Mesh& aMesh,
                          const TopoDS_Shape& aShape,
                          SMESH_Hypothesis::Hypothesis_Status& aStatus)
{
  MESSAGE("GMSHPlugin_GMSH::CheckHypothesis");
  _hypothesis = NULL;
  
  const list<const SMESHDS_Hypothesis*>& hyps = GetUsedHypothesis(aMesh, aShape);
  int nbHyp = hyps.size();
  if (!nbHyp)
  {
    aStatus = SMESH_Hypothesis::HYP_OK;
    return true;  // can work with no hypothesis
  }
  // use only the first hypothesis
  const SMESHDS_Hypothesis* theHyp = hyps.front();
  
  string hypName = theHyp->GetName();
  if ( find( _compatibleHypothesis.begin(), _compatibleHypothesis.end(),
             hypName ) != _compatibleHypothesis.end() )
  {
    _hypothesis = theHyp;
    aStatus = SMESH_Hypothesis::HYP_OK;
  }
  else
  {
    aStatus = SMESH_Hypothesis::HYP_INCOMPATIBLE;
  }

  return aStatus == SMESH_Hypothesis::HYP_OK;
}

//=============================================================================
/*!
 *Here we are going to use the GMSH mesher
 */
//=============================================================================

bool GMSHPlugin_GMSH_3D::Compute(SMESH_Mesh&         aMesh,
                                 const TopoDS_Shape& aShape)
{
  GMSHPlugin_Mesher mesher(&aMesh, aShape,/*2d=*/false, true);
  mesher.SetParameters(dynamic_cast<const GMSHPlugin_Hypothesis*>(_hypothesis));

  return mesher.Compute();
}

#ifdef WITH_SMESH_CANCEL_COMPUTE
void GMSHPlugin_GMSH_3D::CancelCompute()
{}
#endif

//=============================================================================
/*!
 *
 */
//=============================================================================

bool GMSHPlugin_GMSH_3D::Evaluate(SMESH_Mesh&         aMesh,
                                  const TopoDS_Shape& aShape,
                                  MapShapeNbElems& aResMap)
{
  std::vector<smIdType> aResVec(SMDSEntity_Last);
  for(smIdType i=SMDSEntity_Node; i<SMDSEntity_Last; i++) aResVec[i] = 0;
  SMESH_subMesh * sm = aMesh.GetSubMesh(aShape);
  aResMap.insert(std::make_pair(sm,aResVec));
  SMESH_ComputeErrorPtr& smError = sm->GetComputeError();
  smError.reset( new SMESH_ComputeError(COMPERR_ALGO_FAILED,"Evaluation is not implemented",this));
  
  return true;
}
