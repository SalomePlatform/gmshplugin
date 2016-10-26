// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016  EDF R&D
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
#include "GMSHPlugin_GMSH.hxx"
#include "GMSHPlugin_Hypothesis.hxx"
#include "GMSHPlugin_Mesher.hxx"

#include <SMESH_Gen.hxx>
#include <SMESH_Mesh.hxx>
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

GMSHPlugin_GMSH::GMSHPlugin_GMSH(int hypId, int studyId,
                                                   SMESH_Gen* gen)
  : SMESH_3D_Algo(hypId, studyId, gen)
{
  MESSAGE("GMSHPlugin_GMSH::GMSHPlugin_GMSH");
  _name = "GMSH";
  _shapeType = (1 << TopAbs_SHELL) | (1 << TopAbs_SOLID);// 1 bit /shape type
  _compatibleHypothesis.push_back("GMSH_Parameters");
  _requireDiscreteBoundary = false;
  _onlyUnaryInput = false;
  _hypothesis = NULL;
  _supportSubmeshes = true;
}

//=============================================================================
/*!
 *  
 */
//=============================================================================

GMSHPlugin_GMSH::~GMSHPlugin_GMSH()
{
  MESSAGE("GMSHPlugin_GMSH::~GMSHPlugin_GMSH");
}

//=============================================================================
/*!
 *  
 */
//=============================================================================

bool GMSHPlugin_GMSH::CheckHypothesis
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

bool GMSHPlugin_GMSH::Compute(SMESH_Mesh&         aMesh,
                                       const TopoDS_Shape& aShape)
{
  GMSHPlugin_Mesher mesher(&aMesh, aShape);
  mesher.SetParameters(dynamic_cast<const GMSHPlugin_Hypothesis*>(_hypothesis));
  return mesher.Compute();
}

#ifdef WITH_SMESH_CANCEL_COMPUTE
void GMSHPlugin_GMSH::CancelCompute()
{}
#endif

//=============================================================================
/*!
 *
 */
//=============================================================================

bool GMSHPlugin_GMSH::Evaluate(SMESH_Mesh&         aMesh,
                                        const TopoDS_Shape& aShape,
                                        MapShapeNbElems& aResMap)
{
  std::vector<int> aResVec(SMDSEntity_Last);
  for(int i=SMDSEntity_Node; i<SMDSEntity_Last; i++) aResVec[i] = 0;
  SMESH_subMesh * sm = aMesh.GetSubMesh(aShape);
  aResMap.insert(std::make_pair(sm,aResVec));
  SMESH_ComputeErrorPtr& smError = sm->GetComputeError();
  smError.reset( new SMESH_ComputeError(COMPERR_ALGO_FAILED,"Evaluation is not implemented",this));
  
  return true;
}
