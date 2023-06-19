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
#include "GMSHPlugin_Hypothesis_i.hxx"
#include "SMESH_Gen.hxx"
#include "SMESH_PythonDump.hxx"
#include "GEOM_Object.hxx"

#include "Utils_CorbaException.hxx"
#include "utilities.h"

using namespace std;

template<>
bool GMSHPlugin_Hypothesis_i::isToSetParameter<double>(double curValue,
                                                         double newValue,
                                                         /*SettingMethod*/int meth)
{
  return isToSetParameter(true, (fabs(curValue - newValue) < 1e-20), meth);
}

GMSHPlugin_Hypothesis_i::
GMSHPlugin_Hypothesis_i (PortableServer::POA_ptr thePOA,
                           ::SMESH_Gen*            theGenImpl)
  : SALOME::GenericObj_i( thePOA ), 
    SMESH_Hypothesis_i( thePOA ),
    mySetMethodFlags(0)
{
  MESSAGE( "GMSHPlugin_Hypothesis_i::GMSHPlugin_Hypothesis_i" );
  myBaseImpl = new ::GMSHPlugin_Hypothesis (theGenImpl->GetANewId(),
                                              theGenImpl);
}

GMSHPlugin_Hypothesis_i::~GMSHPlugin_Hypothesis_i()
{
  MESSAGE( "GMSHPlugin_Hypothesis_i::~GMSHPlugin_Hypothesis_i" );
}

void GMSHPlugin_Hypothesis_i::SetMaxSize (CORBA::Double theValue)
{
  if ( isToSetParameter( GetMaxSize(), theValue, METH_SetMaxSize ))
  {
    this->GetImpl()->SetMaxSize(theValue);
    SMESH::TPythonDump() << _this() << ".SetMaxSize( " << SMESH::TVar(theValue) << " )";
  }
}

CORBA::Double GMSHPlugin_Hypothesis_i::GetMaxSize()
{
  return this->GetImpl()->GetMaxSize();
}


void GMSHPlugin_Hypothesis_i::SetMeshCurvatureSize (CORBA::Double theMeshCurvatureSize)
{
  if ( isToSetParameter( GetMeshCurvatureSize(), theMeshCurvatureSize, METH_SetMeshCurvatureSize ))
  {
    this->GetImpl()->SetMeshCurvatureSize(theMeshCurvatureSize);
    SMESH::TPythonDump() << _this() << ".SetMeshCurvatureSize( " << SMESH::TVar(theMeshCurvatureSize) << " )";
  }
}

CORBA::Double GMSHPlugin_Hypothesis_i::GetMeshCurvatureSize()
{
  return this->GetImpl()->GetMeshCurvatureSize();
}


void GMSHPlugin_Hypothesis_i::SetMinSize (CORBA::Double theValue)
{
  if ( isToSetParameter( GetMinSize(), theValue, METH_SetMinSize ))
  {
    this->GetImpl()->SetMinSize(theValue);
    SMESH::TPythonDump() << _this() << ".SetMinSize( " << SMESH::TVar(theValue) << " )";
  }
}

CORBA::Double GMSHPlugin_Hypothesis_i::GetMinSize()
{
  return this->GetImpl()->GetMinSize();
}

void GMSHPlugin_Hypothesis_i::SetSecondOrder (CORBA::Boolean theValue)
{
  if ( isToSetParameter( GetSecondOrder(), theValue, METH_SetSecondOrder ))
  {
    this->GetImpl()->SetSecondOrder(theValue);
    SMESH::TPythonDump() << _this() << ".SetSecondOrder( " << theValue << " )";
  }
}

CORBA::Boolean GMSHPlugin_Hypothesis_i::GetSecondOrder()
{
  return this->GetImpl()->GetSecondOrder();
}

void GMSHPlugin_Hypothesis_i::Set2DAlgo (CORBA::Long theValue)
{
  if ( isToSetParameter( Get2DAlgo(), theValue, METH_Set2DAlgo ))
  {
    this->GetImpl()->Set2DAlgo((::GMSHPlugin_Hypothesis::Algo2D)theValue);
    SMESH::TPythonDump() << _this() << ".Set2DAlgo( " << theValue << " )";
  }
}

CORBA::Long GMSHPlugin_Hypothesis_i::Get2DAlgo()
{
  return this->GetImpl()->Get2DAlgo();
}

void GMSHPlugin_Hypothesis_i::Set3DAlgo (CORBA::Long theValue)
{
  if ( isToSetParameter( Get3DAlgo(), theValue, METH_Set3DAlgo ))
  {
    this->GetImpl()->Set3DAlgo((::GMSHPlugin_Hypothesis::Algo3D)theValue);
    SMESH::TPythonDump() << _this() << ".Set3DAlgo( " << theValue << " )";
  }
}

CORBA::Long GMSHPlugin_Hypothesis_i::Get3DAlgo()
{
  return this->GetImpl()->Get3DAlgo();
}

void GMSHPlugin_Hypothesis_i::SetRecomb2DAlgo (CORBA::Long theValue)
{
  if ( isToSetParameter( GetRecomb2DAlgo(), theValue, METH_SetRecomb2DAlgo ))
  {
    this->GetImpl()->SetRecomb2DAlgo((::GMSHPlugin_Hypothesis::Recomb2DAlgo)theValue);
    SMESH::TPythonDump() << _this() << ".SetRecomb2DAlgo( " << theValue << " )";
  }
}

CORBA::Long GMSHPlugin_Hypothesis_i::GetRecomb2DAlgo()
{
  return this->GetImpl()->GetRecomb2DAlgo();
}

void GMSHPlugin_Hypothesis_i::SetRecombineAll (CORBA::Boolean theValue)
{
  if ( isToSetParameter( GetRecombineAll(), theValue, METH_SetRecombineAll ))
  {
    this->GetImpl()->SetRecombineAll(theValue);
    SMESH::TPythonDump() << _this() << ".SetRecombineAll( " << theValue << " )";
  }
}

CORBA::Boolean GMSHPlugin_Hypothesis_i::GetRecombineAll()
{
  return this->GetImpl()->GetRecombineAll();
}

void GMSHPlugin_Hypothesis_i::SetSubdivAlgo (CORBA::Long theValue)
{
  if ( isToSetParameter( GetSubdivAlgo(), theValue, METH_SetSubdivAlgo ))
  {
    this->GetImpl()->SetSubdivAlgo((::GMSHPlugin_Hypothesis::SubdivAlgo)theValue);
    SMESH::TPythonDump() << _this() << ".SetSubdivAlgo( " << theValue << " )";
  }
}

CORBA::Long GMSHPlugin_Hypothesis_i::GetSubdivAlgo()
{
  return this->GetImpl()->GetSubdivAlgo();
}

void GMSHPlugin_Hypothesis_i::SetRemeshAlgo (CORBA::Long theValue)
{
  if ( isToSetParameter( GetRemeshAlgo(), theValue, METH_SetRemeshAlgo ))
  {
    this->GetImpl()->SetRemeshAlgo((::GMSHPlugin_Hypothesis::RemeshAlgo)theValue);
    SMESH::TPythonDump() << _this() << ".SetRemeshAlgo( " << theValue << " )";
  }
}

CORBA::Long GMSHPlugin_Hypothesis_i::GetRemeshAlgo()
{
  return this->GetImpl()->GetRemeshAlgo();
}

void GMSHPlugin_Hypothesis_i::SetRemeshPara (CORBA::Long theValue)
{
  if ( isToSetParameter( GetRemeshPara(), theValue, METH_SetRemeshPara ))
  {
    this->GetImpl()->SetRemeshPara((::GMSHPlugin_Hypothesis::RemeshPara)theValue);
    SMESH::TPythonDump() << _this() << ".SetRemeshAlgo( " << theValue << " )";
  }
}

CORBA::Long GMSHPlugin_Hypothesis_i::GetRemeshPara()
{
  return this->GetImpl()->GetRemeshPara();
}

void GMSHPlugin_Hypothesis_i::SetSmouthSteps (CORBA::Double theValue)
{
  if ( isToSetParameter( GetSmouthSteps(), theValue, METH_SetSmouthSteps ))
  {
    this->GetImpl()->SetSmouthSteps(theValue);
    SMESH::TPythonDump() << _this() << ".SetSmouthSteps( " << theValue << " )";
  }
}

CORBA::Double GMSHPlugin_Hypothesis_i::GetSmouthSteps()
{
  return this->GetImpl()->GetSmouthSteps();
}

void GMSHPlugin_Hypothesis_i::SetSizeFactor (CORBA::Double theValue)
{
  if ( isToSetParameter( GetSizeFactor(), theValue, METH_SetSizeFactor ))
  {
    this->GetImpl()->SetSizeFactor(theValue);
    SMESH::TPythonDump() << _this() << ".SetSizeFactor( " << theValue << " )";
  }
}

CORBA::Double GMSHPlugin_Hypothesis_i::GetSizeFactor()
{
  return this->GetImpl()->GetSizeFactor();
}

void GMSHPlugin_Hypothesis_i::SetUseIncomplElem (CORBA::Boolean theValue)
{
  if ( isToSetParameter( GetUseIncomplElem(), theValue, METH_SetUseIncomplElem ))
  {
    this->GetImpl()->SetUseIncomplElem(theValue);
    SMESH::TPythonDump() << _this() << ".SetUseIncomplElem( " << theValue << " )";
  }
}

CORBA::Boolean GMSHPlugin_Hypothesis_i::GetUseIncomplElem()
{
  return this->GetImpl()->GetUseIncomplElem();
}

void GMSHPlugin_Hypothesis_i::SetIs2d (CORBA::Boolean theValue)
{
  this->GetImpl()->SetIs2d(theValue);
  SMESH::TPythonDump() << _this() << ".SetIs2d( " << theValue << " )";
}

void GMSHPlugin_Hypothesis_i::SetCompoundOnShape(GEOM::GEOM_Object_ptr GeomObj)
{
  string entry;
  entry = GeomObj->GetStudyEntry();
  SetCompoundOnEntry(entry.c_str());
}

void GMSHPlugin_Hypothesis_i::SetCompoundOnEntry(const char*   entry)
{
  //if ( isToSetParameter( GetCompoundOnEntry(entry), METH_SetCompoundOnEntry ))
  {
    this->GetImpl()->SetCompoundOnEntry(entry);
    SMESH::TPythonDump()
      << _this() << ".SetCompoundOnShape(" << entry << ")";
  }
}

GMSHPlugin::string_array* GMSHPlugin_Hypothesis_i::GetCompoundOnEntries()
{
  GMSHPlugin::string_array_var result = new GMSHPlugin::string_array();
  const ::GMSHPlugin_Hypothesis::TCompound compounds =
    this->GetImpl()->GetCompoundOnEntries();
  result->length(compounds.size());
  ::GMSHPlugin_Hypothesis::TCompound::const_iterator it = compounds.begin();
  for (int i=0 ; it != compounds.end() ; i++, it++)
  {
    string entry = *it;
    result[i] = CORBA::string_dup(entry.c_str());
  }
  return result._retn();
}

void GMSHPlugin_Hypothesis_i::UnsetCompoundOnShape(GEOM::GEOM_Object_ptr GeomObj)
{
  string entry;
  entry = GeomObj->GetStudyEntry();
  UnsetCompoundOnEntry(entry.c_str());
}

void GMSHPlugin_Hypothesis_i::UnsetCompoundOnEntry(const char* entry)
{
  this->GetImpl()->UnsetCompoundOnEntry(entry);
  SMESH::TPythonDump() << _this() << ".UnsetCompoundOnShape(" << entry << ")";
}

::GMSHPlugin_Hypothesis* GMSHPlugin_Hypothesis_i::GetImpl()
{
  return (::GMSHPlugin_Hypothesis*)myBaseImpl;
}

//================================================================================
/*!
 * \brief Verify whether hypothesis supports given entity type 
  * \param type - dimension (see SMESH::Dimension enumeration)
  * \retval CORBA::Boolean - TRUE if dimension is supported, FALSE otherwise
 * 
 * Verify whether hypothesis supports given entity type (see SMESH::Dimension enumeration)
 */
//================================================================================  
CORBA::Boolean GMSHPlugin_Hypothesis_i::IsDimSupported( SMESH::Dimension type )
{
  return type == SMESH::DIM_3D;
}

//================================================================================
/*!
 * \brief method intended to remove explicit treatment of Netagen hypotheses from SMESH_NoteBook
 */
//================================================================================

int GMSHPlugin_Hypothesis_i::getParamIndex(const TCollection_AsciiString& method,
                                             int nbVars) const
{
  if ( method == "SetMaxSize"        ) return 0;
  if ( method == "SetGrowthRate"     ) return 1;
  if ( method == "SetNbSegPerEdge"   ) return 2;
  if ( method == "SetNbSegPerRadius" ) return 3;
  if ( method == "SetMinSize" )        return nbVars-1;
  if ( method == "SetMeshCurvatureSize" )        return 5;

  return SMESH_Hypothesis_i::getParamIndex( method, nbVars ); // return default value
}

//================================================================================
/*!
 * \brief Method used to convert variable parameters stored in an old study
 * into myMethod2VarParams. It should return a method name for an index of
 * variable parameters. Index is countered from zero
 */
//================================================================================

std::string GMSHPlugin_Hypothesis_i::getMethodOfParameter(const int paramIndex,
                                                            int nbVars) const
{
  switch ( paramIndex ) {
  case 0: return "SetMaxSize";
  case 1: return nbVars == 2 ? "SetMinSize" : "SetGrowthRate";
  case 2: return "SetNbSegPerEdge";
  case 3: return "SetNbSegPerRadius";
  case 4: return "SetMinSize";
  case 5: return "SetMeshCurvatureSize";
  }
  return "";
}

//================================================================================
/*!
 * \brief Return geometry this hypothesis depends on. Return false if there is no geometry parameter
 */
//================================================================================

bool
GMSHPlugin_Hypothesis_i::getObjectsDependOn( std::vector< std::string > & entryArray,
                                             std::vector< int >         & /*subIDArray*/ ) const
{
  typedef ::GMSHPlugin_Hypothesis THyp;
  const THyp* impl = static_cast<const THyp*>( myBaseImpl );

  const THyp::TCompound& compounds = impl->GetCompoundOnEntries();
  entryArray.assign( compounds.cbegin(), compounds.cend() );

  return true;
}

//================================================================================
/*!
 * \brief Set new geometry instead of that returned by getObjectsDependOn()
 */
//================================================================================

bool
GMSHPlugin_Hypothesis_i::setObjectsDependOn( std::vector< std::string > & entryArray,
                                             std::vector< int >         & /*subIDArray*/ )
{
  typedef ::GMSHPlugin_Hypothesis THyp;
  THyp* impl = static_cast< THyp* >( myBaseImpl );

  size_t iEnt = 0;

  THyp::TCompound& compoundsNew = const_cast< THyp::TCompound& > ( impl->GetCompoundOnEntries() );
  THyp::TCompound compounds;
  compounds.swap( compoundsNew );

  THyp::TCompound::const_iterator entry = compounds.cbegin();
  for ( ; entry != compounds.cend(); ++entry, ++iEnt )
    if ( !entryArray[ iEnt ].empty() )
      compoundsNew.insert( entryArray[ iEnt ]);

  return true;
}
