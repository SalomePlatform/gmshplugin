// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016-2022  EDF R&D
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
#include "GMSHPlugin_Hypothesis_2D_i.hxx"
#include "SMESH_Gen.hxx"
#include "SMESH_PythonDump.hxx"

#include "Utils_CorbaException.hxx"
#include "utilities.h"

using namespace std;

//=============================================================================
/*!
 *  GMSHPlugin_Hypothesis_2D_i::GMSHPlugin_Hypothesis_2D_i
 *
 *  Constructor
 */
//=============================================================================
GMSHPlugin_Hypothesis_2D_i::
GMSHPlugin_Hypothesis_2D_i (PortableServer::POA_ptr thePOA,
                              ::SMESH_Gen*            theGenImpl)
  : SALOME::GenericObj_i( thePOA ),
    SMESH_Hypothesis_i( thePOA ),
    GMSHPlugin_Hypothesis_i( thePOA, theGenImpl )
{
  MESSAGE( "GMSHPlugin_Hypothesis_2D_i::GMSHPlugin_Hypothesis_2D_i" );
  if (myBaseImpl)
    delete (::GMSHPlugin_Hypothesis*)myBaseImpl;
  myBaseImpl = new ::GMSHPlugin_Hypothesis_2D (theGenImpl->GetANewId(),
                                                 theGenImpl);
}

//=============================================================================
/*!
 *  GMSHPlugin_Hypothesis_2D_i::~GMSHPlugin_Hypothesis_2D_i
 *
 *  Destructor
 */
//=============================================================================
GMSHPlugin_Hypothesis_2D_i::~GMSHPlugin_Hypothesis_2D_i()
{
  MESSAGE( "GMSHPlugin_Hypothesis_2D_i::~GMSHPlugin_Hypothesis_2D_i" );
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
CORBA::Boolean GMSHPlugin_Hypothesis_2D_i::IsDimSupported( SMESH::Dimension type )
{
  return type == SMESH::DIM_2D;
}
