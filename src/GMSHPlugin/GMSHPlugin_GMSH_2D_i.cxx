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
#include "GMSHPlugin_GMSH_2D_i.hxx"
#include "SMESH_Gen.hxx"

#include "Utils_CorbaException.hxx"
#include "utilities.h"

using namespace std;

//=============================================================================
/*!
 *  GMSHPlugin_GMSH_2D_i::GMSHPlugin_GMSH_2D_i
 *
 *  Constructor
 */
//=============================================================================

GMSHPlugin_GMSH_2D_i::GMSHPlugin_GMSH_2D_i( PortableServer::POA_ptr thePOA,
                                                    ::SMESH_Gen*            theGenImpl )
     : SALOME::GenericObj_i( thePOA ), 
       SMESH_Hypothesis_i( thePOA ), 
       SMESH_Algo_i( thePOA ),
       SMESH_2D_Algo_i( thePOA )
{
  MESSAGE( "GMSHPlugin_GMSH_2D_i::GMSHPlugin_GMSH_2D_i" );
  myBaseImpl = new ::GMSHPlugin_GMSH_2D( theGenImpl->GetANewId(),
                                             theGenImpl );
}

//=============================================================================
/*!
 *  GMSHPlugin_GMSH_2D_i::~GMSHPlugin_GMSH_2D_i
 *
 *  Destructor
 */
//=============================================================================

GMSHPlugin_GMSH_2D_i::~GMSHPlugin_GMSH_2D_i()
{
  MESSAGE( "GMSHPlugin_GMSH_2D_i::~GMSHPlugin_GMSH_2D_i" );
}

//=============================================================================
/*!
 *  GMSHPlugin_GMSH_2D_i::GetImpl
 *
 *  Get implementation
 */
//=============================================================================

::GMSHPlugin_GMSH_2D* GMSHPlugin_GMSH_2D_i::GetImpl()
{
  MESSAGE( "GMSHPlugin_GMSH_2D_i::GetImpl" );
  return ( ::GMSHPlugin_GMSH_2D* )myBaseImpl;
}
