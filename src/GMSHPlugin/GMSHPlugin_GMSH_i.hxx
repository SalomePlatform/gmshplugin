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
#ifndef _GMSHPlugin_GMSH_I_HXX_
#define _GMSHPlugin_GMSH_I_HXX_

#include "GMSHPlugin_Defs.hxx"

#include <SALOMEconfig.h>
#include CORBA_SERVER_HEADER(GMSHPlugin_Algorithm)

#include "SMESH_3D_Algo_i.hxx"
#include "GMSHPlugin_GMSH.hxx"

// ======================================================
// GMSH 3d algorithm
// ======================================================
class GMSHPLUGIN_EXPORT GMSHPlugin_GMSH_i:
  public virtual POA_GMSHPlugin::GMSHPlugin_GMSH,
  public virtual SMESH_3D_Algo_i
{
public:
  // Constructor
  GMSHPlugin_GMSH_i( PortableServer::POA_ptr thePOA,
                              ::SMESH_Gen*            theGenImpl );
  // Destructor
  virtual ~GMSHPlugin_GMSH_i();
 
  // Get implementation
  ::GMSHPlugin_GMSH* GetImpl();
};

#endif
