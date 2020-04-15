// Copyright (C) 2012-2015  ALNEOS
// Copyright (C) 2016-2020  EDF R&D
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
#ifndef _GMSHPlugin_Hypothesis_2D_i_HXX_
#define _GMSHPlugin_Hypothesis_2D_i_HXX_

#include "GMSHPlugin_Defs.hxx"

#include <SALOMEconfig.h>
#include CORBA_SERVER_HEADER(GMSHPlugin_Algorithm)

#include "GMSHPlugin_Hypothesis_i.hxx"
#include "GMSHPlugin_Hypothesis_2D.hxx"

class SMESH_Gen;

// GMSHPlugin parameters hypothesis (2D case)

class GMSHPLUGIN_EXPORT  GMSHPlugin_Hypothesis_2D_i:
  public virtual POA_GMSHPlugin::GMSHPlugin_Hypothesis_2D,
  public GMSHPlugin_Hypothesis_i
{
 public:
  // Constructor
  GMSHPlugin_Hypothesis_2D_i (PortableServer::POA_ptr thePOA,
                                ::SMESH_Gen*            theGenImpl);
  // Destructor
  virtual ~GMSHPlugin_Hypothesis_2D_i();
  
  // Verify whether hypothesis supports given entity type 
  CORBA::Boolean IsDimSupported( SMESH::Dimension type );
};

#endif
