// Copyright (C) 2012-2013  ALNEOS
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
#include "utilities.h"

#include "GMSHPlugin_Hypothesis_i.hxx"
#include "GMSHPlugin_Hypothesis_2D_i.hxx"
#include "GMSHPlugin_GMSH_i.hxx"
#include "GMSHPlugin_GMSH_2D_i.hxx"

template <class T> class GMSHPlugin_Creator_i:public HypothesisCreator_i<T>
{
  // as we have 'module GMSHPlugin' in GMSHPlugin_Algorithm.idl
  virtual std::string GetModuleName() { return "GMSHPlugin"; }
};

//=============================================================================
/*!
 *
 */
//=============================================================================

extern "C"
{
  GMSHPLUGIN_EXPORT
  GenericHypothesisCreator_i* GetHypothesisCreator (const char* aHypName)
  {
    MESSAGE("GetHypothesisCreator " << aHypName);

    GenericHypothesisCreator_i* aCreator = 0;

    // Hypotheses

    // Algorithms
    if (strcmp(aHypName, "GMSH") == 0)
      aCreator = new GMSHPlugin_Creator_i<GMSHPlugin_GMSH_i>;
    else if (strcmp(aHypName, "GMSH_2D") == 0)
      aCreator = new GMSHPlugin_Creator_i<GMSHPlugin_GMSH_2D_i>;
    // Hypotheses
    else if (strcmp(aHypName, "GMSH_Parameters") == 0)
      aCreator = new GMSHPlugin_Creator_i<GMSHPlugin_Hypothesis_i>;
    else if (strcmp(aHypName, "GMSH_Parameters_2D") == 0)
      aCreator = new GMSHPlugin_Creator_i<GMSHPlugin_Hypothesis_2D_i>;
    else ;

    return aCreator;
  }
}
