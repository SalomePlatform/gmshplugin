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
#ifndef _GMSHPlugin_Hypothesis_i_HXX_
#define _GMSHPlugin_Hypothesis_i_HXX_

#include "GMSHPlugin_Defs.hxx"

#include <SALOMEconfig.h>
#include CORBA_SERVER_HEADER(GMSHPlugin_Algorithm)

#include "SMESH_Hypothesis_i.hxx"
#include "GMSHPlugin_Hypothesis.hxx"

class SMESH_Gen;
//class GEOM_Object;

// GMSHPlugin parameters hypothesis

class GMSHPLUGIN_EXPORT GMSHPlugin_Hypothesis_i:
  public virtual POA_GMSHPlugin::GMSHPlugin_Hypothesis,
  public virtual SMESH_Hypothesis_i
{
 public:
  // Constructor
  GMSHPlugin_Hypothesis_i (PortableServer::POA_ptr thePOA,
                             ::SMESH_Gen*            theGenImpl);
  // Destructor
  virtual ~GMSHPlugin_Hypothesis_i();
  
  // Ajout d'un truc
  
  void SetMaxSize(CORBA::Double theSize);
  CORBA::Double GetMaxSize();
  
  void SetMinSize(CORBA::Double theSize);
  CORBA::Double GetMinSize();
  
  void SetSecondOrder(CORBA::Boolean theVal);
  CORBA::Boolean GetSecondOrder();
  
  void Set2DAlgo(CORBA::Long the2DAlgo);
  CORBA::Long Get2DAlgo();
  void Set3DAlgo(CORBA::Long the3DAlgo);
  CORBA::Long Get3DAlgo();
  void SetRecomb2DAlgo(CORBA::Long theRecomb2DAlgo);
  CORBA::Long GetRecomb2DAlgo();
  void SetRecombineAll(CORBA::Boolean theRecombineAll);
  CORBA::Boolean GetRecombineAll();
  void SetSubdivAlgo(CORBA::Long theSubdivAlgo);
  CORBA::Long GetSubdivAlgo();
  void SetRemeshAlgo(CORBA::Long theRemeshAlgo);
  CORBA::Long GetRemeshAlgo();
  void SetRemeshPara(CORBA::Long theRemeshPara);
  CORBA::Long GetRemeshPara();
  void SetSmouthSteps(CORBA::Double theSmouthSteps);
  CORBA::Double GetSmouthSteps();
  void SetSizeFactor(CORBA::Double theSizeFactor);
  CORBA::Double GetSizeFactor();
  void SetUseIncomplElem(CORBA::Boolean theUseIncomplElem);
  CORBA::Boolean GetUseIncomplElem();
  void SetIs2d(CORBA::Boolean theIs2d);
  
  void SetCompoundOnShape(GEOM::GEOM_Object_ptr GeomObj);
  void SetCompoundOnEntry(const char* entry);
  void UnsetCompoundOnShape(GEOM::GEOM_Object_ptr GeomObj);
  void UnsetCompoundOnEntry(const char* entry);
  GMSHPlugin::string_array* GetCompoundOnEntries();
  
  // fin d'ajout

  void SetGrowthRate(CORBA::Double theRate);
  CORBA::Double GetGrowthRate();

  void SetNbSegPerEdge(CORBA::Double theVal);
  CORBA::Double GetNbSegPerEdge();

  void SetNbSegPerRadius(CORBA::Double theVal);
  CORBA::Double GetNbSegPerRadius();

  // void SetLocalSizeOnShape(GEOM::GEOM_Object_ptr GeomObj, CORBA::Double localSize);
  // void SetLocalSizeOnEntry(const char* entry, CORBA::Double localSize);
  // CORBA::Double GetLocalSizeOnEntry(const char* entry);
  // GMSHPlugin::string_array* GetLocalSizeEntries();
  
  void UnsetLocalSizeOnEntry(const char* entry);

  // Get implementation
  ::GMSHPlugin_Hypothesis* GetImpl();
  
  // Verify whether hypothesis supports given entity type 
  CORBA::Boolean IsDimSupported( SMESH::Dimension type );

 protected:

  // to remember whether a parameter is already set (issue 0021364)
  enum SettingMethod
  {
    METH_SetMaxSize          = 1,
    METH_SetMinSize          = 2,
    METH_SetSecondOrder      = 4,

    //METH_SetFineness         = 16,
    METH_SetGrowthRate       = 32,
    METH_SetNbSegPerEdge     = 64,
    METH_SetNbSegPerRadius   = 128,
    METH_SetLocalSizeOnEntry = 256,
    // ajout ici
    METH_SetCompoundOnEntry    = 299,
    METH_Set2DAlgo           = 300,
    METH_Set3DAlgo           = 301,
    METH_SetRecomb2DAlgo     = 302,
    METH_SetRecombineAll     = 303,
    METH_SetSubdivAlgo       = 304,
    METH_SetRemeshAlgo       = 305,
    METH_SetRemeshPara       = 306,
    METH_SetSmouthSteps      = 307,
    METH_SetSizeFactor       = 308,
    METH_SetUseIncomplElem   = 309,
    METH_SetIs2d             = 310,
    // fin d'ajout
    METH_LAST                = METH_SetLocalSizeOnEntry
  };
  int mySetMethodFlags;

  // Return true if a parameter is not yet set, else return true if a parameter changes.
  // PythonDumping depends on the result of this function.
  // Checking only change of a parameter is not enough because then the default values are
  // not dumped and if the defaults will change then the behaviour of scripts
  // created without dump of the default parameters will also change what is not good.
  template<typename T>
    bool isToSetParameter(T curValue, T newValue, /*SettingMethod*/int meth)
  {
    if ( mySetMethodFlags & meth ) // already set, check if a value is changing
      return ( curValue != newValue );
    else
      return ( mySetMethodFlags |= meth ); // == return true
  }

 public:
  // method intended to remove explicit treatment of Netagen hypotheses from
  // SMESH_NoteBook to assure backward compatibility after implemeneting
  // issue 0021308: Remove hard-coded dependency of the external mesh plugins
  virtual int getParamIndex(const TCollection_AsciiString& method, int nbVars) const;

  // method used to convert variable parameters stored in an old study
  // into myMethod2VarParams. It should return a method name for an index of
  // variable parameters. Index is countered from zero
  virtual std::string getMethodOfParameter(const int paramIndex, int nbVars) const;


  // Methods for copying mesh definition to other geometry

  // Return geometry this hypothesis depends on. Return false if there is no geometry parameter
  virtual bool getObjectsDependOn( std::vector< std::string > & entryArray,
                                   std::vector< int >         & subIDArray ) const;

  // Set new geometry instead of that returned by getObjectsDependOn()
  virtual bool setObjectsDependOn( std::vector< std::string > & entryArray,
                                   std::vector< int >         & subIDArray );
};

#endif
