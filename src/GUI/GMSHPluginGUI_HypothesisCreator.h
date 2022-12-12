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
#ifndef GMSHPLUGINGUI_HypothesisCreator_HeaderFile
#define GMSHPLUGINGUI_HypothesisCreator_HeaderFile

#include "GmshVersion.h"
#include "GMSHPluginGUI.h"

#include <SMESHGUI_Hypotheses.h>

#include <TopAbs_ShapeEnum.hxx>

#include <QSet>

class SMESHGUI_SpinBox;
class GeomSelectionTools;
class QComboBox;
class QCheckBox;
class QLineEdit;
class QTableWidget;

typedef struct
{
  QString             myName;
  int                 my2DAlgo,my3DAlgo,myRecomb2DAlgo;
  bool                myRecombineAll;
  int                 mySubdivAlgo,myRemeshAlgo,myRemeshPara,mySmouthSteps;
  bool                myUseIncomplElem;
  bool                mySecondOrder;
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
  double              mySizeFactor,myMaxSize, myMinSize, myMeshCurvatureSize;
  QString             myMaxSizeVar, myMinSizeVar, mySmouthStepsVar, mySizeFactorVar, myMeshCurvatureSizeVar;
#else
  double              mySizeFactor,myMaxSize, myMinSize;
  QString             myMaxSizeVar, myMinSizeVar, mySmouthStepsVar, mySizeFactorVar;
#endif
  mutable QString     myErrorMsg;
} GmshHypothesisData;

/*!
 * \brief Class for creation of GMSH2D and GMSH3D hypotheses
*/
class GMSHPLUGIN_GUI_EXPORT GMSHPluginGUI_HypothesisCreator : public SMESHGUI_GenericHypothesisCreator
{
  Q_OBJECT

public:
  GMSHPluginGUI_HypothesisCreator( const QString& );
  virtual ~GMSHPluginGUI_HypothesisCreator();

  virtual bool     checkParams(QString& msg) const;
  virtual QString  helpPage() const;

protected:
  virtual QFrame*  buildFrame    ();
  virtual void     retrieveParams() const;
  virtual QString  storeParams   () const;
  
  virtual QString  caption() const;
  virtual QPixmap  icon() const;
  virtual QString  type() const;

protected slots:
  void               updateWidgets();
  virtual void       onAddCompound();
  virtual void       onRemoveCompound();
  
private:
  bool readParamsFromHypo( GmshHypothesisData& ) const;
  bool readParamsFromWidgets( GmshHypothesisData& ) const;
  bool storeParamsToHypo( const GmshHypothesisData& ) const;
  GeomSelectionTools* getGeomSelectionTools();

private:
 QLineEdit*        myName;
 QComboBox*        my2DAlgo;
 QComboBox*        my3DAlgo;
 QComboBox*        myRecomb2DAlgo;
 QCheckBox*        myRecombineAll;
 QComboBox*        mySubdivAlgo;
 QComboBox*        myRemeshAlgo;
 QComboBox*        myRemeshPara;
 SMESHGUI_SpinBox* mySmouthSteps;
 SMESHGUI_SpinBox* mySizeFactor;
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
 SMESHGUI_SpinBox* myMeshCurvatureSize;
#endif
 SMESHGUI_SpinBox* myMaxSize;
 SMESHGUI_SpinBox* myMinSize;
 QCheckBox*        myUseIncomplElem;
 QCheckBox*        mySecondOrder;
 bool myIs2D;
 bool myIs3D;

 QTableWidget* myCompoundTable;
 GeomSelectionTools* myGeomSelectionTools;
 QSet<QString> myCompoundSet;
 QSet<QString> myCompoundToRemove;
};

#endif
