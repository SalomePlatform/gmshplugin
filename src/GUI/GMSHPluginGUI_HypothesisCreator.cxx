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
#include "GMSHPluginGUI_HypothesisCreator.h"

#include <SMESHGUI_Utils.h>
#include <SMESHGUI_HypothesesUtils.h>
#include <SMESHGUI_SpinBox.h>
#include <GeomSelectionTools.h>

#include CORBA_SERVER_HEADER(GMSHPlugin_Algorithm)

#include <SUIT_Session.h>
#include <SUIT_ResourceMgr.h>

#include <SalomeApp_Tools.h>
#include <LightApp_SelectionMgr.h>
#include <SALOME_ListIO.hxx>

#include <QComboBox>
#include <QLabel>
#include <QGroupBox>
#include <QFrame>
#include <QLayout>
#include <QLineEdit>
#include <QCheckBox>
#include <QPixmap>
#include <QTableWidget>
#include <QHeaderView>
#include <QPushButton>

enum Algo2D
  {
   automatic,
   meshadapt,
   delaunay,
   frontal,
   delaunayforquad,
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
   packingparallelograms,
   quadqs
#else
   packingparallelograms
#endif
  };

  enum Algo3D
  {
   delaunay3,
   frontal3,
   mmg3d,
   rtree,
   hxt
  };

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
  enum Recomb2DAlgo
  {
   simple,
   blossom,
   simplefullquads,
   blossomfullquads
  };
#else
  enum Recomb2DAlgo
  {
   standard,
   blossom
  };
#endif

enum SubdivAlgo
  {
   none,
   allquads,
   allhexas
  };

enum RemeshAlgo
  {
   nosplit,
   automaticR,
   automaticmetis
  };

enum RemeshPara
  {
   harmonic,
   conformal,
   rbfharmonic
  };

enum VerbLvl
  {
   silent,
   errors,
   warnings,
   direct,
   information,
   status,
   debug
  };


GMSHPluginGUI_HypothesisCreator::GMSHPluginGUI_HypothesisCreator( const QString& theHypType )
  : SMESHGUI_GenericHypothesisCreator( theHypType )
{
  myGeomSelectionTools = NULL;
  myCompoundSet.clear();
  myIs2D = ( theHypType.endsWith("2D"));
  myIs3D = (theHypType.endsWith("3D"));
}

GMSHPluginGUI_HypothesisCreator::~GMSHPluginGUI_HypothesisCreator()
{
}

bool GMSHPluginGUI_HypothesisCreator::checkParams(QString& msg) const
{
  GmshHypothesisData data_old, data_new;
  readParamsFromHypo( data_old );
  readParamsFromWidgets( data_new );
  bool res = storeParamsToHypo( data_new );
  msg = data_new.myErrorMsg;
  storeParamsToHypo( data_old );
  return res;
}

QFrame* GMSHPluginGUI_HypothesisCreator::buildFrame()
{
  QFrame* fr = new QFrame( 0 );
  fr->setObjectName( "myframe" );
  QVBoxLayout* lay = new QVBoxLayout( fr );
  lay->setMargin( 5 );
  lay->setSpacing( 0 );

  QTabWidget* tab = new QTabWidget( fr );
  tab->setTabShape( QTabWidget::Rounded );
  tab->setTabPosition( QTabWidget::North );
  lay->addWidget( tab );
  QWidget* GroupC1 = new QWidget();
  tab->insertTab( 0, GroupC1, tr( "SMESH_ARGUMENTS" ) );

  QGridLayout* aGroupLayout = new QGridLayout( GroupC1 );
  aGroupLayout->setSpacing( 6 );
  aGroupLayout->setMargin( 11 );

  int row = 0;
  myName = 0;
  if( isCreation() )
  {
    aGroupLayout->addWidget( new QLabel( tr( "SMESH_NAME" ), GroupC1 ), row, 0 );
    myName = new QLineEdit( GroupC1 );
    myName->setMinimumWidth(160);
    aGroupLayout->addWidget( myName, row, 1 );
    row++;
  }

  my2DAlgo = 0;
  if (!myIs3D)
  {
    aGroupLayout->addWidget(new QLabel(tr("GMSH_2D_ALGO"), GroupC1), row, 0);
    my2DAlgo = new QComboBox(GroupC1);
    QStringList types2DAlgo;
    types2DAlgo << tr("GMSH_AUTOMATIC")
      << tr("GMSH_MESH_ADAPT")
      << tr("GMSH_DELAUNAY")
      << tr("GMSH_FRONTAL")
      << tr("GMSH_DELAUNAY_FOR_QUAD")
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
      << tr("GMSH_PACKING_OF_PARALLELOGRAMS")
      << tr("GMSH_QUASI_STRUCTURED_QUAD");
#else
      << tr("GMSH_PACKING_OF_PARALLELOGRAMS");
#endif
    my2DAlgo->addItems(types2DAlgo);
    aGroupLayout->addWidget(my2DAlgo, row, 1);
    row++;
  }

  my3DAlgo = 0;
  if ( !myIs2D )
  {
    aGroupLayout->addWidget( new QLabel( tr( "GMSH_3D_ALGO" ), GroupC1 ), row, 0 );
    my3DAlgo = new QComboBox( GroupC1 );
    QStringList types3DAlgo;
    types3DAlgo << tr( "GMSH_DELAUNAY3" ) << tr( "GMSH_FRONTAL_DELAUNAY" ) << tr( "GMSH_MMG3D" ) <<
	           tr( "GMSH_R_TREE" ) << tr( "GMSH_HXT" );
    my3DAlgo->addItems( types3DAlgo );
    aGroupLayout->addWidget( my3DAlgo, row, 1 );
    row++;
  }

  myRecomb2DAlgo = 0;
  if (!myIs3D)
  {
    aGroupLayout->addWidget(new QLabel(tr("GMSH_2D_RECOMB_ALGO"), GroupC1), row, 0);
    myRecomb2DAlgo = new QComboBox(GroupC1);
    QStringList typesRecomb2DAlgo;
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=8
    typesRecomb2DAlgo << tr("GMSH_SIMPLE") << tr("GMSH_BLOSSOM") << tr("GMSH_SIMPLE_FULL_QUADS") << tr("GMSH_BLOSSOM_FULL_QUADS");
#else
    typesRecomb2DAlgo << tr("GMSH_STANDARD") << tr("GMSH_BLOSSOM");
#endif
    myRecomb2DAlgo->addItems(typesRecomb2DAlgo);
    aGroupLayout->addWidget(myRecomb2DAlgo, row, 1);
    row++;
  }
  myRecombineAll = new QCheckBox( tr( "GMSH_RECOMBINE_ALL" ), GroupC1 );
  aGroupLayout->addWidget( myRecombineAll, row, 0 );
  row++;

  mySubdivAlgo = 0;

  aGroupLayout->addWidget(new QLabel(tr("GMSH_SUBDIV_ALGO"), GroupC1), row, 0);
  mySubdivAlgo = new QComboBox(GroupC1);
  QStringList typesSubdivAlgo;
  typesSubdivAlgo << tr("GMSH_NONE") << tr("GMSH_ALL_QUADS") << tr("GMSH_ALL_HEXAS");
  mySubdivAlgo->addItems(typesSubdivAlgo);
  aGroupLayout->addWidget(mySubdivAlgo, row, 1);
  row++;

  myRemeshAlgo = 0;
  myRemeshPara = 0;
  if (!myIs3D)
  {
    aGroupLayout->addWidget(new QLabel(tr("GMSH_REMESH_ALGO"), GroupC1), row, 0);
    myRemeshAlgo = new QComboBox(GroupC1);
    QStringList typesRemeshAlgo;
    typesRemeshAlgo << tr("GMSH_NO_SPLIT") << tr("GMSH_AUTO") << tr("GMSH_AUTO_ONLY_WITH_METIS");
    myRemeshAlgo->addItems(typesRemeshAlgo);
    aGroupLayout->addWidget(myRemeshAlgo, row, 1);
    row++;

    aGroupLayout->addWidget(new QLabel(tr("GMSH_REMESH_PARA"), GroupC1), row, 0);
    myRemeshPara = new QComboBox(GroupC1);
    QStringList typesRemeshPara;
    typesRemeshPara << tr("GMSH_HARMONIC") << tr("GMSH_CONFORMAL") << tr("GMSH_RBF_HARMONIC");
    myRemeshPara->addItems(typesRemeshPara);
    aGroupLayout->addWidget(myRemeshPara, row, 1);
    row++;
  }

  aGroupLayout->addWidget( new QLabel( tr( "GMSH_SMOOTHING_STEPS" ), GroupC1 ), row, 0 );
  mySmouthSteps = new SMESHGUI_SpinBox( GroupC1 );
  mySmouthSteps->RangeStepAndValidator( 1, 1000, 1, "length_precision" );
  aGroupLayout->addWidget( mySmouthSteps, row, 1 );
  row++;

  aGroupLayout->addWidget( new QLabel( tr( "GMSH_SIZE_FACTOR" ), GroupC1 ), row, 0 );
  mySizeFactor = new SMESHGUI_SpinBox( GroupC1 );
  mySizeFactor->RangeStepAndValidator( 1e-06, 1e+06, 0.1, "length_precision" );
  aGroupLayout->addWidget( mySizeFactor, row, 1 );
  row++;

#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
  aGroupLayout->addWidget( new QLabel( tr( "GMSH_SIZE_FROM_CURVATURE" ), GroupC1 ), row, 0 );
  myMeshCurvatureSize = new SMESHGUI_SpinBox( GroupC1 );
  myMeshCurvatureSize->RangeStepAndValidator( 0.0, 1e+22, 1.0, "length_precision" );
  aGroupLayout->addWidget( myMeshCurvatureSize, row, 1 );
  row++;
#endif

  aGroupLayout->addWidget( new QLabel( tr( "GMSH_MIN_SIZE" ), GroupC1 ), row, 0 );
  myMinSize = new SMESHGUI_SpinBox( GroupC1 );
  myMinSize->RangeStepAndValidator( 0.0, 1e+22, 1., "length_precision" );
  aGroupLayout->addWidget( myMinSize, row, 1 );
  row++;

  aGroupLayout->addWidget( new QLabel( tr( "GMSH_MAX_SIZE" ), GroupC1 ), row, 0 );
  myMaxSize = new SMESHGUI_SpinBox( GroupC1 );
  myMaxSize->RangeStepAndValidator( 0.0, 1e+22, 1e+21, "length_precision" );
  aGroupLayout->addWidget( myMaxSize, row, 1 );
  row++;

  mySecondOrder = new QCheckBox( tr( "GMSH_SECOND_ORDER" ), GroupC1 );
  aGroupLayout->addWidget( mySecondOrder, row, 0 );

  myUseIncomplElem = new QCheckBox( tr( "GMSH_USE_INCOMPLETE_ELEMENT" ), GroupC1 );
  aGroupLayout->addWidget( myUseIncomplElem, row, 1 );
  row++;

  connect( mySecondOrder, SIGNAL( toggled( bool ) ), this, SLOT( updateWidgets() ) );

  aGroupLayout->addWidget( new QLabel( tr( "GMSH_VERB_LVL" ), GroupC1 ), row, 0 );
  myVerbLvl = new QComboBox(GroupC1);
  QStringList typesVerbLvl;
  typesVerbLvl << tr("GMSH_SILENT") << tr("GMSH_ERRORS") << tr("GMSH_WARNINGS") << tr("GMSH_DIRECT") << tr("GMSH_INFORMATION") << tr("GMSH_STATUS") << tr("GMSH_DEBUG");
  myVerbLvl->addItems(typesVerbLvl);
  aGroupLayout->addWidget(myVerbLvl, row, 1);
  row++;

  // Compounds
  if (!myIs3D)
  {
    QWidget* compoundGroup = new QWidget();
    tab->insertTab(1, compoundGroup, tr("GMSH_COMPOUND"));

    myCompoundTable = new QTableWidget(0, 2, compoundGroup);
    QGridLayout* compoundLayout = new QGridLayout(compoundGroup);
    compoundLayout->addWidget(myCompoundTable, 1, 0, 8, 1);

    QStringList compoundHeaders;
    compoundHeaders << tr("GMSH_COMPOUND_ENTRY_COLUMN") << tr("GMSH_COMPOUND_NAME_COLUMN");
    myCompoundTable->setHorizontalHeaderLabels(compoundHeaders);
    myCompoundTable->horizontalHeader()->hideSection(0);
    myCompoundTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
    myCompoundTable->resizeColumnToContents(1);
    myCompoundTable->setAlternatingRowColors(true);
    myCompoundTable->verticalHeader()->hide();

    QPushButton* addCompoundButton = new QPushButton(tr("GMSH_COMPOUND_ADD"), compoundGroup);
    compoundLayout->addWidget(addCompoundButton, 1, 1, 1, 1);
    QFrame *line2 = new QFrame(compoundGroup);

    line2->setFrameShape(QFrame::HLine);
    line2->setFrameShadow(QFrame::Sunken);
    compoundLayout->addWidget(line2, 2, 1, 1, 1);

    QPushButton* removeButton = new QPushButton(tr("GMSH_COMPOUND_REMOVE"), compoundGroup);
    compoundLayout->addWidget(removeButton, 3, 1, 1, 1);

    connect(addCompoundButton, SIGNAL(clicked()), this, SLOT(onAddCompound()));
    connect(removeButton, SIGNAL(clicked()), this, SLOT(onRemoveCompound()));
  }
  return fr;
}

void GMSHPluginGUI_HypothesisCreator::updateWidgets()
{
  myUseIncomplElem->setEnabled(mySecondOrder->isChecked());
}

void GMSHPluginGUI_HypothesisCreator::onAddCompound()
{
  GMSHPlugin::GMSHPlugin_Hypothesis_var h = GMSHPlugin::GMSHPlugin_Hypothesis::_narrow(initParamsHypothesis());
  GeomSelectionTools* geomSelectionTools = getGeomSelectionTools();
  LightApp_SelectionMgr* mySel = geomSelectionTools->selectionMgr();
  SALOME_ListIO ListSelectedObjects;
  mySel->selectedObjects(ListSelectedObjects, NULL, false );
  SALOME_ListIteratorOfListIO Object_It(ListSelectedObjects);
  for ( ; Object_It.More() ; Object_It.Next())
  {
    Handle(SALOME_InteractiveObject) anObject = Object_It.Value();
    std::string entry, shapeName;
    entry = geomSelectionTools->getEntryOfObject(anObject);
    shapeName = anObject->getName();
    TopAbs_ShapeEnum shapeType;
    shapeType = geomSelectionTools->entryToShapeType(entry);
    if ((shapeType == TopAbs_SHAPE) || (shapeType != TopAbs_EDGE && shapeType != TopAbs_FACE))
      continue;
    myCompoundTable->setFocus();
    QString shapeEntry;
    shapeEntry = QString::fromStdString(entry);
    if (myCompoundSet.contains(shapeEntry))
      continue;
    int row = myCompoundTable->rowCount() ;
    myCompoundTable->setRowCount(row+1);
    myCompoundTable->setItem(row, 0, new QTableWidgetItem(shapeEntry));
    myCompoundTable->item(row, 0 )->setFlags(0);
    myCompoundTable->setItem(row, 1, new QTableWidgetItem(QString::fromStdString(shapeName)));
    myCompoundTable->item(row, 1 )->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);
    myCompoundTable->resizeColumnToContents(1);
    myCompoundTable->clearSelection();
    myCompoundTable->scrollToItem( myCompoundTable->item( row, 1 ) );
    myCompoundSet.insert(shapeEntry);
    myCompoundToRemove.remove(shapeEntry);
  }
}

void GMSHPluginGUI_HypothesisCreator::onRemoveCompound()
{
  QList<int> selectedRows;
  QList<QTableWidgetItem*> selected = myCompoundTable->selectedItems();
  QTableWidgetItem* item;
  int row;
  foreach(item, selected)
  {
    row = item->row();
    if (!selectedRows.contains(row))
      selectedRows.append( row );
  }
  qSort( selectedRows );
  QListIterator<int> it( selectedRows );
  it.toBack();
  while (it.hasPrevious())
  {
    row = it.previous();
    QString entry = myCompoundTable->item(row,0)->text();
    if (myCompoundSet.contains(entry))
    {
      myCompoundSet.remove(entry);
      myCompoundToRemove.insert(entry);
    }
    myCompoundTable->removeRow(row );
  }
  myCompoundTable->resizeColumnToContents(1);
}

void GMSHPluginGUI_HypothesisCreator::retrieveParams() const
{
  GmshHypothesisData data;
  readParamsFromHypo( data );

  if( myName )
    myName->setText( data.myName );
  if (!myIs3D)
    my2DAlgo->setCurrentIndex( data.my2DAlgo );
  if ( !myIs2D )
    my3DAlgo->setCurrentIndex( data.my3DAlgo );
  if(!myIs3D)
    myRecomb2DAlgo->setCurrentIndex( data.myRecomb2DAlgo );
  if ( myRecombineAll )
    myRecombineAll->setChecked( data.myRecombineAll );
  if ( mySubdivAlgo )
    mySubdivAlgo->setCurrentIndex( data.mySubdivAlgo );
  if (!myIs3D)
  {
    myRemeshAlgo->setCurrentIndex(data.myRemeshAlgo);
    myRemeshPara->setCurrentIndex(data.myRemeshPara);
  }
  if(data.mySmouthStepsVar.isEmpty())
    mySmouthSteps->setValue( data.mySmouthSteps );
  else
    mySmouthSteps->setText( data.mySmouthStepsVar );
  if(data.mySizeFactorVar.isEmpty())
    mySizeFactor->setValue( data.mySizeFactor );
  else
    mySizeFactor->setText( data.mySizeFactorVar );
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
  if(data.myMeshCurvatureSizeVar.isEmpty())
    myMeshCurvatureSize->setValue( data.myMeshCurvatureSize );
  else
    myMeshCurvatureSize->setText( data.myMeshCurvatureSizeVar );
#endif
  if(data.myMaxSizeVar.isEmpty())
    myMaxSize->setValue( data.myMaxSize );
  else
    myMaxSize->setText( data.myMaxSizeVar );
  if(data.myMinSizeVar.isEmpty())
    myMinSize->setValue( data.myMinSize );
  else
    myMinSize->setText( data.myMinSizeVar );
  if ( mySecondOrder )
    mySecondOrder->setChecked( data.mySecondOrder );
  if ( myUseIncomplElem )
    myUseIncomplElem->setChecked( data.myUseIncomplElem );
  if (myVerbLvl)
    myVerbLvl->setCurrentIndex(data.myVerbLvl);

  GMSHPluginGUI_HypothesisCreator* that = (GMSHPluginGUI_HypothesisCreator*)this;
  that->updateWidgets();

  if (!myIs3D)
  {
    GeomSelectionTools* geomSelectionTools = that->getGeomSelectionTools();
    for (QSet<QString>::const_iterator i = myCompoundSet.begin(); i != myCompoundSet.end(); ++i)
    {
      const QString entry = *i;
      std::string shapeName = geomSelectionTools->getNameFromEntry(entry.toStdString());
      int row = myCompoundTable->rowCount();
      myCompoundTable->setRowCount(row + 1);
      myCompoundTable->setItem(row, 0, new QTableWidgetItem(entry));
      myCompoundTable->item(row, 0)->setFlags(0);
      myCompoundTable->setItem(row, 1, new QTableWidgetItem(QString::fromStdString(shapeName)));
      myCompoundTable->item(row, 1)->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);
    }
    myCompoundTable->resizeColumnToContents(1);
  }
}

QString GMSHPluginGUI_HypothesisCreator::storeParams() const
{
  GmshHypothesisData data;
  readParamsFromWidgets( data );
  storeParamsToHypo( data );

  QString valStr = tr("GMSH_MAX_SIZE") + " = " + QString::number( data.myMaxSize ) + "; ";
  valStr += tr("GMSH_MIN_SIZE") + " = " + QString::number( data.myMinSize ) + "; ";
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
  valStr += tr("GMSH_SIZE_FROM_CURVATURE") + " = " + QString::number( data.myMeshCurvatureSize ) + "; ";
#endif
  if ( data.mySecondOrder )
    valStr +=  tr("GMSH_SECOND_ORDER") + "; ";

  return valStr;
}

bool GMSHPluginGUI_HypothesisCreator::readParamsFromHypo( GmshHypothesisData& h_data ) const
{
  GMSHPlugin::GMSHPlugin_Hypothesis_var h =
    GMSHPlugin::GMSHPlugin_Hypothesis::_narrow( initParamsHypothesis() );

  HypothesisData* data = SMESH::GetHypothesisData( hypType() );
  h_data.myName = isCreation() && data ? data->Label : "";

  h_data.my2DAlgo = (int) h->Get2DAlgo();
  if ( !myIs2D )
    h_data.my3DAlgo = (int) h->Get3DAlgo();
  h_data.myRecomb2DAlgo = (int) h->GetRecomb2DAlgo();
  h_data.myRecombineAll = h->GetRecombineAll();
  h_data.mySubdivAlgo = (int) h->GetSubdivAlgo();
  h_data.myRemeshAlgo = (int) h->GetRemeshAlgo();
  h_data.myRemeshPara = (int) h->GetRemeshPara();
  h_data.mySmouthSteps = h->GetSmouthSteps();
  h_data.mySizeFactor = h->GetSizeFactor();
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
  h_data.myMeshCurvatureSize = h->GetMeshCurvatureSize();
  h_data.myMeshCurvatureSizeVar = getVariableName("SetMeshCurvatureSize");
#endif
  h_data.myMinSize = h->GetMinSize();
  h_data.myMaxSize = h->GetMaxSize();
  h_data.mySmouthStepsVar = getVariableName("SmouthSteps");
  h_data.mySizeFactorVar = getVariableName("SizeFactor");
  h_data.myMinSizeVar = getVariableName("SetMinSize");
  h_data.myMaxSizeVar = getVariableName("SetMaxSize");
  h_data.mySecondOrder = h->GetSecondOrder();
  h_data.myUseIncomplElem = h->GetUseIncomplElem();
  h_data.myVerbLvl = (int) h->GetVerbosityLevel();

  GMSHPluginGUI_HypothesisCreator* that = (GMSHPluginGUI_HypothesisCreator*)this;
  GMSHPlugin::string_array_var myEntries = h->GetCompoundOnEntries();
  for ( CORBA::ULong i=0 ; i<myEntries->length() ; i++ )
    {
      QString entry = myEntries[i].in();
      that->myCompoundSet.insert(entry);
    }

  return true;
}

bool GMSHPluginGUI_HypothesisCreator::storeParamsToHypo( const GmshHypothesisData& h_data ) const
{
  GMSHPlugin::GMSHPlugin_Hypothesis_var h =
    GMSHPlugin::GMSHPlugin_Hypothesis::_narrow( hypothesis() );

  bool ok = true;
  try
  {
    if( isCreation() )
      SMESH::SetName( SMESH::FindSObject( h ), h_data.myName.toLatin1().data() );
    if( !myIs3D )
      h->Set2DAlgo( h_data.my2DAlgo );
    if ( !myIs2D )
      h->Set3DAlgo( h_data.my3DAlgo );
    if (!myIs3D)
      h->SetRecomb2DAlgo( h_data.myRecomb2DAlgo );
    h->SetRecombineAll( h_data.myRecombineAll );
    h->SetSubdivAlgo( h_data.mySubdivAlgo );
    if (!myIs3D)
    {
      h->SetRemeshAlgo(h_data.myRemeshAlgo);
      h->SetRemeshPara(h_data.myRemeshPara);
    }
    h->SetSmouthSteps( h_data.mySmouthSteps );
    h->SetSizeFactor( h_data.mySizeFactor );
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
    h->SetMeshCurvatureSize( h_data.myMeshCurvatureSize );
#endif
    h->SetMinSize( h_data.myMinSize );
    h->SetMaxSize( h_data.myMaxSize );
    h->SetVarParameter( h_data.mySmouthStepsVar.toLatin1().constData(), "SmouthSteps");
    h->SetVarParameter( h_data.mySizeFactorVar.toLatin1().constData(), "SizeFactor");
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
    h->SetVarParameter( h_data.myMeshCurvatureSizeVar.toLatin1().constData(), "SetMeshCurvatureSize");
#endif
    h->SetVarParameter( h_data.myMinSizeVar.toLatin1().constData(), "SetMinSize");
    h->SetVarParameter( h_data.myMaxSizeVar.toLatin1().constData(), "SetMaxSize");
    h->SetSecondOrder( h_data.mySecondOrder );
    h->SetUseIncomplElem( h_data.myUseIncomplElem );
    h->SetIs2d( myIs2D );
    h->SetVerbosityLevel(h_data.myVerbLvl);

    QString mainEntry = getMainShapeEntry();
    for (QSet<QString>::const_iterator i = myCompoundSet.begin(); i != myCompoundSet.end(); ++i)
    {
      QString entry = *i;
      if ( myCompoundToRemove.contains( entry ))
        continue;
      if ( !mainEntry.isEmpty() && !entry.startsWith( mainEntry ))
      {
        h_data.myErrorMsg = "Compound group is not defined on the main geometry";
        ok = false;
        break;
      }
      h->SetCompoundOnEntry(entry.toLatin1().constData());
    }
    for (QSet<QString>::const_iterator i = myCompoundToRemove.begin(); i != myCompoundToRemove.end(); ++i)
    {
      const QString entry = *i;
      h->UnsetCompoundOnEntry(entry.toLatin1().constData());
    }
  }
  catch(const SALOME::SALOME_Exception& ex)
  {
    SalomeApp_Tools::QtCatchCorbaException(ex);
    ok = false;
  }
  return ok;
}

bool GMSHPluginGUI_HypothesisCreator::readParamsFromWidgets( GmshHypothesisData& h_data ) const
{
  h_data.myName           = myName ? myName->text() : "";
  if(my2DAlgo)
    h_data.my2DAlgo       = my2DAlgo->currentIndex();
  if (my3DAlgo)
    h_data.my3DAlgo       = my3DAlgo->currentIndex();
  if(myRecomb2DAlgo)
    h_data.myRecomb2DAlgo   = myRecomb2DAlgo->currentIndex();
  h_data.myRecombineAll   = myRecombineAll->isChecked();
  h_data.mySubdivAlgo     = mySubdivAlgo->currentIndex();
  if(myRemeshAlgo)
    h_data.myRemeshAlgo   = myRemeshAlgo->currentIndex();
  if(myRemeshPara)
    h_data.myRemeshPara   = myRemeshPara->currentIndex();
  h_data.mySmouthSteps    = mySmouthSteps->value();
  h_data.mySizeFactor     = mySizeFactor->value();
#if GMSH_MAJOR_VERSION >=4 && GMSH_MINOR_VERSION >=10
  h_data.myMeshCurvatureSize    = myMeshCurvatureSize->value();
  h_data.myMeshCurvatureSizeVar = myMeshCurvatureSize->text();
#endif
  h_data.myMinSize        = myMinSize->value();
  h_data.myMaxSize        = myMaxSize->value();
  h_data.mySmouthStepsVar = mySmouthSteps->text();
  h_data.mySizeFactorVar  = mySizeFactor->text();
  h_data.myMinSizeVar     = myMinSize->text();
  h_data.myMaxSizeVar     = myMaxSize->text();
  h_data.mySecondOrder    = mySecondOrder->isChecked();
  h_data.myUseIncomplElem = myUseIncomplElem->isChecked();
  h_data.myVerbLvl        = myVerbLvl->currentIndex();

  // ne semble pas utile dans la mesure ou myCompoundSet n'a pas besoin d'etre modifier
  /*
  GMSHPluginGUI_HypothesisCreator* that = (GMSHPluginGUI_HypothesisCreator*)this;
  int nbRows = myCompoundTable->rowCount();
  for(int row=0 ; row < nbRows ; row++)
  {
    QString entry = myLocalSizeTable->item(row, 1)->text();
    that->myCompoundSet.insert(entry);
  }
  */
  return true;
}

// on ne modifie rien à partir de là

GeomSelectionTools* GMSHPluginGUI_HypothesisCreator::getGeomSelectionTools()
{
  if (myGeomSelectionTools == NULL) {
    myGeomSelectionTools = new GeomSelectionTools();
  }
  return myGeomSelectionTools;
}

QString GMSHPluginGUI_HypothesisCreator::caption() const
{
  return tr( QString( "GMSH_%1_TITLE" ).arg(myIs2D?QString("2D"):QString("3D")).toLatin1().data() );
}

QPixmap GMSHPluginGUI_HypothesisCreator::icon() const
{
  QString hypIconName = tr( QString("ICON_DLG_GMSH_PARAMETERS%1").arg(myIs2D?QString("_2D"):QString("")).toLatin1().data() );
  return SUIT_Session::session()->resourceMgr()->loadPixmap( "GMSHPlugin", hypIconName );
}

QString GMSHPluginGUI_HypothesisCreator::type() const
{
  return tr( QString( "GMSH_%1_HYPOTHESIS" ).arg(myIs2D?QString("2D"):QString("3D")).toLatin1().data() );
}

QString GMSHPluginGUI_HypothesisCreator::helpPage() const
{
  return "gmsh_2d_3d_hypo_page.html";
}
