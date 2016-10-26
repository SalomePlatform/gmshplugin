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
   delaunayforquad
  };

enum Algo3D
  {
   frontal3,
   frontaldelaunay,
   fontalhex,
   mmg3d,
   rtree
  };

enum Recomb2DAlgo
  {
   standard,
   blossom
  };

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


GMSHPluginGUI_HypothesisCreator::GMSHPluginGUI_HypothesisCreator( const QString& theHypType )
  : SMESHGUI_GenericHypothesisCreator( theHypType )
{
  myGeomSelectionTools = NULL;
  myCompoundSet.clear();
  myIs2D = ( theHypType.endsWith("2D"));
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

  aGroupLayout->addWidget( new QLabel( tr( "GMSH_2D_ALGO" ), GroupC1 ), row, 0 );
  my2DAlgo = new QComboBox( GroupC1 );
  QStringList types2DAlgo;
  types2DAlgo << tr( "GMSH_AUTOMATIC" ) << tr( "GMSH_MESH_ADAPT" )   << tr( "GMSH_DELAUNAY" ) <<
                 tr( "GMSH_FRONTAL" )       << tr( "GMSH_DELAUNAY_FOR_QUAD" );
  my2DAlgo->addItems( types2DAlgo );
  aGroupLayout->addWidget( my2DAlgo, row, 1 );
  row++;
  
  my3DAlgo = 0;
  if ( !myIs2D )
  {
    aGroupLayout->addWidget( new QLabel( tr( "GMSH_3D_ALGO" ), GroupC1 ), row, 0 );
    my3DAlgo = new QComboBox( GroupC1 );
    QStringList types3DAlgo;
    types3DAlgo << tr( "GMSH_FRONTAL_DELAUNAY" ) << tr( "GMSH_FRONTAL_HEX" )   << tr( "GMSH_MMG3D" ) <<
                   tr( "GMSH_R_TREE" );
    my3DAlgo->addItems( types3DAlgo );
    aGroupLayout->addWidget( my3DAlgo, row, 1 );
    row++;
  }
  
  aGroupLayout->addWidget( new QLabel( tr( "GMSH_2D_RECOMB_ALGO" ), GroupC1 ), row, 0 );
  myRecomb2DAlgo = new QComboBox( GroupC1 );
  QStringList typesRecomb2DAlgo;
  typesRecomb2DAlgo << tr( "GMSH_STANDARD" ) << tr( "GMSH_BLOSSOM" );
  myRecomb2DAlgo->addItems( typesRecomb2DAlgo );
  aGroupLayout->addWidget( myRecomb2DAlgo, row, 1 );
  row++;
  
  myRecombineAll = new QCheckBox( tr( "GMSH_RECOMBINE_ALL" ), GroupC1 );
  aGroupLayout->addWidget( myRecombineAll, row, 0 );
  row++;
  
  aGroupLayout->addWidget( new QLabel( tr( "GMSH_SUBDIV_ALGO" ), GroupC1 ), row, 0 );
  mySubdivAlgo = new QComboBox( GroupC1 );
  QStringList typesSubdivAlgo;
  typesSubdivAlgo << tr( "GMSH_NONE" ) << tr( "GMSH_ALL_QUADS" )   << tr( "GMSH_ALL_HEXAS" );
  mySubdivAlgo->addItems( typesSubdivAlgo );
  aGroupLayout->addWidget( mySubdivAlgo, row, 1 );
  row++;
  
  aGroupLayout->addWidget( new QLabel( tr( "GMSH_REMESH_ALGO" ), GroupC1 ), row, 0 );
  myRemeshAlgo = new QComboBox( GroupC1 );
  QStringList typesRemeshAlgo;
  typesRemeshAlgo << tr( "GMSH_NO_SPLIT" ) << tr( "GMSH_AUTO" )   << tr( "GMSH_AUTO_ONLY_WITH_METIS" );
  myRemeshAlgo->addItems( typesRemeshAlgo );
  aGroupLayout->addWidget( myRemeshAlgo, row, 1 );
  row++;
  
  aGroupLayout->addWidget( new QLabel( tr( "GMSH_REMESH_PARA" ), GroupC1 ), row, 0 );
  myRemeshPara = new QComboBox( GroupC1 );
  QStringList typesRemeshPara;
  typesRemeshPara << tr( "GMSH_HARMONIC" ) << tr( "GMSH_CONFORMAL" )   << tr( "GMSH_RBF_HARMONIC" );
  myRemeshPara->addItems( typesRemeshPara );
  aGroupLayout->addWidget( myRemeshPara, row, 1 );
  row++;
  
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
  
  // Compounds
  QWidget* compoundGroup = new QWidget();
  tab->insertTab(1, compoundGroup, tr("GMSH_COMPOUND"));
  
  myCompoundTable = new QTableWidget(0, 2, compoundGroup);
  QGridLayout* compoundLayout = new QGridLayout(compoundGroup);
  compoundLayout->addWidget(myCompoundTable, 1, 0, 8, 1);
  
  QStringList compoundHeaders;
  compoundHeaders << tr( "GMSH_COMPOUND_ENTRY_COLUMN" ) << tr( "GMSH_COMPOUND_NAME_COLUMN" );
  myCompoundTable->setHorizontalHeaderLabels(compoundHeaders);
  myCompoundTable->horizontalHeader()->hideSection(0);
#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
  myCompoundTable->horizontalHeader()->setResizeMode(QHeaderView::Interactive);
#else
  myCompoundTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
#endif
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

  connect( addCompoundButton, SIGNAL(clicked()), this, SLOT(onAddCompound()));
  connect( removeButton, SIGNAL(clicked()), this, SLOT(onRemoveCompound()));
  
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
  for (Object_It ; Object_It.More() ; Object_It.Next())
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
  my2DAlgo->setCurrentIndex( data.my2DAlgo );
  if ( !myIs2D )
    my3DAlgo->setCurrentIndex( data.my3DAlgo );
  myRecomb2DAlgo->setCurrentIndex( data.myRecomb2DAlgo );
  if ( myRecombineAll )
    myRecombineAll->setChecked( data.myRecombineAll );
  if ( mySubdivAlgo )
  mySubdivAlgo->setCurrentIndex( data.mySubdivAlgo );
  myRemeshAlgo->setCurrentIndex( data.myRemeshAlgo);
  myRemeshPara->setCurrentIndex( data.myRemeshPara);
  if(data.mySmouthStepsVar.isEmpty())
    mySmouthSteps->setValue( data.mySmouthSteps );
  else
    mySmouthSteps->setText( data.mySmouthStepsVar );
  if(data.mySizeFactorVar.isEmpty())
    mySizeFactor->setValue( data.mySizeFactor );
  else
    mySizeFactor->setText( data.mySizeFactorVar );
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
  
  GMSHPluginGUI_HypothesisCreator* that = (GMSHPluginGUI_HypothesisCreator*)this;
  that->updateWidgets();
  
  GeomSelectionTools* geomSelectionTools = that->getGeomSelectionTools();
  for (QSet<QString>::const_iterator i = myCompoundSet.begin(); i != myCompoundSet.end(); ++i)
  {
    const QString entry = *i;
    std::string shapeName = geomSelectionTools->getNameFromEntry(entry.toStdString());
    int row = myCompoundTable->rowCount();
    myCompoundTable->setRowCount(row+1);
    myCompoundTable->setItem(row, 0, new QTableWidgetItem(entry));
    myCompoundTable->item(row, 0)->setFlags(0);
    myCompoundTable->setItem(row, 1, new QTableWidgetItem(QString::fromStdString(shapeName)));
    myCompoundTable->item(row, 1)->setFlags(Qt::ItemIsSelectable|Qt::ItemIsEnabled);
  }
  myCompoundTable->resizeColumnToContents(1);
}

QString GMSHPluginGUI_HypothesisCreator::storeParams() const
{
  GmshHypothesisData data;
  readParamsFromWidgets( data );
  storeParamsToHypo( data );
  
  QString valStr = tr("GMSH_MAX_SIZE") + " = " + QString::number( data.myMaxSize ) + "; ";
  valStr += tr("GMSH_MIN_SIZE") + " = " + QString::number( data.myMinSize ) + "; ";
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
  h_data.myMinSize = h->GetMinSize();
  h_data.myMaxSize = h->GetMaxSize();
  h_data.mySmouthStepsVar = getVariableName("SmouthSteps");
  h_data.mySizeFactorVar = getVariableName("SizeFactor");
  h_data.myMinSizeVar = getVariableName("SetMinSize");
  h_data.myMaxSizeVar = getVariableName("SetMaxSize");
  h_data.mySecondOrder = h->GetSecondOrder();
  h_data.myUseIncomplElem = h->GetUseIncomplElem();
  
  GMSHPluginGUI_HypothesisCreator* that = (GMSHPluginGUI_HypothesisCreator*)this;
  GMSHPlugin::string_array_var myEntries = h->GetCompoundOnEntries();
  for ( int i=0 ; i<myEntries->length() ; i++ )
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

    h->Set2DAlgo( h_data.my2DAlgo );
    if ( !myIs2D )
      h->Set3DAlgo( h_data.my3DAlgo );
    h->SetRecomb2DAlgo( h_data.myRecomb2DAlgo );
    h->SetRecombineAll( h_data.myRecombineAll );
    h->SetSubdivAlgo( h_data.mySubdivAlgo );
    h->SetRemeshAlgo( h_data.myRemeshAlgo );
    h->SetRemeshPara( h_data.myRemeshPara );
    h->SetSmouthSteps( h_data.mySmouthSteps );
    h->SetSizeFactor( h_data.mySizeFactor );
    h->SetMinSize( h_data.myMinSize );
    h->SetMaxSize( h_data.myMaxSize );
    h->SetVarParameter( h_data.mySmouthStepsVar.toLatin1().constData(), "SmouthSteps");
    h->SetVarParameter( h_data.mySizeFactorVar.toLatin1().constData(), "SizeFactor");
    h->SetVarParameter( h_data.myMinSizeVar.toLatin1().constData(), "SetMinSize");
    h->SetVarParameter( h_data.myMaxSizeVar.toLatin1().constData(), "SetMaxSize");
    h->SetSecondOrder( h_data.mySecondOrder );
    h->SetUseIncomplElem( h_data.myUseIncomplElem );
    h->SetIs2d( myIs2D );
    
    for (QSet<QString>::const_iterator i = myCompoundSet.begin(); i != myCompoundSet.end(); ++i)
    {
      const QString entry = *i;
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
  h_data.my2DAlgo         = my2DAlgo->currentIndex();
  if (my3DAlgo)
    h_data.my3DAlgo       = my3DAlgo->currentIndex();
  h_data.myRecomb2DAlgo   = myRecomb2DAlgo->currentIndex();
  h_data.myRecombineAll   = myRecombineAll->isChecked();
  h_data.mySubdivAlgo     = mySubdivAlgo->currentIndex();
  h_data.myRemeshAlgo     = myRemeshAlgo->currentIndex();
  h_data.myRemeshPara     = myRemeshPara->currentIndex();
  h_data.mySmouthSteps    = mySmouthSteps->value();
  h_data.mySizeFactor     = mySizeFactor->value();
  h_data.myMinSize        = myMinSize->value();
  h_data.myMaxSize        = myMaxSize->value();
  h_data.mySmouthStepsVar = mySmouthSteps->text();
  h_data.mySizeFactorVar  = mySizeFactor->text();
  h_data.myMinSizeVar     = myMinSize->text();
  h_data.myMaxSizeVar     = myMaxSize->text();
  h_data.mySecondOrder    = mySecondOrder->isChecked();
  h_data.myUseIncomplElem = myUseIncomplElem->isChecked();
  
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
  _PTR(Study) aStudy = SMESH::GetActiveStudyDocument();
  if (myGeomSelectionTools == NULL || myGeomSelectionTools->getMyStudy() != aStudy) {
    delete myGeomSelectionTools;
    myGeomSelectionTools = new GeomSelectionTools(aStudy);
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
