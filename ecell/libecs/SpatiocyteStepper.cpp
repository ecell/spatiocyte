//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 2006-2009 Keio University
//
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-Cell is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// E-Cell is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with E-Cell -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
// E-Cell Project, Institute for Advanced Biosciences, Keio University.
//


#include <time.h>
#include <gsl/gsl_randist.h>
#include <libecs/Model.hpp>
#include <libecs/System.hpp>
#include <libecs/Process.hpp>
#include <libecs/Stepper.hpp>
#include <libecs/VariableReference.hpp>
#include <libecs/SpatiocyteStepper.hpp>
#include <libecs/SpatiocyteSpecies.hpp>
#include <libecs/SpatiocyteThread.hpp>

LIBECS_DM_INIT_STATIC(SpatiocyteStepper, Stepper);

void SpatiocyteStepper::initialize()
{
  if(isInitialized)
    {
      return;
    }
  isInitialized = true;
  Stepper::initialize(); 
  if(theProcessVector.empty())
    {
      THROW_EXCEPTION(InitializationFailed,
                      getPropertyInterface().getClassName() + 
                      ": at least one Process must be defined in this" +
                      " Stepper.");
    } 
  std::cout << "1. checking model..." << std::endl;
  checkModel();
  initializeThreads();
  //We need a Comp tree to assign the voxels to each Comp
  //and get the available number of vacant voxels. The compartmentalized
  //vacant voxels are needed to randomly place molecules according to the
  //Comp:
  std::cout << "2. creating compartments..." << std::endl;
  registerComps();
  setCompsProperties();
  std::cout << "3. setting up lattice properties..." << std::endl;
  setLatticeProperties(); 
  broadcastLatticeProperties();
  setCompsCenterPoint();
  //All species have been created at this point, we initialize them now:
  std::cout << "4. initializing species..." << std::endl;
  initSpecies();
  initializeFirst();
  std::cout << "5. initializing processes the second time..." << std::endl;
  initializeSecond();
  std::cout << "6. constructing lattice..." << std::endl;
  constructLattice(0);
  concatenateLattice(0);
  setBoundaries();
  //checkLattice();
  std::cout << "7. setting intersecting compartment list..." << std::endl;
  setIntersectingCompartmentList();
  std::cout << "8. compartmentalizing lattice..." << std::endl;
  compartmentalizeLattice();
  std::cout << "9. setting up compartment voxels properties..." << std::endl;
  setCompVoxelProperties();
  resizeProcessLattice();
  std::cout << "10. initializing processes the third time..." << std::endl;
  initializeThird();
  std::cout << "11. printing simulation parameters..." << std::endl;
  updateSpecies();
  storeSimulationParameters();
  printSimulationParameters();
  std::cout << "12. populating compartments with molecules..." << std::endl;
  populateComps();
  std::cout << "13. initializing processes the fourth time..." << std::endl;
  initializeFourth();
  std::cout << "14. initializing the priority queue..." << std::endl;
  initPriorityQueue();
  std::cout << "15. initializing processes the fifth time..." << std::endl;
  initializeFifth();
  std::cout << "16. initializing processes the last time..." << std::endl;
  initializeLastOnce();
  std::cout << "17. finalizing species..." << std::endl;
  finalizeSpecies();
  std::cout << "18. printing final process parameters..." << std::endl <<
    std::endl;
  printProcessParameters();
  std::cout << "19. simulation is started..." << std::endl;
}

SpatiocyteStepper::~SpatiocyteStepper()
{
      std::cout << "destructing" << std::endl;
  for(unsigned i(0); i != theThreads.size(); ++i)
    {
      std::cout << "destructing" << i << std::endl;
      //theThreads[i]->joinChildren();
    }
}

void SpatiocyteStepper::initializeThreads()
{
  nThreadsRunning = 0;
  flagA = FLAG_STOP;
  flagB = FLAG_STOP;
  __sync_synchronize();
  theThreads.resize(ThreadSize);
  for(unsigned i(0); i != ThreadSize; ++i)
    {
      theThreads[i] = new Thread(i, ThreadSize, nThreadsRunning, flagA,
                                 flagB, theSpecies, *this);
      if(i)
        {
          theThreads[i]->create();
        }
    }
}

void SpatiocyteStepper::interrupt(double aTime)
{
  setCurrentTime(aTime); 
  for(std::vector<Process*>::const_iterator 
      i(theExternInterruptedProcesses.begin());
      i != theExternInterruptedProcesses.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->substrateValueChanged(getCurrentTime()); 
    }
  setNextTime(thePriorityQueue.getTop()->getTime());
}

void SpatiocyteStepper::finalizeSpecies()
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      (*i)->finalizeSpecies();
    }
  theThreads[0]->initialize();
  theThreads[0]->initializeLists();
}

unsigned SpatiocyteStepper::getBoxSize()
{
  return theTotalBoxSize;
}

void SpatiocyteStepper::updateSpecies()
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      (*i)->updateSpecies();
    }
}

unsigned SpatiocyteStepper::getRowSize()
{
  return theTotalRows;
}

unsigned SpatiocyteStepper::getLayerSize()
{
  return theTotalLayers;
}

unsigned SpatiocyteStepper::getColSize()
{
  return theTotalCols;
}

unsigned SpatiocyteStepper::getLatticeSize()
{
  return theIDs[0].size()*theTotalBoxSize;
}

Point SpatiocyteStepper::getCenterPoint()
{
  return theCenterPoint;
} 

double SpatiocyteStepper::getNormalizedVoxelRadius()
{
  return nVoxelRadius;
}

void SpatiocyteStepper::reset(int seed)
{
  gsl_rng_set(getRng(), seed); 
  setCurrentTime(0);
  initializeSecond();
  clearComps();
  initializeThird();
  populateComps();
  initializeFourth();
  initPriorityQueue();
  initializeFifth();
  finalizeSpecies();
  //checkLattice();
}

Species* SpatiocyteStepper::addSpecies(Variable* aVariable)
{
  std::vector<Species*>::iterator aSpeciesIter(variable2ispecies(aVariable));
  if(aSpeciesIter == theSpecies.end())
    {
      Species *aSpecies(new Species(this, aVariable, theSpecies.size(),
                          getRng(), VoxelRadius, theIDs, theInfo, theAdjoins,
                          theAdjBoxes, theAdjAdjBoxes));
      theSpecies.push_back(aSpecies);
      return aSpecies;
    }
  return *aSpeciesIter;
}

Species* SpatiocyteStepper::getSpecies(Variable* aVariable)
{
  std::vector<Species*>::iterator aSpeciesIter(variable2ispecies(aVariable));
  if(aSpeciesIter == theSpecies.end())
    {
      return NULL;
    }
  return *aSpeciesIter;
}

std::vector<Species*> SpatiocyteStepper::getSpecies()
{
  return theSpecies;
}

bool SpatiocyteStepper::isBoundaryMol(unsigned aMol,
                                        unsigned aDimension)
{
  //This method is only for checking boundaries on a cube:
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aMol, aRow, aLayer, aCol);
  if(aDimension == 3)
    {
      //If the voxel is on one of the 6 cuboid surfaces:
      if(aRow == 1 || aLayer == 1 || aCol == 1 ||
         aRow == theTotalRows-2 || aLayer == theTotalLayers-2 || 
         aCol == theTotalCols-2)
        {
          return true;
        }
    }
  //If it is a surface voxel:
  else
    {
      //If the voxel is on one of the 12 edges of the cube:
      if((aRow <= 1 && aCol <= 1) ||
         (aRow <= 1 && aLayer <= 1) ||
         (aCol <= 1 && aLayer <= 1) ||
         (aRow <= 1 && aCol >= theTotalCols-2) ||
         (aRow <= 1 && aLayer >= theTotalLayers-2) ||
         (aCol <= 1 && aRow >= theTotalRows-2) ||
         (aCol <= 1 && aLayer >= theTotalLayers-2) ||
         (aLayer <= 1 && aRow >= theTotalRows-2) ||
         (aLayer <= 1 && aCol >= theTotalCols-2) ||
         (aRow >= theTotalRows-2 && aCol >= theTotalCols-2) ||
         (aRow >= theTotalRows-2 && aLayer >= theTotalLayers-2) ||
         (aCol >= theTotalCols-2 && aLayer >= theTotalLayers-2))
        {
          return true;
        }
    }
  return false;
}

unsigned SpatiocyteStepper::getPeriodicMol(unsigned aMol,
                                                 unsigned aDimension,
                                                 Origin* anOrigin)
{
  //This method is only for checking boundaries on a cube:
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aMol, aRow, aLayer, aCol);
  unsigned nextRow(aRow);
  unsigned nextLayer(aLayer);
  unsigned nextCol(aCol);
  unsigned adj(1);
  if(aDimension == 3)
    {
      adj = 0;
    }
  if(aRow == 1+adj)
    {
      nextRow = theTotalRows-(3+adj);
      --anOrigin->row;
    }
  else if(aRow == theTotalRows-(2+adj))
    {
      nextRow = 2+adj;
      ++anOrigin->row;
    }
  if(aLayer == 1+adj)
    {
      nextLayer = theTotalLayers-(3+adj);
      --anOrigin->layer;
    }
  else if(aLayer == theTotalLayers-(2+adj))
    {
      nextLayer = 2+adj;
      ++anOrigin->layer;
    }
  if(aCol == 1+adj)
    {
      nextCol = theTotalCols-(3+adj);
      --anOrigin->col;
    }
  else if(aCol == theTotalCols-(2+adj))
    {
      nextCol = 2+adj;
      ++anOrigin->col;
    }
  if(nextRow != aRow || nextCol != aCol || nextLayer != aLayer)
    {
      return nextRow+theTotalRows*nextLayer+theTotalRows*theTotalLayers*nextCol;
    }
  return theNullMol;
}

Point SpatiocyteStepper::getPeriodicPoint(unsigned aMol,
                                          unsigned aDimension,
                                          Origin* anOrigin)
{
  unsigned adj(1);
  if(aDimension == 3)
    {
      adj = 0;
    }
  unsigned row(theTotalRows-(3+adj)-(1+adj));
  unsigned layer(theTotalLayers-(3+adj)-(1+adj));
  unsigned col(theTotalCols-(3+adj)-(1+adj));
  unsigned aGlobalCol;
  unsigned aGlobalLayer;
  unsigned aGlobalRow;
  coord2global(aMol, aGlobalRow, aGlobalLayer, aGlobalCol);
  int aRow(aGlobalRow+anOrigin->row*row);
  int aLayer(aGlobalLayer+anOrigin->layer*layer);
  int aCol(aGlobalCol+anOrigin->col*col);
  Point aPoint;
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      aPoint.y = (aCol%2)*theHCPl+theHCPy*aLayer;
      aPoint.z = aRow*2*nVoxelRadius+
        ((aLayer+aCol)%2)*nVoxelRadius;
      aPoint.x = aCol*theHCPx;
      break;
    case CUBIC_LATTICE:
      aPoint.y = aLayer*2*nVoxelRadius;
      aPoint.z = aRow*2*nVoxelRadius;
      aPoint.x = aCol*2*nVoxelRadius;
      break;
    }
  return aPoint;
}


std::vector<Species*>::iterator
SpatiocyteStepper::variable2ispecies(Variable* aVariable)
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      if((*i)->getVariable() == aVariable)
        {
          return i;
        }
    }
  return theSpecies.end();
} 

Species* SpatiocyteStepper::variable2species(Variable* aVariable)
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      if((*i)->getVariable() == aVariable)
        {
          return (*i);
        }
    }
  return NULL;
}

void SpatiocyteStepper::checkModel()
{
  //check if nonHD species are being used by non-SpatiocyteProcesses
  Model::StepperMap aStepperMap(getModel()->getStepperMap());  
  for(Model::StepperMap::const_iterator i(aStepperMap.begin());
      i != aStepperMap.end(); ++i )
    {   
      if(i->second != this)
        {
          std::vector<Process*> aProcessVector(i->second->getProcessVector());
          for(std::vector<Process*>::const_iterator j(aProcessVector.begin());
              j != aProcessVector.end(); ++j)
            {
              Process::VariableReferenceVector aVariableReferenceVector( 
                               (*j)->getVariableReferenceVector());
              for(Process::VariableReferenceVector::const_iterator 
                  k(aVariableReferenceVector.begin());
                  k != aVariableReferenceVector.end(); ++k)
                { 
                  const VariableReference& aNewVariableReference(*k);
                  Variable* aVariable(aNewVariableReference.getVariable()); 
                  for(std::vector<Species*>::iterator m(theSpecies.begin());
                      m !=theSpecies.end(); ++m)
                    {
                      if((*m)->getVariable() == aVariable)
                        {
                          THROW_EXCEPTION(ValueError, 
                            getPropertyInterface().getClassName() +
                            ": " + aVariable->getFullID().asString()  +  
                            " is a non-HD species but it is being used " +
                            "by non-SpatiocyteProcess: " +
                            (*j)->getFullID().asString());
                        }
                    } 
                }
            }
        }
    }
}

void SpatiocyteStepper::checkLattice()
{
  unsigned cnt(0);
  for(unsigned i(0); i != theTotalBoxSize; ++i)
    {
      std::vector<unsigned>& anAdjoins(theAdjoins[i]);
      unsigned offset(i*theBoxMaxSize);
      for(unsigned j(0); j != theIDs[i].size();  ++j)
        { 
          for(unsigned l(0); l != theAdjoinSize; ++l)
            { 
              if(anAdjoins[j*theAdjoinSize+l] == j+offset)
                {
                  //theSpecies[1]->addMol(i, j);
                  cnt++;
                }
            }
        }
    }
  //theSpecies[1]->setIsPopulated();
  std::cout << "cnt:" << cnt << std::endl;

  /*
  std::vector<int> list;
  for(unsigned i(0); i!=theSpecies.size(); ++i)
    {
      list.push_back(0);
    }
  for(unsigned i(0); i!=theIDs.size(); ++i)
    {
      ++list[theIDs[i]];
    }
  int volumeCnt(0);
  int surfaceCnt(0);
  for(unsigned i(0); i!=list.size(); ++i)
    {
      std::cout << "i:" << i << " ";
      if(theSpecies[i]->getVariable() != NULL)
        {
          std::cout << theSpecies[i]->getVariable()->getFullID().asString();
          if(theSpecies[i]->getDimension() == 3)
            {
              volumeCnt += list[i];
            }
          else
            {
              surfaceCnt += list[i];
            }
        }
      std::cout << " cnt:" << list[i] << std::endl;
    }
  std::cout << "total volume:" << volumeCnt << std::endl;
  std::cout << "total surface:" << surfaceCnt << std::endl;
  std::cout << "total volume+surface:" << surfaceCnt+volumeCnt << std::endl;
  */
}

void SpatiocyteStepper::initSpecies()
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      (*i)->initialize(theSpecies.size(), theBoxMaxSize, theAdjoinSize,
                       theNullMol, theNullID, theThreads);
    }
  for(std::vector<Comp*>::const_iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      (*i)->interfaceID = theSpecies.size();
    }
}

void SpatiocyteStepper::broadcastLatticeProperties()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->setLatticeProperties(&theIDs, &theInfo, &theAdjoins,
                                     theAdjoinSize, theBoxMaxSize, theNullMol,
                                     theNullID, theThreads[0]);
    }
}

void SpatiocyteStepper::initializeFirst()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeFirst();
    }
}

void SpatiocyteStepper::initializeSecond()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeSecond();
    }
}

void SpatiocyteStepper::printProcessParameters()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->printParameters();
    }
  std::cout << std::endl;
}

void SpatiocyteStepper::resizeProcessLattice()
{
  /*
  unsigned startMol(theIDs.size());
  unsigned endMol(startMol);
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      endMol += aProcess->getLatticeResizeMol(endMol);
    }
  theIDs.resize(endMol);
  theInfo.resize(endMol);
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {    
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->updateResizedLattice();
    }
    */
}

void SpatiocyteStepper::initializeThird()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeThird();
    }
}

void SpatiocyteStepper::initializeFourth()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeFourth();
    }
}

void SpatiocyteStepper::initializeFifth()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeFifth();
    }
  setStepInterval(thePriorityQueue.getTop()->getTime()-getCurrentTime());
}

void SpatiocyteStepper::initializeLastOnce()
{
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      SpatiocyteProcess*
        aProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      aProcess->initializeLastOnce();
    }
}

void SpatiocyteStepper::initPriorityQueue()
{
  const double aCurrentTime(getCurrentTime());
  thePriorityQueue.clear();
  for(std::vector<Process*>::const_iterator i(theProcessVector.begin());
      i != theProcessVector.end(); ++i)
    {      
      Process* const aProcess(*i);
      SpatiocyteProcess*
        aSpatiocyteProcess(dynamic_cast<SpatiocyteProcess*>(*i));
      if(aSpatiocyteProcess)
        {
          aSpatiocyteProcess->setTime(aCurrentTime+
                                      aSpatiocyteProcess->getInterval());
          aSpatiocyteProcess->setPriorityQueue(&thePriorityQueue);
          //Not all SpatiocyteProcesses are inserted into the priority queue.
          //Only the following processes are inserted in the PriorityQueue and
          //executed at simulation steps according to their execution times:
          if(aSpatiocyteProcess->getIsPriorityQueued())
            {
              aSpatiocyteProcess->setQueueID(
                                   thePriorityQueue.push(aSpatiocyteProcess));
              //ExternInterrupted processes are processes which are
              //interrupted by non SpatiocyteStepper processes such as
              //ODE processes:
              if(aSpatiocyteProcess->getIsExternInterrupted())
                {
                  theExternInterruptedProcesses.push_back(aProcess);
                }
            }
        }
      //Interruption between SpatiocyteProcesses:
      //Only ReactionProcesses can interrupt other processes because only 
      //they can change the number of molecules. 
      //This method is called to set the list of processes which will be
      //interrupted by the current ReactionProcess:
      ReactionProcess*
        aReactionProcess(dynamic_cast<ReactionProcess*>(*i));
      if(aReactionProcess)
        {
          aReactionProcess->setInterruption(theProcessVector);
        }
    }
}

void SpatiocyteStepper::populateComps()
{
  for(std::vector<Comp*>::const_iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      populateComp(*i);
    }
}

void SpatiocyteStepper::clearComps()
{
  for(std::vector<Comp*>::const_iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      clearComp(*i);
    }
}


inline void SpatiocyteStepper::step()
{
  do
    {
      thePriorityQueue.getTop()->fire();
    }
  while(thePriorityQueue.getTop()->getTime() == getCurrentTime());
  setNextTime(thePriorityQueue.getTop()->getTime());
  //checkLattice();
} 


void SpatiocyteStepper::registerComps()
{
  System* aRootSystem(getModel()->getRootSystem());
  std::vector<Comp*> allSubs;
  //The root Comp is theComps[0]
  theComps.push_back(registerComp(aRootSystem, &allSubs));
  //After this we will create an species to get an ID to represent
  //NULL Comps. So let us specify the correct
  //size of the biochemical species to be simulated before additional
  //non-biochemical species are created:
  theBioSpeciesSize = theSpecies.size();
  //Create one last species to represent a NULL Comp. This is for
  //voxels that do not belong to any Comps:
  Species* aSpecies(new Species(this, NULL, theSpecies.size(), getRng(),
                                VoxelRadius, theIDs, theInfo, theAdjoins,
                                theAdjBoxes, theAdjAdjBoxes));
  theSpecies.push_back(aSpecies);
  aSpecies->setComp(NULL);
  theNullID = aSpecies->getID(); 
  //Expand the tree of immediate subComps into single list such that
  //the super Comps come first while the subComps 
  //come later in the list:
  std::vector<Comp*> Comps(theComps[0]->immediateSubs);
  while(!Comps.empty())
    {
      std::vector<Comp*> subComps;
      for(unsigned i(0); i != Comps.size(); ++i)
        {
          theComps.push_back(Comps[i]);
          for(unsigned j(0);
              j != Comps[i]->immediateSubs.size(); ++j)
            {
              subComps.push_back(Comps[i]->immediateSubs[j]);
            }
        }
      Comps = subComps;
    }
}

//allSubs contains all the subComps (child, grand child, great grand
//child, etc). Used to calculate the total number of Comp voxels.
Comp* SpatiocyteStepper::registerComp(System* aSystem,
                                            std::vector<Comp*>* allSubs)
{ 
  //We execute this function to register the System, and its subsystems
  //recursively.
  Comp* aComp(new Comp);
  aComp->minRow = UINT_MAX;
  aComp->minCol = UINT_MAX;
  aComp->minLayer = UINT_MAX;
  aComp->maxRow = 0;
  aComp->maxCol = 0;
  aComp->maxLayer = 0;
  aComp->lengthX = 0;
  aComp->lengthY = 0;
  aComp->lengthZ = 0;
  aComp->originX = 0;
  aComp->originY = 0;
  aComp->originZ = 0;
  aComp->rotateX = 0;
  aComp->rotateY = 0;
  aComp->rotateZ = 0;
  aComp->xyPlane = 0;
  aComp->xzPlane = 0;
  aComp->yzPlane = 0;
  aComp->specVolume = 0;
  aComp->system = aSystem;
  aComp->surfaceSub = NULL;
  aComp->diffusiveComp = NULL;
  //Default Comp geometry is Cuboid:
  aComp->geometry = 0;
  //Default is volume Comp:
  aComp->dimension = 3;
  if(getVariable(aSystem, "DIMENSION"))
    { 
      const double aDimension(aSystem->getVariable("DIMENSION")->getValue());
      aComp->dimension = aDimension;
    }
  if(aComp->dimension == 3)
    {
      if(getVariable(aSystem, "GEOMETRY"))
        { 
          aComp->geometry = aSystem->getVariable("GEOMETRY")->getValue();
        }
      if(getVariable(aSystem, "LENGTHX"))
        {
          aComp->lengthX = aSystem->getVariable("LENGTHX")->getValue();
        }
      if(getVariable(aSystem, "LENGTHY"))
        {
          aComp->lengthY = aSystem->getVariable("LENGTHY")->getValue();
        }
      if(getVariable(aSystem, "LENGTHZ"))
        {
          aComp->lengthZ = aSystem->getVariable("LENGTHZ")->getValue();
        }
      if(getVariable(aSystem, "ORIGINX"))
        {
          aComp->originX = aSystem->getVariable("ORIGINX")->getValue();
        }
      if(getVariable(aSystem, "ORIGINY"))
        {
          aComp->originY = aSystem->getVariable("ORIGINY")->getValue();
        }
      if(getVariable(aSystem, "ORIGINZ"))
        {
          aComp->originZ = aSystem->getVariable("ORIGINZ")->getValue();
        }
      if(getVariable(aSystem, "ROTATEX"))
        {
          aComp->rotateX = aSystem->getVariable("ROTATEX")->getValue();
        }
      if(getVariable(aSystem, "ROTATEY"))
        {
          aComp->rotateY = aSystem->getVariable("ROTATEY")->getValue();
        }
      if(getVariable(aSystem, "ROTATEZ"))
        {
          aComp->rotateZ = aSystem->getVariable("ROTATEZ")->getValue();
        }
      if(getVariable(aSystem, "XYPLANE"))
        {
          aComp->xyPlane = aSystem->getVariable("XYPLANE")->getValue();
        }
      if(getVariable(aSystem, "XZPLANE"))
        {
          aComp->xzPlane = aSystem->getVariable("XZPLANE")->getValue();
        }
      if(getVariable(aSystem, "YZPLANE"))
        {
          aComp->yzPlane = aSystem->getVariable("YZPLANE")->getValue();
        }
    }
  registerCompSpecies(aComp);
  //Systems contains all the subsystems of a System.
  //For example /membrane is the subsystem of /:
  FOR_ALL(System::Systems, aSystem->getSystems())
    {
      std::vector<Comp*> currSubs; 
      Comp* aSubComp(registerComp(i->second, &currSubs)); 
      allSubs->push_back(aSubComp);
      for(unsigned j(0); j != currSubs.size(); ++j)
        {
          allSubs->push_back(currSubs[j]);
        }
      aComp->immediateSubs.push_back(aSubComp);
      if(aSubComp->dimension == 2)
      {
          aSubComp->geometry = aComp->geometry;
          aSubComp->specVolume = aComp->specVolume;
          aSubComp->lengthX = aComp->lengthX;
          aSubComp->lengthY = aComp->lengthY;
          aSubComp->lengthZ = aComp->lengthZ;
          aSubComp->originX = aComp->originX;
          aSubComp->originY = aComp->originY;
          aSubComp->originZ = aComp->originZ;
          aSubComp->xyPlane = aComp->xyPlane;
          aSubComp->xzPlane = aComp->xzPlane;
          aSubComp->yzPlane = aComp->yzPlane;
          aSubComp->rotateX = aComp->rotateX;
          aSubComp->rotateY = aComp->rotateY;
          aSubComp->rotateZ = aComp->rotateZ;

          for(unsigned i(0); i != aSubComp->lineSubs.size(); ++i)
          {
              Comp* lineComp(aSubComp->lineSubs[i]);
              lineComp->geometry = aComp->geometry;
              lineComp->specVolume = aComp->specVolume;
              lineComp->lengthX = aComp->lengthX;
              lineComp->lengthY = aComp->lengthY;
              lineComp->lengthZ = aComp->lengthZ;
              lineComp->originX = aComp->originX;
              lineComp->originY = aComp->originY;
              lineComp->originZ = aComp->originZ;
              lineComp->xyPlane = aComp->xyPlane;
              lineComp->xzPlane = aComp->xzPlane;
              lineComp->yzPlane = aComp->yzPlane;
              lineComp->rotateX = aComp->rotateX;
              lineComp->rotateY = aComp->rotateY;
              lineComp->rotateZ = aComp->rotateZ;
          }

          aComp->surfaceSub = aSubComp;
      }
      else if(aSubComp->dimension == 1)
      {
          // Properties of aComp (dimension == 2) is not fully set here yet.
          aComp->lineSubs.push_back(aSubComp);
      }
    }
  aComp->allSubs = *allSubs;
  return aComp;
}

void SpatiocyteStepper::setCompsProperties()
{
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      setCompProperties(theComps[i]);
    }
}

void SpatiocyteStepper::setCompsCenterPoint()
{
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      setCompCenterPoint(theComps[i]);
    }
}

void SpatiocyteStepper::registerCompSpecies(Comp* aComp)
{
  System* aSystem(aComp->system);
  FOR_ALL(System::Variables, aSystem->getVariables())
    {
      Variable* aVariable(i->second);
      if(aVariable->getID() == "VACANT")
        {
          aComp->enclosed = aVariable->getValue();
          //Set the number of vacant molecules to be always 0 because
          //when we populate lattice we shouldn't create more vacant
          //molecules than the ones already created for the Comp:
          aVariable->setValue(0);
          Species* aSpecies(addSpecies(aVariable));
          aComp->vacantSpecies = aSpecies;
          aComp->vacantID = aSpecies->getID(); //remove this
          aSpecies->setVacantSpecies(aSpecies);
          aSpecies->setIsCompVacant();
        }
      std::vector<Species*>::iterator j(variable2ispecies(aVariable));
      if(j != theSpecies.end())
        {
          aComp->species.push_back(*j);
          (*j)->setComp(aComp);
          (*j)->setDimension(aComp->dimension);
        }
    }
}

void SpatiocyteStepper::setLatticeProperties()
{
  Comp* aRootComp(theComps[0]);
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      theAdjoinSize = 12;
      theHCPl = nVoxelRadius/sqrt(3); 
      theHCPx = nVoxelRadius*sqrt(8.0/3); //Lx
      theHCPy = nVoxelRadius*sqrt(3); //Ly
      break;
    case CUBIC_LATTICE:
      theAdjoinSize = 6;
      break;
    }
  if(aRootComp->geometry == CUBOID)
    {
      //We do not give any leeway space between the simulation boundary
      //and the cell boundary if it is CUBOID to support
      //periodic boundary conditions:
      theCenterPoint.z = aRootComp->lengthZ/2; //row
      theCenterPoint.y = aRootComp->lengthY/2; //layer
      theCenterPoint.x = aRootComp->lengthX/2; //column
    }
  else
    {
      switch(LatticeType)
        {
        case HCP_LATTICE: 
          theCenterPoint.z = aRootComp->lengthZ/2+4*
            nVoxelRadius; //row
          theCenterPoint.y = aRootComp->lengthY/2+2*theHCPy; //layer
          theCenterPoint.x = aRootComp->lengthX/2+2*theHCPx; //column
          break;
        case CUBIC_LATTICE:
          theCenterPoint.z = aRootComp->lengthZ/2+8*nVoxelRadius;
          theCenterPoint.y = aRootComp->lengthY/2+8*nVoxelRadius;
          theCenterPoint.x = aRootComp->lengthX/2+8*nVoxelRadius;
          break;
        }
    }
  aRootComp->centerPoint = theCenterPoint; 
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      theTotalRows = (unsigned)rint((theCenterPoint.z)/
                                      (nVoxelRadius));
      theTotalLayers = (unsigned)rint((theCenterPoint.y*2)/theHCPy);
      theTotalCols = (unsigned)rint((theCenterPoint.x*2)/theHCPx);
      break;
    case CUBIC_LATTICE:
      theTotalRows = (unsigned)rint((theCenterPoint.z)/
                                      (nVoxelRadius));
      theTotalLayers = (unsigned)rint((theCenterPoint.y)/
                                      (nVoxelRadius));
      theTotalCols = (unsigned)rint((theCenterPoint.x)/
                                      (nVoxelRadius));
      break;
    }
  //For the CUBOID cell geometry, we need to readjust the size of
  //row, layer and column according to the boundary condition of its surfaces
  //to reflect the correct volume. This is because periodic boundary will
  //consume a layer of the surface voxels:
  if(aRootComp->geometry == CUBOID)
    {
      //We need to increase the row, layer and col size by 2 because
      //the entire volume must be surrounded by nullID voxels to avoid
      //self-homodimerization reaction.
      theTotalRows += 2;
      theTotalCols += 2;
      theTotalLayers += 2;
      readjustSurfaceBoundarySizes();
    }
  theBoxRows = 5;
  theBoxCols = 8;
  theBoxLayers = 6;
  theTotalBoxSize = theBoxRows*theBoxCols*theBoxLayers;
  theBoxSize = theTotalBoxSize/ThreadSize;
  theRows.resize(theTotalBoxSize);
  theCols.resize(theTotalBoxSize);
  theLayers.resize(theTotalBoxSize);
  theBoxMaxSize = 0;
  unsigned aBoxRowSize(theTotalRows/theBoxRows);
  unsigned aBoxColSize(theTotalCols/theBoxCols);
  unsigned aBoxLayerSize(theTotalLayers/theBoxLayers);
  if(theBoxRows > 1 && aBoxRowSize%2 != 0)
    {
      if((aBoxRowSize+1)*(theBoxRows-1) < theTotalRows)
        {
          aBoxRowSize++;
        }
      else
        {
          aBoxRowSize--;
        }
    }
  if(theBoxLayers > 1 && aBoxLayerSize%2 != 0)
    {
      if((aBoxLayerSize+1)*(theBoxLayers-1) < theTotalLayers)
        {
          aBoxLayerSize++;
        }
      else
        {
          aBoxLayerSize--;
        }
    }
  if(theBoxCols > 1 && aBoxColSize%2 != 0)
    {
      if((aBoxColSize+1)*(theBoxCols-1) < theTotalCols)
        {
          aBoxColSize++;
        }
      else
        {
          aBoxColSize--;
        }
    }
  for(unsigned i(0); i != theBoxRows; ++i)
    {
      unsigned rowSize(aBoxRowSize);
      if(i == theBoxRows-1)
        {
          rowSize = theTotalRows-(aBoxRowSize*(theBoxRows-1));
        }
      for(unsigned j(0); j != theBoxLayers; ++j)
        {
          unsigned layerSize(aBoxLayerSize);
          if(j == theBoxLayers-1)
            {
              layerSize = theTotalLayers-(aBoxLayerSize*(theBoxLayers-1));
            }
          for(unsigned k(0); k != theBoxCols; ++k)
            { 
              unsigned colSize(aBoxColSize);
              if(k == theBoxCols-1)
                {
                  colSize = theTotalCols-(aBoxColSize*(theBoxCols-1)); 
                }
              const unsigned box(i+
                                 theBoxRows*j+
                                 theBoxRows*theBoxLayers*k);
              theRows[box] = rowSize;
              theLayers[box] = layerSize;
              theCols[box] = colSize;
            }
        }
    }
  setAdjBoxes();
  setAdjAdjBoxes();
  const unsigned box(theBoxRows-1+
                     theBoxRows*(theBoxLayers-1)+
                     theBoxRows*theBoxLayers*(theBoxCols-1));
  /*
  std::cout << "aBoxRowSize;" << theRows[box] << std::endl;
  std::cout << "aBoxColSize;" << theCols[box] << std::endl;
  std::cout << "aBoxLayerSize;" << theLayers[box] << std::endl;
  */
  theBoxMaxSize = 2*std::max(aBoxRowSize, theRows[box])*
    std::max(aBoxColSize, theCols[box])*std::max(aBoxLayerSize, theLayers[box]);
  double max(theTotalRows*theTotalCols*theTotalLayers);
  if(max*2 > UINT_MAX)
    {
      std::cout << "ERROR: too many voxels (" << max*2 << ") more than " <<
        UINT_MAX << std::endl;
    }

  /*
  unsigned totalVoxels(theTotalRows*theTotalCols*theTotalLayers);
  unsigned boxVoxels(totalVoxels/theTotalBoxSize);
  unsigned idealSides(pow(boxVoxels, 1.0/3.0));
  theRows[i] = idealSides;
  theCols[i] = idealSides;
  theLayers[i] = idealSides;
  if(theTotalRows <= idealSize)
    {
      theRows[i] = theTotalRows;
    }
  if(theTotalCols <= idealSize)
    {
      theCols[i] = theTotalCols;
    }
  if(theTotalLayers <= idealSize)
    {
      theLayers[i] = theTotalLayers;
    }
    */

  theIDs.resize(theTotalBoxSize);
  theAdjoins.resize(theTotalBoxSize);
  theInfo.resize(theTotalBoxSize);
  theNullMol = 0;
  //Initialize the null coord:
}

void SpatiocyteStepper::setAdjAdjBoxes()
{
  theAdjAdjBoxes.resize(theTotalBoxSize);
  for(unsigned i(0); i != theTotalBoxSize; ++i)
    {
      const unsigned boxA(i);
      const std::vector<unsigned> adjBoxesA(theAdjBoxes[boxA]);
      std::vector<unsigned> adjAdjA(theAdjAdjBoxes[boxA]);
      for(unsigned j(0); j != adjBoxesA.size(); ++j)
        {
          const unsigned boxB(adjBoxesA[j]);
          const std::vector<unsigned> adjBoxesB(theAdjBoxes[boxB]);
          std::vector<unsigned> adjAdjB(theAdjAdjBoxes[boxB]);
          for(unsigned k(0); k != adjBoxesB.size(); ++k)
            {
              const unsigned boxC(adjBoxesB[k]);
              if(std::find(adjBoxesA.begin(), adjBoxesA.end(), boxC) !=
                 adjBoxesA.end())
                {
                  if(std::find(adjAdjA.begin(), adjAdjA.end(), boxC) ==
                     adjAdjA.end())
                    {
                      adjAdjA.push_back(boxC);
                    }
                  if(std::find(adjAdjB.begin(), adjAdjB.end(), boxC) ==
                     adjAdjB.end())
                    {
                      adjAdjB.push_back(boxC);
                    }
                }
            }
        }
      std::cout << "box:" << boxA << std::endl;
      for(unsigned j(0); j != adjAdjA.size(); ++j)
        {
          std::cout << "     adjAdj:" << adjAdjA[j] << std::endl;
        }
    }
}


void SpatiocyteStepper::setAdjBoxes()
{
  Comp* aRootComp(theComps[0]);
  theAdjBoxes.resize(theTotalBoxSize);
  for(unsigned i(0); i != theTotalBoxSize; ++i)
    {
      const unsigned bc(i/(theBoxRows*theBoxLayers)); 
      const unsigned bl((i%(theBoxRows*theBoxLayers))/theBoxRows); 
      const unsigned br((i%(theBoxRows*theBoxLayers))%theBoxRows); 
      std::vector<unsigned> cols;
      std::vector<unsigned> rows;
      std::vector<unsigned> layers;
      if(std::find(cols.begin(), cols.end(), bc) == cols.end())
        {
          cols.push_back(bc);
        }
      if(theBoxCols != 1)
        {
          if(!bc)
            {
              if(aRootComp->geometry == CUBOID &&
                 aRootComp->yzPlane == PERIODIC)
                {
                  if(std::find(cols.begin(), cols.end(), theBoxCols-1) ==
                     cols.end())
                    {
                      cols.push_back(theBoxCols-1);
                    }
                }
            }
          else
            {
              if(std::find(cols.begin(), cols.end(), bc-1) == cols.end())
                {
                  cols.push_back(bc-1);
                }
            }
          if(bc == theBoxCols-1)
            {
              if(aRootComp->geometry == CUBOID &&
                 aRootComp->yzPlane == PERIODIC && theBoxCols > 2)
                {
                  if(std::find(cols.begin(), cols.end(), 0) == cols.end())
                    {
                      cols.push_back(0);
                    }
                }
            }
          else
            {
              if(std::find(cols.begin(), cols.end(), bc+1) == cols.end())
                {
                  cols.push_back(bc+1);
                }
            }
        }
      if(std::find(rows.begin(), rows.end(), br) == rows.end())
        {
          rows.push_back(br);
        }
      if(theBoxRows != 1)
        {
          if(!br)
            {
              if(aRootComp->geometry == CUBOID &&
                 aRootComp->xyPlane == PERIODIC)
                {
                  if(std::find(rows.begin(), rows.end(), theBoxRows-1) ==
                     rows.end())
                    {
                      rows.push_back(theBoxRows-1);
                    }
                }
            }
          else
            {
              if(std::find(rows.begin(), rows.end(), br-1) == rows.end())
                {
                  rows.push_back(br-1);
                }
            }
          if(br == theBoxRows-1)
            {
              if(aRootComp->geometry == CUBOID &&
                 aRootComp->xyPlane == PERIODIC && theBoxRows > 2)
                {
                  if(std::find(rows.begin(), rows.end(), 0) == rows.end())
                    {
                      rows.push_back(0);
                    }
                }
            }
          else
            {
              if(std::find(rows.begin(), rows.end(), br+1) == rows.end())
                {
                  rows.push_back(br+1);
                }
            }
        }
      if(std::find(layers.begin(), layers.end(), bl) == layers.end())
        {
          layers.push_back(bl);
        }
      if(theBoxLayers != 1)
        {
          if(!bl)
            {
              if(aRootComp->geometry == CUBOID &&
                 aRootComp->xzPlane == PERIODIC)
                {
                  if(std::find(layers.begin(), layers.end(), theBoxLayers-1) ==
                     layers.end())
                    {
                      layers.push_back(theBoxLayers-1);
                    }
                }
            }
          else
            {
              if(std::find(layers.begin(), layers.end(), bl-1) == layers.end())
                {
                  layers.push_back(bl-1);
                }
            }
          if(bl == theBoxLayers-1)
            {
              if(aRootComp->geometry == CUBOID &&
                 aRootComp->xzPlane == PERIODIC && theBoxLayers > 2)
                {
                  if(std::find(layers.begin(), layers.end(), 0) == layers.end())
                    {
                      layers.push_back(0);
                    }
                }
            }
          else
            {
              if(std::find(layers.begin(), layers.end(), bl+1) == layers.end())
                {
                  layers.push_back(bl+1);
                }
            }
        }
      for(unsigned j(0); j != rows.size(); ++j)
        {
          const unsigned aRow(rows[j]);
          for(unsigned k(0); k != cols.size(); ++k)
            {
              const unsigned aCol(cols[k]);
              for(unsigned l(0); l != layers.size(); ++l)
                {
                  const unsigned aLayer(layers[l]); 
                  const unsigned box(aRow+
                                     theBoxRows*aLayer+
                                     theBoxRows*theBoxLayers*aCol);
                  if(box != i)
                    {
                      theAdjBoxes[i].push_back(box);
                    }
                }
            }
        }
    }
}


void SpatiocyteStepper::storeSimulationParameters()
{
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      Comp* aComp(theComps[i]); 
      if(aComp->dimension == 2)
        {
          aComp->actualArea =  (72*pow(VoxelRadius,2))*
            aComp->vacantSpecies->size()/(6*pow(2,0.5)+4*pow(3,0.5)+
                                         3*pow(6, 0.5));
          setSystemSize(aComp->system, aComp->actualArea*1e+2);
        }
      else // (aComp->dimension == 3)
        { 
          int voxelCnt(aComp->vacantSpecies->size());
          for(unsigned j(0); j != aComp->allSubs.size(); ++j)
            {
              voxelCnt += aComp->allSubs[j]->vacantSpecies->size();
            }
          aComp->actualVolume = (4*pow(2,0.5)*pow(VoxelRadius,3))*
            voxelCnt;
          setSystemSize(aComp->system, aComp->actualVolume*1e+3);
        }
    }
}

void SpatiocyteStepper::setSystemSize(System* aSystem, double aSize)
{
  Variable* aVariable(getVariable(aSystem, "SIZE"));
  if(aVariable)
    {
      aVariable->setValue(aSize);
    }
  else
    {
      createVariable(aSystem, "SIZE")->setValue(aSize);
    }
}

Variable* SpatiocyteStepper::createVariable(System* aSystem, String anID)
{
  String anEntityType("Variable");
  SystemPath aSystemPath(aSystem->getSystemPath());
  aSystemPath.push_back(aSystem->getID());
  FullID aFullID(anEntityType, aSystemPath, anID);
  Variable* aVariable(reinterpret_cast<Variable*>(
                              getModel()->createEntity("Variable", aFullID)));
  aVariable->setValue(0);
  return aVariable;
}

Species* SpatiocyteStepper::createSpecies(System* aSystem, String anID)
{
  Variable* aVariable(createVariable(aSystem, anID));
  return addSpecies(aVariable);
}

void SpatiocyteStepper::printSimulationParameters()
{
  std::cout << std::endl;
  std::cout << "   Voxel radius, r_v:" << VoxelRadius << " m" << std::endl;
  std::cout << "   Simulation height:" << theCenterPoint.y*2*VoxelRadius*2 <<
    " m" << std::endl;
  std::cout << "   Simulation width:" << theCenterPoint.z*2*VoxelRadius*2 << 
    " m" << std::endl;
  std::cout << "   Simulation length:" << theCenterPoint.x*2*VoxelRadius*2 <<
    " m" << std::endl;
  std::cout << "   Layer size:" << theTotalLayers << std::endl;
  std::cout << "   Row size:" << theTotalRows << std::endl;
  std::cout << "   Column size:" << theTotalCols << std::endl;
  std::cout << "   Total allocated voxels:" << 
    theTotalRows*theTotalLayers*theTotalCols << std::endl;
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      Comp* aComp(theComps[i]);
      double aSpecVolume(aComp->specVolume);
      double aSpecArea(aComp->specArea);
      double anActualVolume(aComp->actualVolume);
      double anActualArea(aComp->actualArea);
      switch(aComp->geometry)
        {
        case CUBOID:
          std::cout << "   Cuboid ";
          break;
        case ELLIPSOID:
          std::cout << "   Ellipsoid ";
          break;
        case CYLINDER:
          std::cout << "   Cylinder (radius=" << aComp->lengthY*VoxelRadius
            << "m, length=" << (aComp->lengthX)*VoxelRadius*2 << "m) ";
          break;
        case ROD:
          std::cout << "   Rod (radius=" << aComp->lengthY*VoxelRadius << 
            "m, cylinder length=" <<
            (aComp->lengthX-aComp->lengthY*2)*VoxelRadius*2 << "m) ";
          break;
        }
      std::cout << aComp->system->getFullID().asString();
      switch (aComp->dimension)
      { 
      case 1:
          std::cout << " Line compartment:" << std::endl;
          break;
      case 2:
          std::cout << " Surface compartment:" << std::endl;
          std::cout << "     [" << int(aSpecArea*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))/
                              (72*VoxelRadius*VoxelRadius)) << 
            "] Specified surface voxels {n_s = S_specified*"
            << "(6*2^0.5+4*3^0.5+3*6^0.5)/(72*r_v^2}" << std::endl;
          std::cout << "     [" << aComp->vacantSpecies->size() <<
            "] Actual surface voxels {n_s}" << std::endl;
          std::cout << "     [" << aSpecArea << " m^2] Specified surface area " <<
            "{S_specified}" << std::endl;
          std::cout << "     [" << anActualArea << " m^2] Actual surface area " <<
            "{S = (72*r_v^2)*n_s/(6*2^0.5+4*3^0.5+3*6^0.5)}" << std::endl;
          break;
      case 3:
      default:
          std::cout << " Volume compartment:" << std::endl;
          int voxelCnt(aComp->vacantSpecies->size());
          for(unsigned j(0); j != aComp->allSubs.size(); ++j)
            {
              //Don't include the comp surface voxels when you count the
              //total volume voxels:
              if(aComp->surfaceSub != aComp->allSubs[j])
                {
                  voxelCnt += aComp->allSubs[j]->vacantSpecies->size();
                }
            }
          std::cout << "     [" << int(aSpecVolume/(4*sqrt(2)*pow(VoxelRadius, 3))) << 
            "] Specified volume voxels {n_v = V_specified/(4*2^0.5*r_v^3)}" <<
          std::endl;  
          std::cout << "     [" << voxelCnt << "] Actual volume voxels {n_v}"  << std::endl;
          std::cout << "     [" << aSpecVolume << " m^3] Specified volume {V_specified}"
            << std::endl; 
          std::cout << "     [" << anActualVolume << " m^3] Actual volume " <<
            "{V = (4*2^0.5*r_v^3)*n_v}" << std::endl; 
      }
    }
  std::cout << std::endl;
}

void SpatiocyteStepper::readjustSurfaceBoundarySizes()
{
  Comp* aRootComp(theComps[0]);
  //[XY, XZ, YZ]PLANE: the boundary type of the surface when 
  //the geometry of the root Comp is CUBOID.
  //Where the root cuboid compartment is enclosed with a surface compartment,
  //the boundary type can be either REFLECTIVE, REMOVE_UPPER, REMOVE_LOWER or
  //REMOVE_BOTH. To make the actualVolume of the root compartment equivalent
  //to the specVolume we need to increase the size of the row, col and layer
  //according to the additional voxels required to occupy the surface voxels.
  if(aRootComp->surfaceSub)
    {
      if(aRootComp->xyPlane == REFLECTIVE)
        {
          theTotalRows += 2;
        }
      else if(aRootComp->xyPlane == REMOVE_UPPER ||
              aRootComp->xyPlane == REMOVE_LOWER)
        {
          theTotalRows += 1;
        }
      if(aRootComp->yzPlane == REFLECTIVE)
        {
          theTotalCols += 2;
        }
      else if(aRootComp->yzPlane == REMOVE_UPPER ||
              aRootComp->yzPlane == REMOVE_LOWER)
        {
          theTotalCols += 1;
        }
      if(aRootComp->xzPlane == REFLECTIVE)
        {
          theTotalLayers += 2;
        }
      else if(aRootComp->xzPlane == REMOVE_UPPER ||
              aRootComp->xzPlane == REMOVE_LOWER)
        {
          theTotalLayers += 1;
        }
    }
  else
    {
      //Boundary type can also be either PERIODIC or REFLECTIVE when there is
      //no surface compartment for the root compartment.
      //Increase the size of [row,layer,col] by one voxel and make them odd
      //sized if the system uses periodic boundary conditions.
      if(aRootComp->yzPlane == PERIODIC)
        { 
          if(theTotalCols%2 != 1)
            {
              theTotalCols += 1;
            }
          else
            {
              theTotalCols += 2;
            }
        }
      if(aRootComp->xzPlane == PERIODIC)
        {
          if(theTotalLayers%2 != 1)
            {
              theTotalLayers +=1;
            }
          else
            {
              theTotalLayers += 2;
            }
        }
      if(aRootComp->xyPlane == PERIODIC)
        {
          if(theTotalRows%2 != 1)
            {
              theTotalRows += 1;
            }
          else
            {
              theTotalRows += 2;
            }
        }
    }
  if(isPeriodicEdge)
    {
      if(theTotalCols%2 == 1)
        {
          theTotalCols += 1;
        }
      if(theTotalLayers%2 == 1)
        {
          theTotalLayers +=1;
        }
      if(theTotalRows%2 == 1)
        {
          theTotalRows += 1;
        }
    }
}

void SpatiocyteStepper::constructLattice()
{ 
  Comp* aRootComp(theComps[0]);
  const unsigned short rootID(aRootComp->vacantSpecies->getID());
  for(unsigned i(0); i != theTotalBoxSize; ++i)
    {
      const unsigned aSize(theRows[i]*theCols[i]*theLayers[i]);
      theIDs[i].resize(aSize);
      theInfo[i].resize(aSize);
      theAdjoins[i].resize(aSize*theAdjoinSize);
      theIDs[i][theNullMol] = theNullID;
    }
  for(unsigned i(0); i != theTotalBoxSize; ++i)
    {
      std::vector<unsigned short>& anIDs(theIDs[i]);
      std::vector<VoxelInfo>& anInfo(theInfo[i]);
      std::vector<unsigned>& anAdjoins(theAdjoins[i]);
      const unsigned bc(i/(theBoxRows*theBoxLayers)); 
      const unsigned bl((i%(theBoxRows*theBoxLayers))/theBoxRows); 
      const unsigned br((i%(theBoxRows*theBoxLayers))%theBoxRows); 
      const unsigned offset(i*theBoxMaxSize);
      unsigned aCol(0);
      unsigned aLayer(0);
      unsigned aRow(0);
      for(unsigned m(0); m != bc; ++m)
        {
          aCol += theCols[m];
        }
      for(unsigned m(0); m != br; ++m)
        {
          aRow += theRows[m];
        }
      for(unsigned m(0); m != bl; ++m)
        {
          aLayer += theLayers[m];
        }
      for(unsigned j(0); j != anIDs.size();  ++j)
        { 
          unsigned col(aCol+j/(theRows[i]*theLayers[i])); 
          unsigned layer(aLayer+(j%(theRows[i]*theLayers[i]))/theRows[i]); 
          unsigned row(aRow+(j%(theRows[i]*theLayers[i]))%theRows[i]); 
          const unsigned k(row+ 
                           theTotalRows*layer+ 
                           theTotalRows*theTotalLayers*col);
          /*
          anInfo[j].point.y = (aCol%2)*theHCPl+theHCPy*aLayer;
          anInfo[j].point.z = aRow*2*nVoxelRadius+((aLayer+aCol)%2)*nVoxelRadius;
          anInfo[j].point.x = aCol*theHCPx;
          */
          anInfo[j].diffuseSize = theAdjoinSize;
          anInfo[j].adjoinSize = theAdjoinSize;
          anInfo[j].coord = k;
          if(aRootComp->geometry == CUBOID || isInsideMol(k, aRootComp, 0))
            {
              //By default, the voxel is vacant and we set it to the root id:
              anIDs[j] = rootID;
              for(unsigned l(0); l != theAdjoinSize; ++l)
                { 
                  // By default let the adjoin voxel pointer point to the 
                  // source voxel (i.e., itself)
                  anAdjoins[j*theAdjoinSize+l] = j+offset;
                  //anAdjoins[j*theAdjoinSize+l] = j;
                }
            }
          else
            {
              //We set id = theNullID if it is an invalid voxel,
              //i.e., no molecules will occupy it:
              anIDs[j] = theNullID;
            }
        }
    }
  for(unsigned i(0); i != theTotalBoxSize; ++i)
    {
      for(unsigned j(0); j != theIDs[i].size();  ++j)
        { 
          concatenateVoxel(i, j);
        }
    }
}

void SpatiocyteStepper::constructLattice(unsigned anID)
{
  theThreads[anID]->runChildren();
  Comp* aRootComp(theComps[0]);
  const unsigned short rootID(aRootComp->vacantSpecies->getID());
  const unsigned startID(anID*theBoxSize);
  const unsigned endID(anID*theBoxSize+theBoxSize);
  for(unsigned i(startID); i != endID; ++i)
    {
      const unsigned aSize(theRows[i]*theCols[i]*theLayers[i]);
      theIDs[i].resize(aSize);
      theInfo[i].resize(aSize);
      theAdjoins[i].resize(aSize*theAdjoinSize);
      theIDs[i][theNullMol] = theNullID;
    }
  for(unsigned i(startID); i != endID; ++i)
    {
      std::vector<unsigned short>& anIDs(theIDs[i]);
      std::vector<VoxelInfo>& anInfo(theInfo[i]);
      std::vector<unsigned>& anAdjoins(theAdjoins[i]);
      const unsigned bc(i/(theBoxRows*theBoxLayers)); 
      const unsigned bl((i%(theBoxRows*theBoxLayers))/theBoxRows); 
      const unsigned br((i%(theBoxRows*theBoxLayers))%theBoxRows); 
      const unsigned offset(i*theBoxMaxSize);
      unsigned aCol(0);
      unsigned aLayer(0);
      unsigned aRow(0);
      for(unsigned m(0); m != bc; ++m)
        {
          aCol += theCols[m];
        }
      for(unsigned m(0); m != br; ++m)
        {
          aRow += theRows[m];
        }
      for(unsigned m(0); m != bl; ++m)
        {
          aLayer += theLayers[m];
        }
      for(unsigned j(0); j != anIDs.size();  ++j)
        { 
          unsigned col(aCol+j/(theRows[i]*theLayers[i])); 
          unsigned layer(aLayer+(j%(theRows[i]*theLayers[i]))/theRows[i]); 
          unsigned row(aRow+(j%(theRows[i]*theLayers[i]))%theRows[i]); 
          const unsigned k(row+ 
                           theTotalRows*layer+ 
                           theTotalRows*theTotalLayers*col);
          /*
          anInfo[j].point.y = (aCol%2)*theHCPl+theHCPy*aLayer;
          anInfo[j].point.z = aRow*2*nVoxelRadius+((aLayer+aCol)%2)*nVoxelRadius;
          anInfo[j].point.x = aCol*theHCPx;
          */
          anInfo[j].diffuseSize = theAdjoinSize;
          anInfo[j].adjoinSize = theAdjoinSize;
          anInfo[j].coord = k;
          if(aRootComp->geometry == CUBOID || isInsideMol(k, aRootComp, 0))
            {
              //By default, the voxel is vacant and we set it to the root id:
              anIDs[j] = rootID;
              for(unsigned l(0); l != theAdjoinSize; ++l)
                { 
                  // By default let the adjoin voxel pointer point to the 
                  // source voxel (i.e., itself)
                  anAdjoins[j*theAdjoinSize+l] = j+offset;
                  //anAdjoins[j*theAdjoinSize+l] = j;
                }
            }
          else
            {
              //We set id = theNullID if it is an invalid voxel,
              //i.e., no molecules will occupy it:
              anIDs[j] = theNullID;
            }
        }
    }
  theThreads[anID]->waitChildren();
}

void SpatiocyteStepper::concatenateLattice(unsigned anID)
{ 
  theThreads[anID]->runChildren();
  /*
  if(anID)
    {
      return;
    }
    */
  const unsigned startID(anID*theBoxSize);
  const unsigned endID(anID*theBoxSize+theBoxSize);
  for(unsigned i(startID); i != endID; ++i)
    {
      for(unsigned j(0); j != theIDs[i].size();  ++j)
        {
          concatenateVoxel(i, j);
        }
    }
  theThreads[anID]->waitChildren();
}

void SpatiocyteStepper::setBoundaries()
{
  //if the root comp is cuboid:
  if(theComps[0]->geometry == CUBOID)
    {
      concatenatePeriodicSurfaces();
    }
}

void SpatiocyteStepper::setPeriodicEdge()
{
  isPeriodicEdge = true;
}

bool SpatiocyteStepper::isPeriodicEdgeMol(unsigned aMol, Comp* aComp)
{
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aMol, aRow, aLayer, aCol);
  if(aComp->system->getSuperSystem()->isRootSystem() &&
     ((aRow <= 1 && aCol <= 1) ||
      (aRow <= 1 && aLayer <= 1) ||
      (aCol <= 1 && aLayer <= 1) ||
      (aRow <= 1 && aCol >= theTotalCols-2) ||
      (aRow <= 1 && aLayer >= theTotalLayers-2) ||
      (aCol <= 1 && aRow >= theTotalRows-2) ||
      (aCol <= 1 && aLayer >= theTotalLayers-2) ||
      (aLayer <= 1 && aRow >= theTotalRows-2) ||
      (aLayer <= 1 && aCol >= theTotalCols-2) ||
      (aRow >= theTotalRows-2 && aCol >= theTotalCols-2) ||
      (aRow >= theTotalRows-2 && aLayer >= theTotalLayers-2) ||
      (aCol >= theTotalCols-2 && aLayer >= theTotalLayers-2)))
    {
      return true;
    }
  return false;
}

bool SpatiocyteStepper::isRemovableEdgeMol(unsigned aMol, Comp* aComp)
{
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aMol, aRow, aLayer, aCol);
  int sharedCnt(0);
  int removeCnt(0);
  //Minus 1 to maxRow to account for surfaces that use two rows to envelope the
  //volume:
  if(aRow >= aComp->maxRow-1)
    {
      ++sharedCnt;
      if(aComp->xyPlane == REMOVE_UPPER || aComp->xyPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  //Add 1 to minRow to account for surfaces that use two rows to envelope the
  //volume:
  if(aRow <= aComp->minRow+1)
    {
      ++sharedCnt;
      if(aComp->xyPlane == REMOVE_LOWER || aComp->xyPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aLayer >= aComp->maxLayer-1)
    {
      ++sharedCnt;
      if(aComp->xzPlane == REMOVE_UPPER || aComp->xzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aLayer <= aComp->minLayer+1)
    {
      ++sharedCnt;
      if(aComp->xzPlane == REMOVE_LOWER || aComp->xzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aCol >= aComp->maxCol)
    {
      ++sharedCnt;
      if(aComp->yzPlane == REMOVE_UPPER || aComp->yzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aCol <= aComp->minCol) 
    {
      ++sharedCnt;
      if(aComp->yzPlane == REMOVE_LOWER || aComp->yzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(!removeCnt)
    {
      return false;
    }
  else
    {
      return sharedCnt == removeCnt;
    }
}

void SpatiocyteStepper::concatenateVoxel(const unsigned aBox,
                                         const unsigned aMol)
{
  const unsigned aCol(aMol/(theRows[aBox]*theLayers[aBox])); 
  const unsigned aLayer((aMol%(theRows[aBox]*theLayers[aBox]))/theRows[aBox]); 
  const unsigned aRow((aMol%(theRows[aBox]*theLayers[aBox]))%theRows[aBox]); 
  const unsigned bCol(aBox/(theBoxRows*theBoxLayers)); 
  const unsigned bLayer((aBox%(theBoxRows*theBoxLayers))/theBoxRows); 
  const unsigned bRow((aBox%(theBoxRows*theBoxLayers))%theBoxRows); 
  if(aRow > 0)
    { 
      concatenateRows(aBox, aBox, aMol, aRow-1, aLayer, aCol);
    } 
  else if(bRow > 0)
    {
      const unsigned tar(bRow-1+
                         theBoxRows*bLayer+
                         theBoxRows*theBoxLayers*bCol);
      concatenateRows(aBox, tar, aMol, theRows[tar]-1, aLayer, aCol);
    }
  if(aLayer > 0)
    {
      concatenateLayers(aBox, aBox, aMol, aRow, aLayer-1, aCol); 
    }
  else if(bLayer > 0)
    {
      const unsigned tar(bRow+
                         theBoxRows*(bLayer-1)+
                         theBoxRows*theBoxLayers*bCol);
      concatenateLayers(aBox, tar, aMol, aRow, theLayers[tar]-1, aCol); 
    }
  if(aCol > 0)
    { 
      concatenateCols(aBox, aBox, aMol, aRow, aLayer, aCol-1); 
    }
  else if(bCol > 0)
    {
      const unsigned tar(bRow+
                         theBoxRows*bLayer+
                         theBoxRows*theBoxLayers*(bCol-1));
      concatenateCols(aBox, tar, aMol, aRow, aLayer, theCols[tar]-1); 
    }
}

void SpatiocyteStepper::concatenateRows(const unsigned src,
                                        const unsigned tar,
                                        const unsigned a,
                                        const unsigned aRow,
                                        const unsigned aLayer,
                                        const unsigned aCol)
{
  const unsigned b(aRow+ 
                   theRows[tar]*aLayer+ 
                   theRows[tar]*theLayers[tar]*aCol);
  theAdjoins[src][a*theAdjoinSize+NORTH] = b+tar*theBoxMaxSize;
  theAdjoins[tar][b*theAdjoinSize+SOUTH] = a+src*theBoxMaxSize;
}

void SpatiocyteStepper::concatenateLayers(const unsigned src,
                                          const unsigned tar,
                                          const unsigned a,
                                          const unsigned aRow,
                                          const unsigned aLayer,
                                          const unsigned aCol)
{
  const unsigned b(aRow+
                   theRows[tar]*aLayer+
                   theRows[tar]*theLayers[tar]*aCol);
  const unsigned bCol(tar/(theBoxRows*theBoxLayers)); 
  const unsigned bLayer((tar%(theBoxRows*theBoxLayers))/theBoxRows); 
  const unsigned bRow((tar%(theBoxRows*theBoxLayers))%theBoxRows); 
  if((aLayer+1)%2+(aCol)%2 == 1)
    {
      theAdjoins[src][a*theAdjoinSize+VENTRALN] = b+tar*theBoxMaxSize;
      theAdjoins[tar][b*theAdjoinSize+DORSALS] = a+src*theBoxMaxSize;
      if(aRow < theRows[tar]-1)
        {
          const unsigned c(aRow+1+ 
                           theRows[tar]*aLayer+
                           theRows[tar]*theLayers[tar]*aCol);
          theAdjoins[src][a*theAdjoinSize+VENTRALS] = c+tar*theBoxMaxSize;
          theAdjoins[tar][c*theAdjoinSize+DORSALN] = a+src*theBoxMaxSize;
        }
      else if(bRow < theBoxRows-1)
        {
          const unsigned tarc(bRow+1+
                              theBoxRows*bLayer+
                              theBoxRows*theBoxLayers*bCol);
          const unsigned c(0+ 
                           theRows[tarc]*aLayer+
                           theRows[tarc]*theLayers[tarc]*aCol);
          theAdjoins[src][a*theAdjoinSize+VENTRALS] = c+tarc*theBoxMaxSize;
          theAdjoins[tarc][c*theAdjoinSize+DORSALN] = a+src*theBoxMaxSize;
        }
    }
  else
    {
      theAdjoins[src][a*theAdjoinSize+VENTRALS] = b+tar*theBoxMaxSize;
      theAdjoins[tar][b*theAdjoinSize+DORSALN] = a+src*theBoxMaxSize;
      if(aRow > 0)
        {
          const unsigned c(aRow-1+ 
                           theRows[tar]*aLayer+
                           theRows[tar]*theLayers[tar]*aCol);
          theAdjoins[src][a*theAdjoinSize+VENTRALN] = c+tar*theBoxMaxSize;
          theAdjoins[tar][c*theAdjoinSize+DORSALS] = a+src*theBoxMaxSize;
        }
      else if(bRow > 0)
        {
          const unsigned tarc(bRow-1+
                              theBoxRows*bLayer+
                              theBoxRows*theBoxLayers*bCol);
          const unsigned c(theRows[tarc]-1+ 
                           theRows[tarc]*aLayer+
                           theRows[tarc]*theLayers[tarc]*aCol);
          theAdjoins[src][a*theAdjoinSize+VENTRALN] = c+tarc*theBoxMaxSize;
          theAdjoins[tarc][c*theAdjoinSize+DORSALS] = a+src*theBoxMaxSize;
        }
    }
}


void SpatiocyteStepper::concatenateCols(const unsigned src,
                                        const unsigned tar,
                                        const unsigned a,
                                        const unsigned aRow,
                                        const unsigned aLayer,
                                        const unsigned aCol)
{
  const unsigned b(aRow+
                   theRows[tar]*aLayer+
                   theRows[tar]*theLayers[tar]*aCol);
  const unsigned bCol(tar/(theBoxRows*theBoxLayers)); 
  const unsigned bLayer((tar%(theBoxRows*theBoxLayers))/theBoxRows); 
  const unsigned bRow((tar%(theBoxRows*theBoxLayers))%theBoxRows); 
  if(aLayer%2 == 0)
    {
      if((aCol+1)%2 == 1)
        {
          theAdjoins[src][a*theAdjoinSize+NW] = b+tar*theBoxMaxSize;
          theAdjoins[tar][b*theAdjoinSize+SE] = a+src*theBoxMaxSize;
          if(aRow < theRows[tar] - 1)
            {
              const unsigned c(aRow+1+ 
                               theRows[tar]*aLayer+
                               theRows[tar]*theLayers[tar]*aCol);
              theAdjoins[src][a*theAdjoinSize+SW] = c+tar*theBoxMaxSize;
              theAdjoins[tar][c*theAdjoinSize+NE] = a+src*theBoxMaxSize;
            }
          else if(bRow < theBoxRows-1)
            {
              const unsigned tarc(bRow+1+
                                  theBoxRows*bLayer+
                                  theBoxRows*theBoxLayers*bCol);
              const unsigned c(0+ 
                               theRows[tarc]*aLayer+
                               theRows[tarc]*theLayers[tarc]*aCol);
              theAdjoins[src][a*theAdjoinSize+SW] = c+tarc*theBoxMaxSize;
              theAdjoins[tarc][c*theAdjoinSize+NE] = a+src*theBoxMaxSize;
            }
          if(aLayer < theLayers[tar]-1)
            {
              const unsigned c(aRow+ 
                               theRows[tar]*(aLayer+1)+
                               theRows[tar]*theLayers[tar]*aCol);
              theAdjoins[src][a*theAdjoinSize+WEST] = c+tar*theBoxMaxSize;
              theAdjoins[tar][c*theAdjoinSize+EAST] = a+src*theBoxMaxSize;
            }
          else if(bLayer < theBoxLayers-1)
            {
              const unsigned tarc(bRow+
                                  theBoxRows*(bLayer+1)+
                                  theBoxRows*theBoxLayers*bCol);
              const unsigned c(aRow+
                               0+
                               theRows[tarc]*theLayers[tarc]*aCol);
              theAdjoins[src][a*theAdjoinSize+WEST] = c+tarc*theBoxMaxSize;
              theAdjoins[tarc][c*theAdjoinSize+EAST] = a+src*theBoxMaxSize;
            }
        }
      else
        {
          theAdjoins[src][a*theAdjoinSize+SW] = b+tar*theBoxMaxSize;
          theAdjoins[tar][b*theAdjoinSize+NE] = a+src*theBoxMaxSize;
          if(aRow > 0)
            {
              const unsigned c(aRow-1+ 
                               theRows[tar]*aLayer+
                               theRows[tar]*theLayers[tar]*aCol);
              theAdjoins[src][a*theAdjoinSize+NW] = c+tar*theBoxMaxSize;
              theAdjoins[tar][c*theAdjoinSize+SE] = a+src*theBoxMaxSize;
            }
          else if(bRow > 0)
            {
              const unsigned tarc(bRow-1+
                                  theBoxRows*bLayer+
                                  theBoxRows*theBoxLayers*bCol);
              const unsigned c(theRows[tarc]-1+ 
                               theRows[tarc]*aLayer+
                               theRows[tarc]*theLayers[tarc]*aCol);
              theAdjoins[src][a*theAdjoinSize+NW] = c+tarc*theBoxMaxSize;
              theAdjoins[tarc][c*theAdjoinSize+SE] = a+src*theBoxMaxSize;
            }
          if(aLayer > 0)
            {
              const unsigned c(aRow+ 
                               theRows[tar]*(aLayer-1)+
                               theRows[tar]*theLayers[tar]*aCol);
              theAdjoins[src][a*theAdjoinSize+WEST] = c+tar*theBoxMaxSize;
              theAdjoins[tar][c*theAdjoinSize+EAST] = a+src*theBoxMaxSize;
            }
          else if(bLayer > 0)
            {
              const unsigned tarc(bRow+
                                  theBoxRows*(bLayer-1)+
                                  theBoxRows*theBoxLayers*bCol);
              const unsigned c(aRow+ 
                               theRows[tarc]*(theLayers[tarc]-1)+
                               theRows[tarc]*theLayers[tarc]*aCol);
              theAdjoins[src][a*theAdjoinSize+WEST] = c+tarc*theBoxMaxSize;
              theAdjoins[tarc][c*theAdjoinSize+EAST] = a+src*theBoxMaxSize;
            }
        }
    }
  else
    {
      if((aCol+1)%2 == 1)
        {
          theAdjoins[src][a*theAdjoinSize+SW] = b+tar*theBoxMaxSize;
          theAdjoins[tar][b*theAdjoinSize+NE] = a+src*theBoxMaxSize;
          if(aRow > 0)
            {
              const unsigned c(aRow-1+ 
                               theRows[tar]*aLayer+
                               theRows[tar]*theLayers[tar]*aCol);
              theAdjoins[src][a*theAdjoinSize+NW] = c+tar*theBoxMaxSize;
              theAdjoins[tar][c*theAdjoinSize+SE] = a+src*theBoxMaxSize;
            }
          else if(bRow > 0)
            {
              const unsigned tarc(bRow-1+
                                  theBoxRows*bLayer+
                                  theBoxRows*theBoxLayers*bCol);
              const unsigned c(theRows[tarc]-1+ 
                               theRows[tarc]*aLayer+
                               theRows[tarc]*theLayers[tarc]*aCol);
              theAdjoins[src][a*theAdjoinSize+NW] = c+tarc*theBoxMaxSize;
              theAdjoins[tarc][c*theAdjoinSize+SE] = a+src*theBoxMaxSize;
            }
          if(aLayer < theLayers[tar]-1)
            {
              const unsigned c(aRow+ 
                               theRows[tar]*(aLayer+1)+
                               theRows[tar]*theLayers[tar]*aCol);
              theAdjoins[src][a*theAdjoinSize+WEST] = c+tar*theBoxMaxSize;
              theAdjoins[tar][c*theAdjoinSize+EAST] = a+src*theBoxMaxSize;
            }
          else if(bLayer < theBoxLayers-1)
            {
              const unsigned tarc(bRow+
                                  theBoxRows*(bLayer+1)+
                                  theBoxRows*theBoxLayers*bCol);
              const unsigned c(aRow+
                               0+
                               theRows[tarc]*theLayers[tarc]*aCol);
              theAdjoins[src][a*theAdjoinSize+WEST] = c+tarc*theBoxMaxSize;
              theAdjoins[tarc][c*theAdjoinSize+EAST] = a+src*theBoxMaxSize;
            }
        }
      else
        {
          theAdjoins[src][a*theAdjoinSize+NW] = b+tar*theBoxMaxSize;
          theAdjoins[tar][b*theAdjoinSize+SE] = a+src*theBoxMaxSize;
          if(aRow < theRows[tar] - 1)
            {
              const unsigned c(aRow+1+ 
                               theRows[tar]*aLayer+
                               theRows[tar]*theLayers[tar]*aCol);
              theAdjoins[src][a*theAdjoinSize+SW] = c+tar*theBoxMaxSize;
              theAdjoins[tar][c*theAdjoinSize+NE] = a+src*theBoxMaxSize;
            }
          else if(bRow < theBoxRows-1)
            {
              const unsigned tarc(bRow+1+
                                  theBoxRows*bLayer+
                                  theBoxRows*theBoxLayers*bCol);
              const unsigned c(0+ 
                               theRows[tarc]*aLayer+
                               theRows[tarc]*theLayers[tarc]*aCol);
              theAdjoins[src][a*theAdjoinSize+SW] = c+tarc*theBoxMaxSize;
              theAdjoins[tarc][c*theAdjoinSize+NE] = a+src*theBoxMaxSize;
            }
          if(aLayer > 0)
            {
              const unsigned c(aRow+ 
                               theRows[tar]*(aLayer-1)+
                               theRows[tar]*theLayers[tar]*aCol);
              theAdjoins[src][a*theAdjoinSize+WEST] = c+tar*theBoxMaxSize;
              theAdjoins[tar][c*theAdjoinSize+EAST] = a+src*theBoxMaxSize;
            }
          else if(bLayer > 0)
            {
              const unsigned tarc(bRow+
                                  theBoxRows*(bLayer-1)+
                                  theBoxRows*theBoxLayers*bCol);
              const unsigned c(aRow+ 
                               theRows[tarc]*(theLayers[tarc]-1)+
                               theRows[tarc]*theLayers[tarc]*aCol);
              theAdjoins[src][a*theAdjoinSize+WEST] = c+tarc*theBoxMaxSize;
              theAdjoins[tarc][c*theAdjoinSize+EAST] = a+src*theBoxMaxSize;
            }
        }
    }
}

Variable* SpatiocyteStepper::getVariable(System* aSystem, String const& anID)
{
  FOR_ALL(System::Variables, aSystem->getVariables())
    {
      Variable* aVariable(i->second);
      if(aVariable->getID() == anID)
        {
          return aVariable;
        }
    }
  return NULL;
}

void SpatiocyteStepper::setCompProperties(Comp* aComp)
{
  System* aSystem(aComp->system);
  double aRadius(0);
  switch(aComp->geometry)
    {
    case CUBOID:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHX of the Cuboid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHY of the Cuboid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aComp->lengthZ)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHZ of the Cuboid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      aComp->specVolume = aComp->lengthX*aComp->lengthY*
        aComp->lengthZ;
      aComp->specArea = getCuboidSpecArea(aComp);
      break;
    case ELLIPSOID:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHX of the Ellipsoid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHY of the Ellipsoid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      if(!aComp->lengthZ)
        {
          THROW_EXCEPTION(NotFound,
                          "Property LENGTHZ of the Ellipsoid Comp " +
                          aSystem->getFullID().asString() + " is not defined.");
        }
      aComp->specVolume = 4*M_PI*aComp->lengthX*
        aComp->lengthY* aComp->lengthZ/24;
      aComp->specArea = 4*M_PI*
        pow((pow(aComp->lengthX/2, 1.6075)*
             pow(aComp->lengthY/2, 1.6075)+
             pow(aComp->lengthX/2, 1.6075)*
             pow(aComp->lengthZ/2, 1.6075)+
             pow(aComp->lengthY/2, 1.6075)*
             pow(aComp->lengthZ/2, 1.6075))/ 3, 1/1.6075); 
      break;
    case CYLINDER:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHX of the Cylinder Comp "
                          + aSystem->getFullID().asString() + " not defined." );
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHY of the Cylinder Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies its diameter, not defined." );
        }
      aRadius = aComp->lengthY/2;
      aComp->specVolume = (aComp->lengthX)*(M_PI*aRadius*aRadius);
      aComp->lengthZ = aComp->lengthY;
      aComp->specArea = 2*M_PI*aRadius*(aComp->lengthX);
      break;
    case ROD:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHX of the Rod Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the rod length (cylinder length + "
                          + "2*hemisphere radius), not defined." );
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHY of the Rod Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies its diameter, not defined." );
        }
      aRadius = aComp->lengthY/2;
      aComp->specVolume = (aComp->lengthX+(4*aRadius/3)-(2*aRadius))*
        (M_PI*aRadius*aRadius);
      aComp->lengthZ = aComp->lengthY;
      aComp->specArea = 4*M_PI*aRadius*aRadius+
        2*M_PI*aRadius*(aComp->lengthX-2*aRadius);
      break;
    case PYRAMID:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHX of the Pyramid Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the length of the pyramid, "
                          + "not defined." );
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHY of the Pyramid Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the height of the pyramid, "
                          + "not defined." );
        }
      if(!aComp->lengthZ)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHZ of the Pyramid Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the depth of the pyramid, "
                          + "not defined." );
        }
      aComp->specVolume = aComp->lengthX*aComp->lengthY*aComp->lengthZ/3;
      aRadius = sqrt(pow(aComp->lengthX/2,2)+pow(aComp->lengthY,2));
      aComp->specArea =  aComp->lengthZ*aComp->lengthX+
        2*(aComp->lengthZ*aRadius)+2*(aComp->lengthX*aRadius);
      break;
    case ERYTHROCYTE:
      if(!aComp->lengthX)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHX of the Erythrocyte Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the length of the Erythrocyte, "
                          + "not defined." );
        }
      if(!aComp->lengthY)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHY of the Erythrocyte Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the height of the Erythrocyte, "
                          + "not defined." );
        }
      if(!aComp->lengthZ)
        {
          THROW_EXCEPTION(NotFound, "Property LENGTHZ of the Erythrocyte Comp "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the depth of the Erythrocyte, "
                          + "not defined." );
        }
      aComp->specVolume = 4*M_PI*aComp->lengthX*
        aComp->lengthY* aComp->lengthZ/24;
      aComp->specArea = 4*M_PI*
        pow((pow(aComp->lengthX/2, 1.6075)*
             pow(aComp->lengthY/2, 1.6075)+
             pow(aComp->lengthX/2, 1.6075)*
             pow(aComp->lengthZ/2, 1.6075)+
             pow(aComp->lengthY/2, 1.6075)*
             pow(aComp->lengthZ/2, 1.6075))/ 3, 1/1.6075); 
      break;
    }
  aComp->lengthX /= VoxelRadius*2;
  aComp->lengthY /= VoxelRadius*2;
  aComp->lengthZ /= VoxelRadius*2;
}

void SpatiocyteStepper::setCompCenterPoint(Comp* aComp)
{
  System* aSystem(aComp->system);
  System* aSuperSystem(aSystem->getSuperSystem());
  if(aComp->dimension == 2)
    {
      aSystem = aComp->system->getSuperSystem();
      aSuperSystem = aSystem->getSuperSystem();
    }
  else if(aComp->dimension == 1)
    {
      aSystem = aComp->system->getSuperSystem()->getSuperSystem();
      aSuperSystem = aSystem->getSuperSystem();
    }
  Comp* aSuperComp(system2Comp(aSuperSystem));
  //The center with reference to the immediate super system:
  aComp->centerPoint = aSuperComp->centerPoint;
  aComp->centerPoint.x += aComp->originX*aSuperComp->lengthX/2;
  aComp->centerPoint.y += aComp->originY*aSuperComp->lengthY/2;
  aComp->centerPoint.z += aComp->originZ*aSuperComp->lengthZ/2;
}

double SpatiocyteStepper::getCuboidSpecArea(Comp* aComp)
{
  double anArea(0);
  if(aComp->xyPlane == UNIPERIODIC || 
     aComp->xyPlane == REMOVE_UPPER ||
     aComp->xyPlane == REMOVE_LOWER)
    { 
      anArea += aComp->lengthX*aComp->lengthY; 
    }
  else if(aComp->xyPlane == REFLECTIVE)
    {
      anArea += 2*aComp->lengthX*aComp->lengthY; 
    }
  if(aComp->xzPlane == UNIPERIODIC || 
     aComp->xzPlane == REMOVE_UPPER ||
     aComp->xzPlane == REMOVE_LOWER)
    { 
      anArea += aComp->lengthX*aComp->lengthZ; 
    }
  else if(aComp->xzPlane == REFLECTIVE)
    {
      anArea += 2*aComp->lengthX*aComp->lengthZ; 
    }
  if(aComp->yzPlane == UNIPERIODIC || 
     aComp->yzPlane == REMOVE_UPPER ||
     aComp->yzPlane == REMOVE_LOWER)
    { 
      anArea += aComp->lengthY*aComp->lengthZ; 
    }
  else if(aComp->yzPlane == REFLECTIVE)
    {
      anArea += 2*aComp->lengthY*aComp->lengthZ; 
    }
  return anArea;
}

Point SpatiocyteStepper::coord2point(unsigned aMol)
{
  unsigned aGlobalCol;
  unsigned aGlobalLayer;
  unsigned aGlobalRow;
  coord2global(aMol, aGlobalRow, aGlobalLayer, aGlobalCol);
  //the center point of a voxel 
  Point aPoint;
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      aPoint.y = (aGlobalCol%2)*theHCPl+theHCPy*aGlobalLayer;
      aPoint.z = aGlobalRow*2*nVoxelRadius+
        ((aGlobalLayer+aGlobalCol)%2)*nVoxelRadius;
      aPoint.x = aGlobalCol*theHCPx;
      break;
    case CUBIC_LATTICE:
      aPoint.y = aGlobalLayer*2*nVoxelRadius;
      aPoint.z = aGlobalRow*2*nVoxelRadius;
      aPoint.x = aGlobalCol*2*nVoxelRadius;
      break;
    }
  return aPoint;
};

double SpatiocyteStepper::getRowLength()
{
  return nVoxelRadius*2;
}

double SpatiocyteStepper::getColLength()
{
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      return theHCPx;
    case CUBIC_LATTICE:
      return nVoxelRadius*2;
    }
  return nVoxelRadius*2;
}

double SpatiocyteStepper::getLayerLength()
{
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      return theHCPy;
    case CUBIC_LATTICE:
      return nVoxelRadius*2;
    }
  return nVoxelRadius*2;
}

double SpatiocyteStepper::getMinLatticeSpace()
{
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      return theHCPl;
    case CUBIC_LATTICE:
      return nVoxelRadius*2;
    }
  return nVoxelRadius*2;
}

unsigned SpatiocyteStepper::point2coord(Point& aPoint)
{
  unsigned aGlobalRow(0);
  unsigned aGlobalLayer(0);
  unsigned aGlobalCol(0);
  point2global(aPoint, aGlobalRow, aGlobalLayer, aGlobalCol);
  return global2coord(aGlobalRow, aGlobalLayer, aGlobalCol);
}

unsigned SpatiocyteStepper::global2coord(unsigned aGlobalRow,
                                             unsigned aGlobalLayer,
                                             unsigned aGlobalCol)
{
  return aGlobalRow+theTotalRows*aGlobalLayer+theTotalRows*theTotalLayers*aGlobalCol;
}

void SpatiocyteStepper::point2global(Point aPoint, 
                                     unsigned& aGlobalRow,
                                     unsigned& aGlobalLayer,
                                     unsigned& aGlobalCol)
{
  double row(0);
  double layer(0);
  double col(0);
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      col = rint(aPoint.x/theHCPx);
      layer = rint((aPoint.y-(aGlobalCol%2)*theHCPl)/
                                        theHCPy);
      row = rint((aPoint.z-((aGlobalLayer+aGlobalCol)%2)*
          nVoxelRadius)/(2*nVoxelRadius));
      break;
    case CUBIC_LATTICE:
      col = rint(aPoint.x/(2*nVoxelRadius));
      layer = rint(aPoint.y/(2*nVoxelRadius));
      row = rint(aPoint.z/(2*nVoxelRadius));
      break;
    }
  if(row < 0)
    {
      row = 0;
    }
  if(layer < 0)
    {
      layer = 0;
    }
  if(col < 0)
    {
      col = 0;
    }
  aGlobalRow = (unsigned)row;
  aGlobalLayer = (unsigned)layer;
  aGlobalCol = (unsigned)col;
  if(aGlobalCol >= theTotalCols)
    {
      aGlobalCol = theTotalCols-1;
    }
  if(aGlobalRow >= theTotalRows)
    {
      aGlobalRow = theTotalRows-1;
    }
  if(aGlobalLayer >= theTotalLayers)
    {
      aGlobalLayer = theTotalLayers-1;
    }
}

void SpatiocyteStepper::coord2global(unsigned aMol,
                                     unsigned& aGlobalRow,
                                     unsigned& aGlobalLayer,
                                     unsigned& aGlobalCol) 
{
  aGlobalCol = (aMol)/(theTotalRows*theTotalLayers);
  aGlobalLayer = ((aMol)%(theTotalRows*theTotalLayers))/theTotalRows;
  aGlobalRow = ((aMol)%(theTotalRows*theTotalLayers))%theTotalRows;
}





void SpatiocyteStepper::concatenatePeriodicSurfaces()
{
  Comp* aRootComp(theComps[0]);
  //concatenate periodic rows
  for(unsigned i(0); i != theBoxLayers; ++i)
    {
      for(unsigned j(0); j != theBoxCols; ++j)
        {
          const unsigned boxA(0+
                              theBoxRows*i+
                              theBoxRows*theBoxLayers*j);
          const unsigned boxB(theBoxRows-1+
                              theBoxRows*i+
                              theBoxRows*theBoxLayers*j);
          for(unsigned k(0); k != theLayers[boxA]; ++k)
            {
              for(unsigned l(0); l != theCols[boxA]; ++l)
                { 
                  unsigned A(0+ 
                             theRows[boxA]*k+ 
                             theRows[boxA]*theLayers[boxA]*l);
                  unsigned B(theRows[boxB]-1+ 
                             theRows[boxB]*k+ 
                             theRows[boxB]*theLayers[boxB]*l);
                  if(aRootComp->xyPlane == UNIPERIODIC)
                    { 
                      replaceUniVoxel(boxA, boxB, A, B);
                    }
                  else if(aRootComp->xyPlane == PERIODIC)
                    { 
                      replaceVoxel(boxA, boxB, A, B);
                    }
                  else if(!isPeriodicEdge)
                    {
                      //We cannot have valid voxels pointing to itself it it
                      //is not periodic o avoid incorrect homodimerization
                      //reaction. So we set such molecules to null ID.
                      theIDs[boxA][A] = theNullID;
                      theIDs[boxB][B] = theNullID;
                    }
                }
            }
        }
    }
  //concatenate periodic layers
  for(unsigned i(0); i != theBoxRows; ++i)
    {
      for(unsigned j(0); j != theBoxCols; ++j)
        {
          const unsigned boxA(i+
                              theBoxRows*0+
                              theBoxRows*theBoxLayers*j);
          const unsigned boxB(i+
                              theBoxRows*(theBoxLayers-1)+
                              theBoxRows*theBoxLayers*j);
          for(unsigned k(0); k != theRows[boxA]; ++k)
            {
              for(unsigned l(0); l != theCols[boxA]; ++l)
                { 
                  unsigned A(k+ 
                             theRows[boxA]*0+ 
                             theRows[boxA]*theLayers[boxA]*l);
                  unsigned B(k+ 
                             theRows[boxB]*(theLayers[boxB]-1)+ 
                             theRows[boxB]*theLayers[boxB]*l);
                  if(aRootComp->xzPlane == UNIPERIODIC)
                    { 
                      replaceUniVoxel(boxA, boxB, A, B);
                    }
                  else if(aRootComp->xzPlane == PERIODIC)
                    { 
                      replaceVoxel(boxA, boxB, A, B);
                    }
                  else if(!isPeriodicEdge)
                    {
                      //We cannot have valid voxels pointing to itself it it
                      //is not periodic o avoid incorrect homodimerization
                      //reaction. So we set such molecules to null ID.
                      theIDs[boxA][A] = theNullID;
                      theIDs[boxB][B] = theNullID;
                    }
                }
            }
        }
    }
  //concatenate periodic cols
  for(unsigned i(0); i != theBoxRows; ++i)
    {
      for(unsigned j(0); j != theBoxLayers; ++j)
        {
          const unsigned boxA(i+
                              theBoxRows*j+
                              theBoxRows*theBoxLayers*0);
          const unsigned boxB(i+
                              theBoxRows*j+
                              theBoxRows*theBoxLayers*(theBoxCols-1));
          for(unsigned k(0); k != theRows[boxA]; ++k)
            {
              for(unsigned l(0); l != theLayers[boxA]; ++l)
                { 
                  unsigned A(k+ 
                             theRows[boxA]*l+ 
                             theRows[boxA]*theLayers[boxA]*0);
                  unsigned B(k+ 
                             theRows[boxB]*l+ 
                             theRows[boxB]*theLayers[boxB]*(theCols[boxB]-1));
                  if(aRootComp->yzPlane == UNIPERIODIC)
                    { 
                      replaceUniVoxel(boxA, boxB, A, B);
                    }
                  else if(aRootComp->yzPlane == PERIODIC)
                    { 
                      replaceVoxel(boxA, boxB, A, B);
                    }
                  else if(!isPeriodicEdge)
                    {
                      //We cannot have valid voxels pointing to itself it it
                      //is not periodic o avoid incorrect homodimerization
                      //reaction. So we set such molecules to null ID.
                      theIDs[boxA][A] = theNullID;
                      theIDs[boxB][B] = theNullID;
                    }
                }
            }
        }
    }
}

unsigned SpatiocyteStepper::coord2row(unsigned aMol)
{
  return (aMol%(theTotalRows*theTotalLayers))%theTotalRows;
}

unsigned SpatiocyteStepper::coord2layer(unsigned aMol)
{
  return (aMol%(theTotalRows*theTotalLayers))/theTotalRows;
}

unsigned SpatiocyteStepper::coord2col(unsigned aMol)
{
  return aMol/(theTotalRows*theTotalLayers);
}

void SpatiocyteStepper::replaceVoxel(const unsigned boxA,
                                     const unsigned boxB,
                                     const unsigned A,
                                     const unsigned B)
{
  const unsigned short& idA(theIDs[boxA][A]);
  unsigned short& idB(theIDs[boxB][B]);
  const unsigned absA(A+boxA*theBoxMaxSize);
  const unsigned absB(B+boxB*theBoxMaxSize);
  if(idA != theNullID && idB != theNullID)
    {
      for(unsigned j(0); j != theAdjoinSize; ++j)
        {
          if(theAdjoins[boxA][A*theAdjoinSize+j] == absA &&
             theAdjoins[boxB][B*theAdjoinSize+j] != absB)
            {
              theAdjoins[boxA][A*theAdjoinSize+j] = 
                theAdjoins[boxB][B*theAdjoinSize+j];
              for(unsigned k(0); k != theAdjoinSize; ++k)
                {
                  const unsigned box(theAdjoins[boxB][B*theAdjoinSize+j]/
                                     theBoxMaxSize);
                  const unsigned adj(theAdjoins[boxB][B*theAdjoinSize+j]%
                                      theBoxMaxSize);
                  if(theAdjoins[box][adj*theAdjoinSize+k] == absB)
                    {
                      theAdjoins[box][adj*theAdjoinSize+k] = absA;
                    }
                }
            }
        }
      idB = theNullID;
    }
}

void SpatiocyteStepper::replaceUniVoxel(const unsigned boxA,
                                        const unsigned boxB,
                                        const unsigned A,
                                        const unsigned B)
{
  const unsigned short& idA(theIDs[boxA][A]);
  unsigned short& idB(theIDs[boxB][B]);
  const unsigned absA(A+boxA*theBoxMaxSize);
  const unsigned absB(B+boxB*theBoxMaxSize);
  if(idA != theNullID && idB != theNullID)
    {
      for(unsigned j(0); j != theAdjoinSize; ++j)
        {
          if(theAdjoins[boxA][A*theAdjoinSize+j] == absA)
            {
              for(unsigned k(0); k != theAdjoinSize; ++k)
                {
                  const unsigned box(theAdjoins[boxB][B*theAdjoinSize+j]/
                                     theBoxMaxSize);
                  const unsigned adj(theAdjoins[boxB][B*theAdjoinSize+j]%
                                      theBoxMaxSize);
                  if(theAdjoins[box][adj*theAdjoinSize+k] == absB)
                    {
                      theAdjoins[box][adj*theAdjoinSize+k] = absA;
                    }
                }
            }
        }
      idB = theNullID;
    }
}

void SpatiocyteStepper::shuffleAdjoins()
{
  /*
  for(unsigned i(0); i != theIDs.size(); ++i)
    {
      if(theIDs[i] != theNullID)
        { 
          gsl_ran_shuffle(getRng(), &theAdjoins[i*theAdjoinSize], theAdjoinSize,
                          sizeof(unsigned));
        }
    }
    */
}

void SpatiocyteStepper::setCompVoxelProperties()
{
  for(std::vector<Comp*>::iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      switch ((*i)->dimension)
        {
        case 1:
          setLineCompProperties(*i);
          setLineVoxelProperties(*i);
          break;
        case 2:
          setSurfaceCompProperties(*i);
          setSurfaceVoxelProperties(*i);
          break;
        case 3:
        default:
          setVolumeCompProperties(*i);
        }
    }
}

void SpatiocyteStepper::setLineCompProperties(Comp* aComp)
{
    setSurfaceCompProperties(aComp);
}

void SpatiocyteStepper::setLineVoxelProperties(Comp* aComp)
{
    setSurfaceVoxelProperties(aComp);
}

void SpatiocyteStepper::setSurfaceCompProperties(Comp* aComp)
{
  aComp->vacantSpecies->removePeriodicEdgeVoxels();
  aComp->vacantSpecies->removeSurfaces();
  setDiffusiveComp(aComp);
}

void SpatiocyteStepper::setVolumeCompProperties(Comp* aComp)
{
  setDiffusiveComp(aComp);
}

void SpatiocyteStepper::setSurfaceVoxelProperties(Comp* aComp)
{
  /*
  if(!aComp->diffusiveComp)
    {
      Species* aVacantSpecies(aComp->vacantSpecies);
      for(unsigned i(0); i != aVacantSpecies->size(); ++i)
        {
          unsigned aMol(aVacantSpecies->getMol(i));
          optimizeSurfaceVoxel(aMol, aComp);
          setSurfaceSubunit(aMol, aComp);
        }
    }
    */
}

void SpatiocyteStepper::setDiffusiveComp(Comp* aComp)
{
  /*
  FOR_ALL(System::Variables, aComp->system->getVariables())
    {
      Variable* aVariable(i->second);
      if(aVariable->getID() == "DIFFUSIVE")
        {
          String aStringID(aVariable->getName()); 
          aStringID = "System:" + aStringID;
          FullID aFullID(aStringID);
          System* aSystem(static_cast<System*>(getModel()->getEntity(aFullID)));
          aComp->diffusiveComp = system2Comp(aSystem);
        }
    }
  if(aComp->diffusiveComp)
    {
      Species* aVacantSpecies(aComp->vacantSpecies);
      for(unsigned i(0); i != aVacantSpecies->size(); ++i)
        {
          unsigned aMol(aVacantSpecies->getMol(i));
          aComp->diffusiveComp->vacantSpecies->addCompVoxel(aMol);
        }
      aVacantSpecies->clearMols();
    }
    */
}

void SpatiocyteStepper::optimizeSurfaceVoxel(unsigned aMol, Comp* aComp)
{
}

Species* SpatiocyteStepper::id2species(unsigned short id)
{
  return theSpecies[id];
}

unsigned short SpatiocyteStepper::getNullID()
{
  return theNullID;
}

Comp* SpatiocyteStepper::id2Comp(unsigned short id)
{
  return theSpecies[id]->getComp();
}

Comp* SpatiocyteStepper::system2Comp(System* aSystem)
{
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      if(theComps[i]->system == aSystem)
        {
          return theComps[i];
        }
    }
  return NULL;
}

void SpatiocyteStepper::setSurfaceSubunit(unsigned aMol, Comp* aComp)
{
  /*
  unsigned short& anID(theIDs[aMol]);
  // The subunit is only useful for a cylindrical surface
  // and for polymerization on it.
  anID.subunit = new Subunit;
  anID.subunit->coord = aMol;
  Point& aPoint(anID.subunit->surfacePoint);
  aPoint = coord2point(aMol);
  double aRadius(aComp->lengthY/2);
  Point aCenterPoint(aComp->centerPoint);
  Point aWestPoint(aComp->centerPoint);
  Point anEastPoint(aComp->centerPoint); 
  aWestPoint.x = aComp->centerPoint.x-aComp->lengthX/2+aComp->lengthY/2;
  anEastPoint.x = aComp->centerPoint.x+aComp->lengthX/2-aComp->lengthY/2;
  if(aPoint.x < aWestPoint.x)
    {
      aCenterPoint.x = aWestPoint.x;
    }
  else if(aPoint.x > anEastPoint.x)
    {
      aCenterPoint.x = anEastPoint.x;
    }
  else
    {
      aCenterPoint.x = aPoint.x;
    }
  double X(aPoint.x-aCenterPoint.x);
  double Y(aPoint.y-aCenterPoint.y);
  double Z(aPoint.z-aCenterPoint.z);
  double f(atan2(X, Z));
  double d(atan2(sqrt(X*X+Z*Z),Y));
  aPoint.x = aCenterPoint.x + sin(f)*aRadius*sin(d);
  aPoint.y = aCenterPoint.y + aRadius*cos(d);
  aPoint.z = aCenterPoint.z + cos(f)*aRadius*sin(d);
  */
}

void SpatiocyteStepper::setIntersectingCompartmentList() 
{
  setIntersectingPeers();
  setIntersectingParent();
}

void SpatiocyteStepper::setIntersectingPeers()
{
  for(std::vector<Comp*>::reverse_iterator i(theComps.rbegin());
      i != theComps.rend(); ++i)
    {
      //Only proceed if the volume is not enclosed:
      if(!(*i)->system->isRootSystem() && (*i)->dimension == 3)
        {
          for(std::vector<Comp*>::reverse_iterator j(theComps.rbegin());
              j != theComps.rend(); ++j)
            {
              //Only proceed if the volume of the peer is enclosed:
              if((*i) != (*j) && !(*j)->system->isRootSystem() &&
                 (*j)->dimension == 3)
                {
                  //If i and j are peer volume compartments:
                  if((*i)->system->getSuperSystem() == 
                     (*j)->system->getSuperSystem())
                    {
                      Point a((*i)->centerPoint);
                      Point b((*j)->centerPoint);
                      if(((a.x+(*i)->lengthX/2 > b.x-(*j)->lengthX/2) ||
                          (a.x-(*i)->lengthX/2 < b.x+(*j)->lengthX/2)) && 
                         ((a.y+(*i)->lengthY/2 > b.y-(*j)->lengthY/2) ||
                          (a.y-(*i)->lengthY/2 < b.y+(*j)->lengthY/2)) &&
                         ((a.z+(*i)->lengthZ/2 > b.z-(*j)->lengthZ/2) ||
                          (a.z-(*i)->lengthZ/2 < b.z+(*j)->lengthZ/2)))
                        {
                          if((*i)->enclosed <= (*j)->enclosed) 
                            {
                              (*i)->intersectPeers.push_back(*j);
                            }
                          else
                            {
                              (*i)->intersectLowerPeers.push_back(*j);
                            }
                        }
                    }
                }
            }
        }
    }
}

void SpatiocyteStepper::setIntersectingParent() 
{
  for(std::vector<Comp*>::iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      (*i)->isIntersectRoot = false;
      (*i)->isIntersectParent = false;
      //only proceed if the volume is not enclosed:
      if(!(*i)->system->isRootSystem() && (*i)->dimension == 3)
        {
          Comp* aParentComp(system2Comp((*i)->system->getSuperSystem())); 
          Point a((*i)->centerPoint);
          Point b(aParentComp->centerPoint);
          if((a.x+(*i)->lengthX/2 > b.x+aParentComp->lengthX/2) ||
             (a.x-(*i)->lengthX/2 < b.x-aParentComp->lengthX/2) || 
             (a.y+(*i)->lengthY/2 > b.y+aParentComp->lengthY/2) ||
             (a.y-(*i)->lengthY/2 < b.y-aParentComp->lengthY/2) ||
             (a.z+(*i)->lengthZ/2 > b.z+aParentComp->lengthZ/2) ||
             (a.z-(*i)->lengthZ/2 < b.z-aParentComp->lengthZ/2))
            {
              if(aParentComp->system->isRootSystem())
                {
                  (*i)->isIntersectRoot = true;
                }
              else
                {
                  (*i)->isIntersectParent = true;
                }
            }
        }
    }
}

void SpatiocyteStepper::compartmentalizeLattice() 
{
  for(unsigned i(0); i != theTotalBoxSize; ++i)
    {
      std::vector<unsigned short>& anIDs(theIDs[i]);
      for(unsigned j(0); j != anIDs.size(); ++j)
        {
          if(anIDs[j] != theNullID)
            { 
              compartmentalizeVoxel(i, j, theComps[0]);
            }
        }
    }
}

bool SpatiocyteStepper::compartmentalizeVoxel(unsigned aBox, unsigned aMol,
                                              Comp* aComp)
{
  aComp->vacantSpecies->addCompVoxel(aBox, aMol);
  return true;
  /*
  unsigned short& anID(anIDs[aMol]);
  if(aComp->dimension == 3)
    {
      if(aComp->system->isRootSystem() ||
         isInsideMol(aMol, aComp, 0))
        {
          //Check if the voxel belongs to an intersecting peer compartment:
          if(!aComp->intersectPeers.empty())
            {
              //The function isPeerMol also checks if the voxel is a future
              //surface voxel of the peer compartment: 
              if(isPeerMol(aMol, aComp))
                { 
                  return false;
                }
              if(aComp->surfaceSub && aComp->surfaceSub->enclosed)
                {
                  //Check if the voxel is neighbor of a peer voxel (either
                  //a future surface voxel)
                  if(isEnclosedSurfaceVoxel(anID, aMol, aComp))
                    {
                      aComp->surfaceSub->vacantSpecies->addCompVoxel(aMol);
                      setMinMaxSurfaceDimensions(aMol, aComp);
                      return true;
                    }
                }
            }
          if(aComp->isIntersectRoot)
            {
              Comp* aRootComp(system2Comp(aComp->system->getSuperSystem())); 
              if(aRootComp->surfaceSub && 
                 aRootComp->surfaceSub->enclosed && 
                 isRootSurfaceVoxel(anID, aMol, aRootComp))
                {
                  return false;
                }
              if(aComp->surfaceSub && 
                 isEnclosedRootSurfaceVoxel(anID, aMol, aComp, aRootComp))
                {
                  aComp->surfaceSub->vacantSpecies->addCompVoxel(aMol);
                  setMinMaxSurfaceDimensions(aMol, aComp);
                  return true;
                }
            }
          else if(aComp->isIntersectParent)
            {
              Comp* aParentComp(system2Comp(aComp->system->getSuperSystem())); 
              if(aComp->surfaceSub && aComp->surfaceSub->enclosed &&
                 isParentSurfaceVoxel(anID, aMol, aParentComp))
                {
                  aComp->surfaceSub->vacantSpecies->addCompVoxel(aMol);
                  setMinMaxSurfaceDimensions(aMol, aComp);
                  return true;
                }
            }
          for(unsigned i(0); i != aComp->immediateSubs.size(); ++i)
            {
              if(compartmentalizeVoxel(aMol, aComp->immediateSubs[i]))
                {
                  return true;
                }
            }
          aComp->vacantSpecies->addCompVoxel(aMol);
          return true;
        }
      if(aComp->surfaceSub)
        { 
          if(isInsideMol(aMol, aComp, 4) &&
             isSurfaceVoxel(anID, aMol, aComp))
            {
              aComp->surfaceSub->vacantSpecies->addCompVoxel(aMol);
              setMinMaxSurfaceDimensions(aMol, aComp);
              return true;
            }
        }
    }
  else if(aComp->system->getSuperSystem()->isRootSystem())
    {
      Comp* aRootComp(system2Comp(aComp->system->getSuperSystem())); 
      if(!isInsideMol(aMol, aRootComp, -4) &&
         isRootSurfaceVoxel(anID, aMol, aRootComp))
        {
          aComp->vacantSpecies->addCompVoxel(aMol);
          setMinMaxSurfaceDimensions(aMol, aRootComp);
          return true;
        }
    }
  return false;
  */
}

bool SpatiocyteStepper::isRootSurfaceVoxel(unsigned short& anID,
                                           unsigned aMol, Comp* aComp)
{
  /*
  for(unsigned i(0); i != theAdjoinSize; ++i)
    {
      if(theIDs[theAdjoins[aMol*theAdjoinSize+i]] == theNullID ||
         theAdjoins[aMol*theAdjoinSize+i] == aMol)
        {
          return true;
        }
    }
    */
  return false;
}

bool SpatiocyteStepper::isParentSurfaceVoxel(unsigned short& anID, unsigned aMol,
                                             Comp* aComp)
{
  /*
  if(!isInsideMol(aMol, aComp, -4))
    {
      for(unsigned i(0); i != theAdjoinSize; ++i)
        {
          if(!isInsideMol(theAdjoins[aMol*theAdjoinSize+i], aComp, 0))
            {
              return true;
            }
        }
    }
    */
  return false;
}

bool SpatiocyteStepper::isEnclosedRootSurfaceVoxel(unsigned short& anID, 
                                                   unsigned aMol,
                                                   Comp* aComp,
                                                   Comp* aRootComp)
{ 
  /*
  if(aComp->surfaceSub->enclosed || aRootComp->enclosed <= aComp->enclosed)
    { 
      if(!isInsideMol(aMol, aRootComp, -4))
        {
          if(!aRootComp->surfaceSub || 
             (aRootComp->surfaceSub && !aRootComp->surfaceSub->enclosed))
            {
              if(isRootSurfaceVoxel(anID, aMol, aComp))
                {
                  return true;
                }
            }
          else
            {
              for(unsigned i(0); i != theAdjoinSize; ++i)
                {
                  unsigned adjoinMol(theAdjoins[aMol*theAdjoinSize+i]);
                  unsigned short& adjoin(theIDs[adjoinMol]);
                  if(isRootSurfaceVoxel(adjoin, adjoinMol, aComp))
                    {
                      return true;
                    }
                }
            }
        }
    }
    */
  return false;
}

bool SpatiocyteStepper::isSurfaceVoxel(unsigned short& anID, unsigned aMol,
                                       Comp* aComp)
{
  /*
  if(isPeerMol(aMol, aComp))
    { 
      return false;
    }
  if(aComp->surfaceSub && !aComp->surfaceSub->enclosed &&
     isLowerPeerMol(aMol, aComp))
    {
      return false;
    }
  if(aComp->isIntersectRoot)
    {
      Comp* aRootComp(system2Comp(aComp->system->getSuperSystem())); 
      if(aRootComp->surfaceSub && aRootComp->surfaceSub->enclosed && 
         isRootSurfaceVoxel(anID, aMol, aRootComp))
        {
          return false;
        }
    }
  for(unsigned i(0); i != theAdjoinSize; ++i)
    {
      if(isInsideMol(theAdjoins[aMol*theAdjoinSize+i], aComp, 0))
        {
          return true;
        }
    }
    */
  return false;
}

bool SpatiocyteStepper::isEnclosedSurfaceVoxel(unsigned short& anID, 
                                               unsigned aMol,
                                               Comp* aComp)
{
  /*
  for(unsigned i(0); i != theAdjoinSize; ++i)
    {
      if(isPeerMol(theAdjoins[aMol*theAdjoinSize+i], aComp))
        {
          return true;
        }
    }
    */
  return false;
}

bool SpatiocyteStepper::isLowerPeerMol(unsigned aMol, Comp* aComp)
{
  /*
  for(std::vector<Comp*>::iterator i(aComp->intersectLowerPeers.begin());
      i != aComp->intersectLowerPeers.end(); ++i)
    {
      if(isInsideMol(aMol, *i, 0))
        {
          return true;
        }
    }
    */
  return false;
}

bool SpatiocyteStepper::isPeerMol(unsigned aMol, Comp* aComp)
{
  /*
  for(std::vector<Comp*>::iterator i(aComp->intersectPeers.begin());
      i != aComp->intersectPeers.end(); ++i)
    {
      if(isInsideMol(aMol, *i, 0))
        {
          return true;
        }
    }
  //The voxel is not inside one of the intersecting peer compartments,
  //we are now checking if the one of voxel's neighbor belongs to the
  //peer comp, to determine if anID is a surface voxel of the peer:
  for(std::vector<Comp*>::iterator i(aComp->intersectPeers.begin());
      i != aComp->intersectPeers.end(); ++i)
    {
      if((*i)->surfaceSub && (*i)->surfaceSub->enclosed)
        { 
          if(isInsideMol(aMol, *i, 4))
            {
              for(unsigned j(0); j != theAdjoinSize; ++j)
                {
                  if(isInsideMol(theAdjoins[aMol*theAdjoinSize+j], *i, 0))
                    {
                      return true;
                    }
                }
            }
        }
    }
    */
  return false;
}

bool SpatiocyteStepper::isLineVoxel(unsigned short& anID, unsigned aMol,
                                    Comp* aComp)
{
  /*
    const double safety(2.0);
    const Point aPoint(coord2point(aMol));

    double distance(aPoint.x - aComp->centerPoint.x);
    if (-safety < distance && distance <= 0)
    {
        // This is not efficient because we don't need to check volume voxels.
        // However, at this time, anID.adjoinigVoxels is not properly 
        // aligned yet.
        for(unsigned i(0); i != theAdjoinSize; ++i)
        {
            const Point& adjoinPoint(theInfo[theAdjoins[aMol*theAdjoinSize+i]].point);
            const double distance_i(adjoinPoint.x - aComp->centerPoint.x);
            if (distance_i > 0)
            {
                return true;
            }
        }
    }
                       
    distance = aPoint.y - aComp->centerPoint.y;
    if (-safety < distance && distance <= 0)
    {
        // This is not efficient because we don't need to check volume voxels.
        // However, at this time, anID.adjoinigVoxels is not properly 
        // aligned yet.
        for(unsigned i(0); i != theAdjoinSize; ++i)
        {
            const Point& adjoinPoint(theInfo[theAdjoins[aMol*theAdjoinSize+i]].point);
            const double distance_i(adjoinPoint.y - aComp->centerPoint.y);
            if (distance_i > 0)
            {
                return true;
            }
        }
    }
    */
    return false;
}

unsigned SpatiocyteStepper::getStartMol()
{
  return 0;
}


void SpatiocyteStepper::setMinMaxSurfaceDimensions(unsigned aMol, 
                                                   Comp* aComp)
{
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aMol, aRow, aLayer, aCol);
  if(aRow < aComp->minRow)
    {
      aComp->minRow = aRow;
      aComp->surfaceSub->minRow = aRow;
    }
  else if(aRow > aComp->maxRow)
    {
      aComp->maxRow = aRow;
      aComp->surfaceSub->maxRow = aRow;
    }
  if(aCol < aComp->minCol)
    {
      aComp->minCol = aCol;
      aComp->surfaceSub->minCol = aCol;
    }
  else if(aCol > aComp->maxCol)
    {
      aComp->maxCol = aCol;
      aComp->surfaceSub->maxCol = aCol;
    }
  if(aLayer < aComp->minLayer)
    {
      aComp->minLayer = aLayer;
      aComp->surfaceSub->minLayer = aLayer;
    }
  else if(aLayer > aComp->maxLayer)
    {
      aComp->maxLayer = aLayer;
      aComp->surfaceSub->maxLayer = aLayer;
    }
}



void SpatiocyteStepper::rotateX(double angle, Point* aPoint, int sign)
{ 
  if(angle)
    {
      angle *= sign;
      double y(aPoint->y);
      double z(aPoint->z);
      aPoint->y = y*cos(angle)-z*sin(angle);
      aPoint->z = y*sin(angle)+z*cos(angle);
    }
}

void SpatiocyteStepper::rotateY(double angle, Point* aPoint, int sign)
{ 
  if(angle)
    {
      angle *= sign;
      double x(aPoint->x);
      double z(aPoint->z);
      aPoint->x = x*cos(angle)+z*sin(angle);
      aPoint->z = z*cos(angle)-x*sin(angle);
    }
}

void SpatiocyteStepper::rotateZ(double angle, Point* aPoint, int sign)
{ 
  if(angle)
    {
      angle *= sign;
      double x(aPoint->x);
      double y(aPoint->y);
      aPoint->x = x*cos(angle)-y*sin(angle);
      aPoint->y = x*sin(angle)+y*cos(angle);
    }
}

bool SpatiocyteStepper::isInsideMol(unsigned aMol,
                                      Comp* aComp, double delta)
{
  Point aPoint(coord2point(aMol));
  Point aCenterPoint(aComp->centerPoint);
  Point aWestPoint(aComp->centerPoint);
  Point anEastPoint(aComp->centerPoint); 
  aPoint.x -= aCenterPoint.x;
  aPoint.y -= aCenterPoint.y;
  aPoint.z -= aCenterPoint.z;
  rotateX(aComp->rotateX, &aPoint);
  rotateY(aComp->rotateY, &aPoint);
  rotateZ(aComp->rotateZ, &aPoint);
  aPoint.x += aCenterPoint.x;
  aPoint.y += aCenterPoint.y;
  aPoint.z += aCenterPoint.z;
  double aRadius(0);
  switch(aComp->geometry)
    {
    case CUBOID:
      if(sqrt(pow(aPoint.x-aCenterPoint.x, 2)) <= 
         aComp->lengthX/2+nVoxelRadius+delta &&
         sqrt(pow(aPoint.y-aCenterPoint.y, 2)) <= 
         aComp->lengthY/2+nVoxelRadius+delta &&
         sqrt(pow(aPoint.z-aCenterPoint.z, 2)) <= 
         aComp->lengthZ/2+nVoxelRadius+delta)
        {
          return true;
        }
      break;
    case ELLIPSOID:
      //If the distance between the voxel and the center point is less than 
      //or equal to radius-2, then the voxel cannot be a surface voxel:
      if(pow(aPoint.x-aCenterPoint.x, 2)/pow((aComp->lengthX+delta)/2, 2)+ 
         pow(aPoint.y-aCenterPoint.y, 2)/pow((aComp->lengthY+delta)/2, 2)+ 
         pow(aPoint.z-aCenterPoint.z, 2)/pow((aComp->lengthZ+delta)/2, 2) <= 1)
        {
          return true;
        }
      break;
    case CYLINDER: 
      //The axial point of the cylindrical portion of the rod:
      aCenterPoint.x = aPoint.x; 
      aWestPoint.x = aComp->centerPoint.x-aComp->lengthX/2;
      anEastPoint.x = aComp->centerPoint.x+aComp->lengthX/2;
      aRadius = aComp->lengthY/2+nVoxelRadius;
      //If the distance between the voxel and the center point is less than 
      //or equal to the radius, then the voxel must be inside the Comp:
      if((aPoint.x >= aWestPoint.x && aPoint.x <= anEastPoint.x &&
          distance(aPoint, aCenterPoint) <= aRadius+delta))
        { 
          return true;
        }
      break;
    case ROD: 
      //The axial point of the cylindrical portion of the rod:
      aCenterPoint.x = aPoint.x; 
      aWestPoint.x = aComp->centerPoint.x-aComp->lengthX/2+aComp->lengthY/2;
      anEastPoint.x = aComp->centerPoint.x+aComp->lengthX/2-aComp->lengthY/2;
      aRadius = aComp->lengthY/2+nVoxelRadius;
      //If the distance between the voxel and the center point is less than 
      //or equal to the radius, then the voxel must be inside the Comp:
      if((aPoint.x >= aWestPoint.x && aPoint.x <= anEastPoint.x &&
          distance(aPoint, aCenterPoint) <= aRadius+delta) ||
         (aPoint.x < aWestPoint.x &&
          distance(aPoint, aWestPoint) <= aRadius+delta) ||
         (aPoint.x > anEastPoint.x &&
          distance(aPoint, anEastPoint) <= aRadius+delta))
        { 
          return true;
        }
      break;
    case PYRAMID: 
      aRadius = ((aCenterPoint.y+aComp->lengthY/2)-aPoint.y)/aComp->lengthY;
      if(sqrt(pow(aPoint.y-aCenterPoint.y, 2)) <= aComp->lengthY/2+delta &&
         sqrt(pow(aPoint.x-aCenterPoint.x, 2)) <= 
         aComp->lengthX*aRadius/2+delta &&
         sqrt(pow(aPoint.z-aCenterPoint.z, 2)) <=
         aComp->lengthZ*aRadius/2+delta)
        {
          return true;
        }
      break;
    case ERYTHROCYTE: 
      if(delta > 0)
        {
          return true;
        }
      else if(delta < 0)
        {
          return false;
        }
      const double Rsq(pow(aPoint.x-aCenterPoint.x, 2)/
                       pow((aComp->lengthX)/2, 2)+ 
                       pow(aPoint.y-aCenterPoint.y, 2)/
                       pow((aComp->lengthY)/2, 2));
      if(Rsq > 1)
        {
          return false;
        }
      const double a(0.5);
      const double b(0.1);
      const double R(sqrt(Rsq));
      const double thickness(((1-cos(M_PI*0.5*R))*(a-b)+b)*sqrt(1-Rsq));
      const double height((aPoint.z-aCenterPoint.z)/(2*(aComp->lengthZ)));
      if(thickness*thickness >= height*height)
        {
          return true;
        }
      break;
    }
  return false;
}

void SpatiocyteStepper::populateComp(Comp* aComp)
{
  std::vector<unsigned> populationSize;
  std::vector<Species*> prioritySpecies;
  std::vector<Species*> multiscaleSpecies;
  std::vector<Species*> diffusiveSpecies;
  std::vector<Species*> normalSpecies;
  populationSize.resize(theSpecies.size());
  for(unsigned i(0); i != theSpecies.size(); ++i)
    {
      populationSize[i] = 0;
    }
  for(std::vector<Species*>::const_iterator i(aComp->species.begin());
      i != aComp->species.end(); ++i)
    {
      if((*i)->getVacantSpecies()->getIsMultiscale())
        {
          multiscaleSpecies.push_back(*i);
        }
      else if((*i)->getVacantSpecies()->getIsDiffusiveVacant())
        {
          diffusiveSpecies.push_back(*i);
        }
      else if((*i)->getIsPopulateSpecies())
        {
          populationSize[(*i)->getVacantSpecies()->getID()] += 
            (*i)->getTotalPopulateMolSize();
          bool isPushed(false);
          std::vector<Species*> temp;
          std::vector<Species*>::const_iterator j(prioritySpecies.begin());
          while(j != prioritySpecies.end())
            {
              //Put high priority species and diffuse vacant species
              //in the high priority populate list
              if((*j)->getPopulatePriority() > (*i)->getPopulatePriority() ||
                 ((*j)->getPopulatePriority() == (*i)->getPopulatePriority() &&
                  ((*j)->getIsDiffusiveVacant() || (*j)->getIsMultiscale())))
                {
                  temp.push_back(*j);
                }
              else
                {
                  temp.push_back(*i); 
                  while(j != prioritySpecies.end())
                    {
                      temp.push_back(*j);
                      ++j;
                    }
                  isPushed = true;
                  break;
                }
              ++j;
            }
          if(!isPushed)
            {
              temp.push_back(*i);
            }
          prioritySpecies = temp;
        }
    }
  for(unsigned i(0); i != populationSize.size(); ++i)
    {
      if(populationSize[i])
        {
          //theSpecies[i] is a vacant species 
          unsigned available(theSpecies[i]->getPopulatableSize());
          if(populationSize[i] > available)
            {
              THROW_EXCEPTION(ValueError, String(
                          getPropertyInterface().getClassName()) +
                          "There are " + int2str(populationSize[i]) + 
                          " total molecules that must be uniformly " +
                          "populated,\nbut there are only "
                          + int2str(available) + " vacant voxels of [" + 
                          theSpecies[i]->getVariable()->getFullID().asString() +
                          "] that can be populated on.");
            } 
          if(double(populationSize[i])/available > 0.2)
            { 
              populateSpeciesDense(prioritySpecies, populationSize[i],
                                   available);
            }
          else
            {
              populateSpeciesSparse(prioritySpecies);
            }
        }
    }
  for(std::vector<Species*>::const_iterator i(multiscaleSpecies.begin());
      i != multiscaleSpecies.end(); ++i)
    {
      (*i)->populateUniformOnMultiscale();
    }
  for(std::vector<Species*>::const_iterator i(diffusiveSpecies.begin());
      i != diffusiveSpecies.end(); ++i)
    {
      (*i)->populateUniformOnDiffusiveVacant();
    }
}

void SpatiocyteStepper::populateSpeciesDense(std::vector<Species*>&
                                             aSpeciesList, unsigned aSize,
                                             unsigned availableVoxelSize)
{
  unsigned count(0);
  unsigned* populateVoxels(new unsigned[aSize]);
  unsigned* availableVoxels(new unsigned [availableVoxelSize]); 
  for(unsigned i(0); i != availableVoxelSize; ++i)
    {
      availableVoxels[i] = i;
    }
  gsl_ran_choose(getRng(), populateVoxels, aSize, availableVoxels,
                 availableVoxelSize, sizeof(unsigned));
  //gsl_ran_choose arranges the position ascending, so we need
  //to shuffle the order of voxel positions:
  gsl_ran_shuffle(getRng(), populateVoxels, aSize, sizeof(unsigned)); 
  for(std::vector<Species*>::const_iterator i(aSpeciesList.begin());
      i != aSpeciesList.end(); ++i)
    {
      (*i)->populateCompUniform(populateVoxels, &count);
    }
  delete[] populateVoxels;
  delete[] availableVoxels;
}

void SpatiocyteStepper::populateSpeciesSparse(std::vector<Species*>&
                                              aSpeciesList)
{
  for(std::vector<Species*>::const_iterator i(aSpeciesList.begin());
      i != aSpeciesList.end(); ++i)
    {
      (*i)->populateCompUniformSparse();
    }
}


void SpatiocyteStepper::clearComp(Comp* aComp)
{
  for(std::vector<Species*>::const_iterator i(aComp->species.begin());
      i != aComp->species.end(); ++i)
    {
      (*i)->removeMols();
      (*i)->updateMols();
    }
}


std::vector<Comp*> const& SpatiocyteStepper::getComps() const
{
  return theComps;
}


