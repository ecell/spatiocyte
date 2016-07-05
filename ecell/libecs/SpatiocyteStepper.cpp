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
#include <unistd.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#endif
#include <gsl/gsl_randist.h>
#include <libecs/Model.hpp>
#include <libecs/System.hpp>
#include <libecs/Process.hpp>
#include <libecs/Stepper.hpp>
#include <libecs/VariableReference.hpp>
#include <libecs/SpatiocyteStepper.hpp>
#include <libecs/SpatiocyteSpecies.hpp>
#include <libecs/CompartmentProcess.hpp>

LIBECS_DM_INIT_STATIC(SpatiocyteStepper, Stepper);

size_t get_total_system_memory()
{
  long mem_size;
#ifdef __APPLE__
  size_t len(sizeof(mem_size));
  int ret(-1);
  /* Note hw.memsize is in bytes, so no need to multiply by page size. */
  sysctlbyname("hw.memsize", &mem_size, &len, NULL, 0);
#else
  long pages = sysconf(_SC_PHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  mem_size = pages * page_size;
#endif
  return mem_size;
}

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
  cout.setLevel(DebugLevel);
  cout << "1. checking model..." << std::endl;
  checkModel();
  //We need a Comp tree to assign the voxels to each Comp
  //and get the available number of vacant voxels. The compartmentalized
  //vacant voxels are needed to randomly place molecules according to the
  //Comp:
  cout << "2. creating compartments..." << std::endl;
  registerComps();
  setCompsProperties();
  cout << "3. setting up lattice properties..." << std::endl;
  setLatticeProperties(); 
  broadcastLatticeProperties();
  setCompsCenterPoint();
  //All species have been created at this point, we initialize them now:
  cout << "4. initializing species..." << std::endl;
  initSpecies();
  cout << "5. initializing processes the first time..." << std::endl;
  initializeFirst();
  cout << "6. initializing processes the second time..." << std::endl;
  initializeSecond();
  cout << "7. constructing lattice..." << std::endl;
  constructLattice();
  cout << "8. setting intersecting compartment list..." << std::endl;
  setIntersectingCompartmentList();
  cout << "9. compartmentalizing lattice..." << std::endl;
  compartmentalizeLattice();
  setCuboidCompsCenterPoint();
  cout << "10. setting up compartment voxels properties..." << std::endl;
  setCompVoxelProperties();
  resizeProcessLattice();
  cout << "11. initializing compartments..." << std::endl;
  initializeCompartments();
  cout << "12. storing simulation parameters..." << std::endl;
  storeSimulationParameters();
  initializeCompartmentsOnce();
  cout << "13. initializing processes the third time..." << std::endl;
  initializeThird();
  cout << "14. updating species..." << std::endl;
  //Cannot put updateSpecies before initializeThird because SNRP declares
  //a species isReactiveVacant in initializeThird
  updateSpecies();
  cout << "15. printing simulation parameters..." << std::endl;
  printSimulationParameters();
  cout << "16. initializing before populate..." << std::endl;
  initializeBeforePopulate();
  cout << "17. populating compartments with molecules..." << std::endl;
  populateComps();
  cout << "18. initializing processes the fourth time..." << std::endl;
  initializeFourth();
  cout << "19. initializing the priority queue..." << std::endl;
  initPriorityQueue();
  cout << "20. initializing processes the fifth time..." << std::endl;
  initializeFifth();
  cout << "21. initializing processes the last time..." << std::endl;
  initializeLastOnce();
  cout << "22. finalizing species..." << std::endl;
  finalizeSpecies();
  cout << "23. printing final process parameters..." << std::endl <<
    std::endl;
  printProcessParameters();
  cout << "24. simulation is started..." << std::endl;
}


unsigned SpatiocyteStepper::getNewMoleculeID()
{
  ++theMoleculeID;
  return theMoleculeID-1;
}

void SpatiocyteStepper::interrupt(double aTime)
{
  if(theExternInterruptedProcesses.size())
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
}

void SpatiocyteStepper::finalizeSpecies()
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      (*i)->finalizeSpecies();
    }
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
  return theRowSize;
}

unsigned SpatiocyteStepper::getLayerSize()
{
  return theLayerSize;
}

unsigned SpatiocyteStepper::getColSize()
{
  return theColSize;
}

unsigned SpatiocyteStepper::getLatticeSize()
{
  return theLattice.size();
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
  theRan.Reseed();
  setCurrentTime(0);
  resetVariables();
  initializeSecond();
  clearComps();
  initializeCompartments();
  initializeThird();
  initializeBeforePopulate();
  populateComps();
  initializeFourth();
  initPriorityQueue();
  initializeFifth();
  finalizeSpecies();
  //checkLattice();
}

void SpatiocyteStepper::resetVariables()
{
  for(unsigned i(0); i != theVariables.size(); ++i)
    {
      if(variable2species(theVariables[i]) == NULL)
        {
          theVariables[i]->setValue(initVariableValues[i]);
        }
    }
}

Species* SpatiocyteStepper::addSpecies(Variable* aVariable)
{
  std::vector<Species*>::iterator aSpeciesIter(variable2ispecies(aVariable));
  if(aSpeciesIter == theSpecies.end())
    {
      Species *aSpecies(new Species(this, aVariable, theSpecies.size(),
                          (int)aVariable->getValue(), theRan, VoxelRadius,
                          theLattice, theSpecies));
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

bool SpatiocyteStepper::isBoundaryCoord(unsigned aCoord,
                                        unsigned aDimension)
{
  //This method is only for checking boundaries on a cube:
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aCoord, aRow, aLayer, aCol);
  if(aDimension == 3)
    {
      //If the voxel is on one of the 6 cuboid surfaces:
      if(aRow == 1 || aLayer == 1 || aCol == 1 ||
         aRow == theRowSize-2 || aLayer == theLayerSize-2 || 
         aCol == theColSize-2)
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
         (aRow <= 1 && aCol >= theColSize-2) ||
         (aRow <= 1 && aLayer >= theLayerSize-2) ||
         (aCol <= 1 && aRow >= theRowSize-2) ||
         (aCol <= 1 && aLayer >= theLayerSize-2) ||
         (aLayer <= 1 && aRow >= theRowSize-2) ||
         (aLayer <= 1 && aCol >= theColSize-2) ||
         (aRow >= theRowSize-2 && aCol >= theColSize-2) ||
         (aRow >= theRowSize-2 && aLayer >= theLayerSize-2) ||
         (aCol >= theColSize-2 && aLayer >= theLayerSize-2))
        {
          return true;
        }
    }
  return false;
}

unsigned SpatiocyteStepper::getPeriodicCoord(unsigned aCoord,
                                                 unsigned aDimension,
                                                 Origin* anOrigin)
{
  //This method is only for checking boundaries on a cube:
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aCoord, aRow, aLayer, aCol);
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
      nextRow = theRowSize-(3+adj);
      --anOrigin->row;
    }
  else if(aRow == theRowSize-(2+adj))
    {
      nextRow = 2+adj;
      ++anOrigin->row;
    }
  if(aLayer == 1+adj)
    {
      nextLayer = theLayerSize-(3+adj);
      --anOrigin->layer;
    }
  else if(aLayer == theLayerSize-(2+adj))
    {
      nextLayer = 2+adj;
      ++anOrigin->layer;
    }
  if(aCol == 1+adj)
    {
      nextCol = theColSize-(3+adj);
      --anOrigin->col;
    }
  else if(aCol == theColSize-(2+adj))
    {
      nextCol = 2+adj;
      ++anOrigin->col;
    }
  if(nextRow != aRow || nextCol != aCol || nextLayer != aLayer)
    {
      return nextRow+theRowSize*nextLayer+theRowSize*theLayerSize*nextCol;
    }
  return theNullCoord;
}

Point SpatiocyteStepper::getPeriodicPoint(unsigned aCoord,
                                          unsigned aDimension,
                                          Origin* anOrigin)
{
  unsigned adj(1);
  if(aDimension == 3)
    {
      adj = 0;
    }
  unsigned row(theRowSize-(3+adj)-(1+adj));
  unsigned layer(theLayerSize-(3+adj)-(1+adj));
  unsigned col(theColSize-(3+adj)-(1+adj));
  unsigned aGlobalCol;
  unsigned aGlobalLayer;
  unsigned aGlobalRow;
  coord2global(aCoord, aGlobalRow, aGlobalLayer, aGlobalCol);
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

Species* SpatiocyteStepper::variable2species(const Variable* aVariable) const
{
  for(std::vector<Species*>::const_iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      if((*i)->getVariable() == aVariable)
        {
          return (*i);
        }
    }
  return NULL;
}

void SpatiocyteStepper::enlistModelVariables()
{
  Model::StepperMap aStepperMap(getModel()->getStepperMap());  
  for(Model::StepperMap::const_iterator i(aStepperMap.begin());
      i != aStepperMap.end(); ++i )
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
              const VariableReference& aVariableReference(*k);
              Variable* aVariable(aVariableReference.getVariable()); 
              if(std::find(theVariables.begin(), theVariables.end(),
                           aVariable) == theVariables.end())
                {
                  theVariables.push_back(aVariable);
                  initVariableValues.push_back(aVariable->getValue());
                }
            }
        }
    }
}

void SpatiocyteStepper::checkModel()
{
  std::vector<Process*> aProcessList;
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
              if(dynamic_cast<SpatiocyteProcess*>(*j))
                {
                  aProcessList.push_back(*j);
                }
              else
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
  for(std::vector<Process*>::const_iterator i(aProcessList.begin());
      i != aProcessList.end(); ++i)
    {
      theProcessVector.push_back(*i);
    }
  enlistModelVariables();
}

void SpatiocyteStepper::checkLattice()
{
  std::vector<unsigned> list;
  for(unsigned i(0); i!=theSpecies.size(); ++i)
    {
      list.push_back(0);
    }
  for(unsigned i(0); i!=theLattice.size(); ++i)
    {
      ++list[getID(theLattice[i])];
    }
  for(unsigned i(0); i!=list.size(); ++i)
    {
      if(!theSpecies[i]->getIsCompVacant() && 
         theSpecies[i]->getVariable() != NULL)
        {
          if(theSpecies[i]->size() != list[i])
            {
              cout << theSpecies[i]->getIDString() << " in lattice:" <<
                list[i] << ", but actual size:" << theSpecies[i]->size() <<
                std::endl;
            }
        }
    }
}

void SpatiocyteStepper::initSpecies()
{
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i != theSpecies.end(); ++i)
    {
      (*i)->initialize(theAdjoiningCoordSize, theNullCoord, theNullID);
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
      theSpatiocyteProcesses.push_back(aProcess);
      aProcess->setLatticeProperties(&theLattice, theAdjoiningCoordSize,
                                     theNullCoord, theNullID, &theRan);
    }
}

void SpatiocyteStepper::initializeFirst()
{
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      theSpatiocyteProcesses[i]->initializeFirst();
    }
}

void SpatiocyteStepper::initializeSecond()
{
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      theSpatiocyteProcesses[i]->initializeSecond();
    }
}

void SpatiocyteStepper::initializeCompartments()
{
  for(std::vector<SpatiocyteProcess*>::const_iterator 
      i(theSpatiocyteProcesses.begin()); i != theSpatiocyteProcesses.end(); ++i)
    {      
      CompartmentProcess* aCompProcess(dynamic_cast<CompartmentProcess*>(*i));
      if(aCompProcess)
        {
          aCompProcess->initializeCompartment();
        }
    }
}

void SpatiocyteStepper::initializeCompartmentsOnce()
{
  for(std::vector<SpatiocyteProcess*>::const_iterator 
      i(theSpatiocyteProcesses.begin()); i != theSpatiocyteProcesses.end(); ++i)
    {      
      CompartmentProcess* aCompProcess(dynamic_cast<CompartmentProcess*>(*i));
      if(aCompProcess)
        {
          aCompProcess->initializeCompartmentOnce();
        }
    }
}


void SpatiocyteStepper::initializeThird()
{
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      theSpatiocyteProcesses[i]->initializeThird();
    }
}

void SpatiocyteStepper::initializeBeforePopulate()
{
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      theSpatiocyteProcesses[i]->initializeBeforePopulate();
    }
}

void SpatiocyteStepper::initializeFourth()
{
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      theSpatiocyteProcesses[i]->initializeFourth();
    }
}

void SpatiocyteStepper::initializeFifth()
{
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      theSpatiocyteProcesses[i]->initializeFifth();
    }
  setStepInterval(thePriorityQueue.getTop()->getTime()-getCurrentTime());
}

void SpatiocyteStepper::initializeLastOnce()
{
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      theSpatiocyteProcesses[i]->initializeLastOnce();
    }
}


void SpatiocyteStepper::printProcessParameters()
{
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      theSpatiocyteProcesses[i]->printParameters();
    }
  cout << std::endl;
}

void SpatiocyteStepper::resizeProcessLattice()
{
  unsigned startCoord(theLattice.size());
  unsigned endCoord(startCoord);
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      endCoord += theSpatiocyteProcesses[i]->getLatticeResizeCoord(endCoord);
    }
  //Save the coords of molecules before resizing theLattice
  //because the voxel address pointed by molecule list will become
  //invalid once it is resized:
  for(unsigned i(0); i != theSpecies.size(); ++i)
    {
      theSpecies[i]->saveCoords();
    }
  theLattice.resize(endCoord);
  for(unsigned i(startCoord); i != endCoord; ++i)
    {
      theLattice[i].coord = i;
    }
  //Update the molecule list with the new voxel addresses:
  for(unsigned i(0); i != theSpecies.size(); ++i)
    {
      theSpecies[i]->updateMoleculePointers();
    }
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {      
      theSpatiocyteProcesses[i]->updateResizedLattice();
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
                                      aSpatiocyteProcess->getNewInterval());
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
      ReactionProcess* aReactionProcess(dynamic_cast<ReactionProcess*>(*i));
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

SpatiocyteProcess* SpatiocyteStepper::getCurrentProcess()
{
  return thePriorityQueue.getTop();
}

inline void SpatiocyteStepper::step()
{
  do
    {
      //cout << "before:" << thePriorityQueue.getTop()->getIDString() << " "
      //<< getCurrentTime() << std::endl;
      thePriorityQueue.getTop()->fire();
      //checkSpecies();
      //checkLattice();
      if(thePriorityQueue.getTop()->getTime() != getCurrentTime())
        {
          interruptProcesses(getCurrentTime());
        }
    }
  while(thePriorityQueue.getTop()->getTime() == getCurrentTime());
  setNextTime(thePriorityQueue.getTop()->getTime());
}

void SpatiocyteStepper::addInterruptedProcess(SpatiocyteProcess* aProcess)
{
  if(std::find(theInterruptedProcesses.begin(), theInterruptedProcesses.end(),
               aProcess) == theInterruptedProcesses.end())
    {
      theInterruptedProcesses.push_back(aProcess);
    }
}

void SpatiocyteStepper::interruptProcesses(const double aCurrentTime)
{
  for(unsigned i(0); i != theInterruptedProcesses.size(); ++i)
    {
      theInterruptedProcesses[i]->substrateValueChanged(aCurrentTime);
    }
  theInterruptedProcesses.resize(0);
}

void SpatiocyteStepper::interruptAllProcesses(const double aCurrentTime)
{
  for(unsigned i(0); i != theSpatiocyteProcesses.size(); ++i)
    {
      theSpatiocyteProcesses[i]->substrateValueChanged(aCurrentTime);
    }
}

void SpatiocyteStepper::checkSpecies()
{
  for(unsigned i(0); i != theSpecies.size(); ++i)
    {
      Species* aSpecies(theSpecies[i]);
      if(!aSpecies->getIsCompVacant())
        {
          std::vector<unsigned> boundCnts;
          if(aSpecies->getIsDeoligomerize())
            {
              boundCnts.resize(7);
            }
          for(unsigned j(0); j != aSpecies->size(); ++j)
            {
              Voxel* aVoxel(aSpecies->getMolecule(j));
              if(aVoxel->idx/theStride != aSpecies->getID())
                {
                  cout << aSpecies->getIDString() << " size:" << 
                    aSpecies->size() << std::endl;
                  cout << "wrong id, mine:" << aSpecies->getID() <<
                    " but in the list (idx/stride):" << 
                    aSpecies->getIDString(aVoxel->idx/theStride) << " idx:" <<
                    aVoxel->idx  << " stride:" << theStride << std::endl;
                }
              if(aVoxel->idx-aSpecies->getID()*theStride != j)
                {
                  cout << aSpecies->getIDString() << " size:" << 
                    aSpecies->size() << std::endl;
                  cout << "wrong index, mine:" << j <<
                    " but in the list (idx-id*stride):" << aVoxel->idx-
                    aSpecies->getID()*theStride << " idx:" << aVoxel->idx
                    << " stride:" << theStride << std::endl;
                }
              if(aSpecies->getIsOnMultiscale())
                {
                  if(!aSpecies->getTag(j)[0].multiIdx)
                    {
                      cout << aSpecies->getIDString() << " size:" << 
                        aSpecies->size() << " vacidx wrong:" << j << std::endl;
                    }
                }
              if(aSpecies->getTag(j)[0].boundCnt)
                {
                  boundCnts[aSpecies->getTag(j)[0].boundCnt]++;
                  unsigned cnt(0);
                  Voxel* aVoxel(aSpecies->getMolecule(j));
                  for(unsigned k(0); k != aVoxel->diffuseSize; ++k)
                    {
                      unsigned coord(aVoxel->adjoiningCoords[k]);
                      if(theSpecies[aSpecies->getID(theLattice[
                                      coord])]->getIsDeoligomerize())
                        {
                          cnt++;
                        }
                    }
                  if(cnt != aSpecies->getTag(j)[0].boundCnt)
                    {
                      cout << aSpecies->getIDString() << " j:" << j <<
                        "error in deoligomerize cnt:" << cnt <<
                        " expected:" << aSpecies->getTag(j)[0].boundCnt <<
                        std::endl;
                    }
                }
            }
          for(unsigned j(0); j < boundCnts.size(); ++j)
            {
              if(boundCnts[j] != aSpecies->getBoundCnts()[j])
                {
                  cout << aSpecies->getIDString() << " j:" << j <<
                    "error in deoligomerize total boundCnts:" << 
                    aSpecies->getBoundCnts()[j] <<
                    " expected:" << boundCnts[j] << std::endl;
                }
            }
        }
      else
        {
          for(unsigned j(0); j != aSpecies->size(); ++j)
            {
              Voxel* aVoxel(aSpecies->getMolecule(j));
              unsigned anID(aVoxel->idx/theStride);
              if(anID != aSpecies->getID() && 
                 !theSpecies[anID]->getIsMultiscale())
                {
                  unsigned anIndex(aVoxel->idx%theStride);
                  if(theSpecies[anID]->getMolecule(anIndex) != aVoxel)
                    {
                      cout << "error in index j:" << j << " " <<
                        theSpecies[anID]->getIDString() << " size:"
                        << theSpecies[anID]->size() << " index:" << anIndex
                        << " idx:" << aVoxel->idx << std::endl;
                    }
                }
            }
        }
    }
}

void SpatiocyteStepper::registerComps()
{
  theRan.Reseed();
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
  Species* aSpecies(new Species(this, NULL, theSpecies.size(), 0, theRan,
                                VoxelRadius, theLattice, theSpecies));
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
  aComp->isProcessComp = false;
  aComp->minCoord.row = UINT_MAX;
  aComp->minCoord.col = UINT_MAX;
  aComp->minCoord.layer = UINT_MAX;
  aComp->maxCoord.row = 0;
  aComp->maxCoord.col = 0;
  aComp->maxCoord.layer = 0;
  aComp->lengthVector = Point(1, 0, 0);
  aComp->widthVector = Point(0, 1, 0);
  aComp->heightVector = Point(0, 0, 1);
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
  aComp->vacantSpecies = NULL;
  //Default Comp geometry is Cuboid:
  aComp->geometry = 0;
  //Default is volume Comp:
  aComp->dimension = 3;
  if(getVariable(aSystem, "PROCESS_COMPARTMENT"))
    {
      aComp->isProcessComp = true;
    }
  else
    {
      setCompAttributes(aSystem, aComp);
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
          aSubComp->lengthVector = aComp->lengthVector;
          aSubComp->widthVector = aComp->widthVector;
          aSubComp->heightVector = aComp->heightVector;
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
              lineComp->lengthVector = aComp->lengthVector;
              lineComp->widthVector = aComp->widthVector;
              lineComp->heightVector = aComp->heightVector;
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

void SpatiocyteStepper::setCompAttributes(System* aSystem, Comp* aComp)
{
  if(getVariable(aSystem, "DIMENSION"))
    { 
      const double aDimension(
          aSystem->getVariable("DIMENSION")->getValue());
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
}

void SpatiocyteStepper::setCompsProperties()
{
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      if(!theComps[i]->isProcessComp)
        {
          setCompProperties(theComps[i]);
        }
    }
}

void SpatiocyteStepper::setCompsCenterPoint()
{
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      if(!theComps[i]->isProcessComp)
        {
          setCompCenterPoint(theComps[i]);
        }
    }
}

void SpatiocyteStepper::setCuboidCompsCenterPoint()
{
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      Comp* aComp(theComps[i]); 
      if(aComp->geometry != CUBOID || aComp->rotateX != 0 || 
         aComp->rotateY != 0 || aComp->rotateZ != 0)
        {
          continue;
        }
      const unsigned minCoord(global2coord(aComp->minCoord.row,
                                           aComp->minCoord.layer,
                                           aComp->minCoord.col));
      const unsigned maxCoord(global2coord(aComp->maxCoord.row,
                                           aComp->maxCoord.layer,
                                           aComp->maxCoord.col));
      Point maxPoint(coord2point(maxCoord));
      Point minPoint(coord2point(minCoord)); 
      aComp->centerPoint = add(minPoint, subDivide(maxPoint, minPoint));
    }
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
  if(aComp->vacantSpecies == NULL) {
    THROW_EXCEPTION(NotFound, "Variable VACANT of the compartment "
                          + aSystem->getFullID().asString() + ", which "
                          + "specifies the vacant species, "
                          + "not defined." );
  }
}

void SpatiocyteStepper::setLatticeProperties()
{
  theStride = UINT_MAX/theSpecies.size();
  Comp* aRootComp(theComps[0]);
  rotateCompartment(aRootComp);
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      theAdjoiningCoordSize = 12;
      theHCPl = nVoxelRadius/sqrt(3); 
      theHCPx = nVoxelRadius*sqrt(8.0/3); //Lx
      theHCPy = nVoxelRadius*sqrt(3); //Ly
      break;
    case CUBIC_LATTICE:
      theAdjoiningCoordSize = 6;
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
      theRowSize = (unsigned)rint((theCenterPoint.z)/
                                      (nVoxelRadius));
      theLayerSize = (unsigned)rint((theCenterPoint.y*2)/theHCPy);
      theColSize = (unsigned)rint((theCenterPoint.x*2)/theHCPx);
      break;
    case CUBIC_LATTICE:
      theRowSize = (unsigned)rint((theCenterPoint.z)/
                                      (nVoxelRadius));
      theLayerSize = (unsigned)rint((theCenterPoint.y)/
                                      (nVoxelRadius));
      theColSize = (unsigned)rint((theCenterPoint.x)/
                                      (nVoxelRadius));
      break;
    }
  //For the CUBOID cell geometry, we need to readjust the size of
  //row, layer and column according to the boundary condition of its surfaces
  //to reflect the correct volume. This is because periodic boundary will
  //consume a layer of the surface voxels:
  if(aRootComp->geometry == CUBOID || aRootComp->geometry == OPEN_CYLINDER)
    {
      //We need to increase the row, layer and col size by 2 because
      //the entire volume must be surrounded by nullID voxels to avoid
      //self-homodimerization reaction.
      theRowSize += 2;
      theColSize += 2;
      theLayerSize += 2;
      readjustSurfaceBoundarySizes(); 
    }
  double mem_required(((theRowSize*theLayerSize*theColSize+1)*108.0)/
                      (1024.0*1024.0*1024.0));
  double mem_available(get_total_system_memory()/(1024.0*1024.0*1024.0));
  cout << "   Estimated total memory required (GB): " << mem_required 
    << std::endl;
  cout << "   Total memory available (GB): " << mem_available <<
    std::endl;
  if (mem_required > mem_available) {
    THROW_EXCEPTION(InitializationFailed,
                    getPropertyInterface().getClassName() + 
                    ": Not enough system memory to support the total " +
                    "number of voxels. Try increasing the VoxelRadius to " +
                    "reduce the number of voxels.");
  }

  //You should only resize theLattice once here. You should never
  //push or readjust the capacity of theLattice since all the pointers
  //to the voxel will become invalid once you do that:
  //Also add one more voxel for the nullVoxel:
  theLattice.resize(theRowSize*theLayerSize*theColSize+1);
  //Initialize the null coord:
  theNullCoord = theRowSize*theLayerSize*theColSize;
  theLattice[theNullCoord].idx = theNullID*theStride;
}

void SpatiocyteStepper::rotateCompartment(Comp* aComp)
{
  rotate(aComp, aComp->lengthVector);
  rotate(aComp, aComp->widthVector);
  rotate(aComp, aComp->heightVector);
  Point min(-aComp->lengthX/2, -aComp->lengthY/2, -aComp->lengthZ/2);
  Point max(aComp->lengthX/2, aComp->lengthY/2, aComp->lengthZ/2);
  rotate(aComp, min);
  rotate(aComp, max);
  aComp->lengthX = std::max(min.x, max.x)-std::min(min.x, max.x);
  aComp->lengthY = std::max(min.y, max.y)-std::min(min.y, max.y);
  aComp->lengthZ = std::max(min.z, max.z)-std::min(min.z, max.z);
  aComp->nLength = aComp->lengthX;
  aComp->nWidth = aComp->lengthY;
  aComp->nHeight = aComp->lengthZ;
  //for cylinder and rod compartment need to create lengthVector and
  //rotate it, only then it would work.
}

void SpatiocyteStepper::rotate(Comp* aComp, Point& V, int sign)
{
  rotateX(aComp->rotateX, &V, sign);
  rotateY(aComp->rotateY, &V, sign);
  rotateZ(aComp->rotateZ, &V, sign);
}

void SpatiocyteStepper::storeSimulationParameters()
{
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      Comp* aComp(theComps[i]); 
      if(aComp->isProcessComp)
        {
          continue;
        }
      if(aComp->dimension == 2)
        {
          /*
          aComp->actualArea =  (72*pow(VoxelRadius,2))*
            aComp->vacantSpecies->size()/(6*pow(2,0.5)+4*pow(3,0.5)+
                                         3*pow(6, 0.5));
                                         */
          aComp->actualArea =  2*pow(VoxelRadius,2)*sqrt(3)*
            aComp->vacantSpecies->size();
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
      createVariable("SIZE", aSystem, getModel())->setValue(aSize);
    }
}

Variable* SpatiocyteStepper::createVariable(const String& anID,
                                            const System* aSystem,
                                            Model* aModel)
{
  String anEntityType("Variable");
  SystemPath aSystemPath(aSystem->getSystemPath());
  aSystemPath.push_back(aSystem->getID());
  FullID aFullID(anEntityType, aSystemPath, anID);
  Variable* aVariable(reinterpret_cast<Variable*>(
                              aModel->createEntity("Variable", aFullID)));
  aVariable->setValue(0);
  return aVariable;
}

Process* SpatiocyteStepper::createProcess(const String& aClassname,
                                          const String& anID,
                                          const System* aSystem,
                                          Model* aModel)
{
  String anEntityType("Process");
  SystemPath aSystemPath(aSystem->getSystemPath());
  aSystemPath.push_back(aSystem->getID());
  FullID aFullID(anEntityType, aSystemPath, anID);
  return reinterpret_cast<Process*>(aModel->createEntity(aClassname, aFullID));
}


Species* SpatiocyteStepper::createSpecies(System* aSystem, String anID)
{
  Variable* aVariable(createVariable(anID, aSystem, getModel()));
  return addSpecies(aVariable);
}

void SpatiocyteStepper::printSimulationParameters()
{
  cout << std::endl;
  cout << "   Voxel radius, r_v:" << VoxelRadius << " m" << std::endl;
  cout << "   Simulation height:" << theCenterPoint.y*2*VoxelRadius*2 <<
    " m" << std::endl;
  cout << "   Simulation width:" << theCenterPoint.z*2*VoxelRadius*2 << 
    " m" << std::endl;
  cout << "   Simulation length:" << theCenterPoint.x*2*VoxelRadius*2 <<
    " m" << std::endl;
  cout << "   Layer size:" << theLayerSize << std::endl;
  cout << "   Row size:" << theRowSize << std::endl;
  cout << "   Column size:" << theColSize << std::endl;
  cout << "   Total allocated voxels:" << 
    theRowSize*theLayerSize*theColSize << std::endl;
  for(unsigned i(0); i != theComps.size(); ++i)
    {
      Comp* aComp(theComps[i]);
      if(aComp->isProcessComp)
        {
          continue;
        }
      double aSpecVolume(aComp->specVolume);
      double aSpecArea(aComp->specArea);
      double aSpecLength(aComp->specLength);
      double anActualVolume(aComp->actualVolume);
      double anActualArea(aComp->actualArea);
      double anActualLength(aComp->actualLength);
      switch(aComp->geometry)
        {
        case CUBOID:
          cout << "   Cuboid ";
          break;
        case ELLIPSOID:
          cout << "   Ellipsoid ";
          break;
        case CYLINDER:
        case OPEN_CYLINDER:
          cout << "   Cylinder (radius=" << aComp->lengthY*VoxelRadius
            << "m, length=" << (aComp->lengthX)*VoxelRadius*2 << "m) ";
          break;
        case ROD:
          cout << "   Rod (radius=" << aComp->lengthY*VoxelRadius << 
            "m, cylinder length=" <<
            (aComp->lengthX-aComp->lengthY)*VoxelRadius*2 << "m) ";
          break;
        }
      cout << aComp->system->getFullID().asString();
      switch (aComp->dimension)
      { 
      case 1:
          cout << " Filament compartment:" << std::endl;
          cout << "     [" << aSpecLength << " m] Specified filament length " <<
            std::endl;
          cout << "     [" << anActualLength << " m] Actual filament area " <<
            std::endl;
          break;
      case 2:
          cout << " Surface compartment:" << std::endl;
          cout << "     [" << int(aSpecArea*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))/
                              (72*VoxelRadius*VoxelRadius)) << 
            "] Specified surface voxels {n_s = S_specified*"
            << "(6*2^0.5+4*3^0.5+3*6^0.5)/(72*r_v^2}" << std::endl;
          cout << "     [" << aComp->vacantSpecies->size() <<
            "] Actual surface voxels {n_s}" << std::endl;
          cout << "     [" << aSpecArea << " m^2] Specified surface area " <<
            "{S_specified}" << std::endl;
          cout << "     [" << anActualArea << " m^2] Actual surface area " <<
            "{S = (72*r_v^2)*n_s/(6*2^0.5+4*3^0.5+3*6^0.5)}" << std::endl;
          break;
      case 3:
      default:
          cout << " Volume compartment:" << std::endl;
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
          cout << "     [" << int(aSpecVolume/(4*sqrt(2)*pow(VoxelRadius, 3))) << 
            "] Specified volume voxels {n_v = V_specified/(4*2^0.5*r_v^3)}" <<
          std::endl;  
          cout << "     [" << voxelCnt << "] Actual volume voxels {n_v}"  << std::endl;
          cout << "     [" << aSpecVolume << " m^3] Specified volume {V_specified}"
            << std::endl; 
          cout << "     [" << anActualVolume << " m^3] Actual volume " <<
            "{V = (4*2^0.5*r_v^3)*n_v}" << std::endl; 
      }
    }
  cout << std::endl;
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
  if(aRootComp->geometry == CUBOID && aRootComp->surfaceSub)
    {
      if(aRootComp->xyPlane == REFLECTIVE)
        {
          theRowSize += 2;
        }
      else if(aRootComp->xyPlane == REMOVE_UPPER ||
              aRootComp->xyPlane == REMOVE_LOWER)
        {
          theRowSize += 1;
        }
      if(aRootComp->yzPlane == REFLECTIVE)
        {
          theColSize += 2;
        }
      else if(aRootComp->yzPlane == REMOVE_UPPER ||
              aRootComp->yzPlane == REMOVE_LOWER)
        {
          theColSize += 1;
        }
      if(aRootComp->xzPlane == REFLECTIVE)
        {
          theLayerSize += 2;
        }
      else if(aRootComp->xzPlane == REMOVE_UPPER ||
              aRootComp->xzPlane == REMOVE_LOWER)
        {
          theLayerSize += 1;
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
          if(theColSize%2 != 1)
            {
              theColSize += 1;
            }
          else
            {
              theColSize += 2;
            }
        }
      if(aRootComp->xzPlane == PERIODIC)
        {
          if(theLayerSize%2 != 1)
            {
              theLayerSize +=1;
            }
          else
            {
              theLayerSize += 2;
            }
        }
      if(aRootComp->xyPlane == PERIODIC)
        {
          if(theRowSize%2 != 1)
            {
              theRowSize += 1;
            }
          else
            {
              theRowSize += 2;
            }
        }
    }
  if(isPeriodicEdge)
    {
      if(theColSize%2 == 1)
        {
          theColSize += 1;
        }
      if(theLayerSize%2 == 1)
        {
          theLayerSize +=1;
        }
      if(theRowSize%2 == 1)
        {
          theRowSize += 1;
        }
    }
}

void SpatiocyteStepper::constructLattice()
{ 
  Comp* aRootComp(theComps[0]);
  unsigned aSize(theRowSize*theLayerSize*theColSize);
  unsigned a(0);
  unsigned short rootID(aRootComp->vacantSpecies->getID());
  for(std::vector<Voxel>::iterator i(theLattice.begin()); a != aSize; ++i, ++a)
    { 
      (*i).coord = a;
      (*i).adjoiningCoords = new unsigned[theAdjoiningCoordSize];
      //(*i).initAdjoins = NULL;
      (*i).diffuseSize = theAdjoiningCoordSize;
      (*i).adjoiningSize = theAdjoiningCoordSize;
      (*i).point = NULL;
      unsigned aCol(a/(theRowSize*theLayerSize)); 
      unsigned aLayer((a%(theRowSize*theLayerSize))/theRowSize); 
      unsigned aRow((a%(theRowSize*theLayerSize))%theRowSize); 
      if(aRootComp->geometry == CUBOID || isInsideCoord(a, aRootComp, 0))
        {
          //By default, the voxel is vacant and we set it to the root id:
          (*i).idx = rootID*theStride;
          for(unsigned j(0); j != theAdjoiningCoordSize; ++j)
            { 
              // By default let the adjoining voxel pointer point to the 
              // source voxel (i.e., itself)
              (*i).adjoiningCoords[j] = a;
            } 
          concatenateVoxel(*i, aRow, aLayer, aCol);
        }
      else
        {
          //We set id = theNullID if it is an invalid voxel, i.e., no molecules
          //will occupy it:
          (*i).idx = theNullID*theStride;
          //Concatenate some of the null voxels close to the surface:
          if(isInsideCoord(a, aRootComp, 4))
            {
              concatenateVoxel(*i, aRow, aLayer, aCol);
            }
        }
    }
  if(aRootComp->geometry == CUBOID || aRootComp->geometry == OPEN_CYLINDER)
    {
      concatenatePeriodicSurfaces();
    }
}

void SpatiocyteStepper::setPeriodicEdge()
{
  isPeriodicEdge = true;
}

bool SpatiocyteStepper::isPeriodicEdgeCoord(unsigned aCoord, Comp* aComp)
{
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aCoord, aRow, aLayer, aCol);
  if(aComp->system->getSuperSystem()->isRootSystem() &&
     ((aRow <= 1 && aCol <= 1) ||
      (aRow <= 1 && aLayer <= 1) ||
      (aCol <= 1 && aLayer <= 1) ||
      (aRow <= 1 && aCol >= theColSize-2) ||
      (aRow <= 1 && aLayer >= theLayerSize-2) ||
      (aCol <= 1 && aRow >= theRowSize-2) ||
      (aCol <= 1 && aLayer >= theLayerSize-2) ||
      (aLayer <= 1 && aRow >= theRowSize-2) ||
      (aLayer <= 1 && aCol >= theColSize-2) ||
      (aRow >= theRowSize-2 && aCol >= theColSize-2) ||
      (aRow >= theRowSize-2 && aLayer >= theLayerSize-2) ||
      (aCol >= theColSize-2 && aLayer >= theLayerSize-2)))
    {
      return true;
    }
  return false;
}

bool SpatiocyteStepper::isRemovableEdgeCoord(unsigned aCoord, Comp* aComp)
{
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aCoord, aRow, aLayer, aCol);
  int sharedCnt(0);
  int removeCnt(0);
  //x => col
  //y => layer
  //z => row
  if(aRow >= aComp->maxCoord.row)
    {
      ++sharedCnt;
      if(aComp->xyPlane == REMOVE_UPPER || aComp->xyPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aRow <= aComp->minCoord.row)
    {
      ++sharedCnt;
      if(aComp->xyPlane == REMOVE_LOWER || aComp->xyPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aLayer >= aComp->maxCoord.layer)
    {
      ++sharedCnt;
      if(aComp->xzPlane == REMOVE_UPPER || aComp->xzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aLayer <= aComp->minCoord.layer)
    {
      ++sharedCnt;
      if(aComp->xzPlane == REMOVE_LOWER || aComp->xzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aCol >= aComp->maxCoord.col)
    {
      ++sharedCnt;
      if(aComp->yzPlane == REMOVE_UPPER || aComp->yzPlane == REMOVE_BOTH)
        {
          ++removeCnt;
        }
    }
  if(aCol <= aComp->minCoord.col) 
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

void SpatiocyteStepper::concatenateVoxel(Voxel& aVoxel,
                                         unsigned aRow,
                                         unsigned aLayer,
                                         unsigned aCol)
{
  unsigned aCoord(aRow+
                      theRowSize*aLayer+
                      theRowSize*theLayerSize*aCol);
  if(aRow > 0)
    { 
      concatenateRows(aVoxel, aCoord, aRow-1, aLayer, aCol); 
    } 
  if(aLayer > 0)
    {
      concatenateLayers(aVoxel, aCoord, aRow, aLayer-1, aCol); 
    }
  if(aCol > 0)
    { 
      concatenateCols(aVoxel, aCoord, aRow, aLayer, aCol-1); 
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
    case OPEN_CYLINDER:
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
  aComp->nLength = aComp->lengthX;
  aComp->nWidth = aComp->lengthY;
  aComp->nHeight = aComp->lengthZ;
  rotate(aComp, aComp->lengthVector, -1);
  rotate(aComp, aComp->widthVector, -1);
  rotate(aComp, aComp->heightVector, -1);
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

Point SpatiocyteStepper::coord2point(unsigned aCoord)
{
  unsigned aGlobalRow;
  unsigned aGlobalLayer;
  unsigned aGlobalCol;
  coord2global(aCoord, aGlobalRow, aGlobalLayer, aGlobalCol);
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
  /*std::cout<<"inrow "<<aGlobalRow<<std::endl;
  std::cout<<"inlayer "<<aGlobalLayer<<std::endl;
  std::cout<<"incol "<<aGlobalCol<<std::endl;
  std::cout<<"theRowsize "<<theRowSize<<std::endl;
  std::cout<<"theLayersize "<<theLayerSize<<std::endl;
  std::cout<<"returnvalue "<<aGlobalRow+theRowSize*aGlobalLayer+theRowSize*theLayerSize*aGlobalCol<<std::endl;*/
  return aGlobalRow+theRowSize*aGlobalLayer+theRowSize*theLayerSize*aGlobalCol;
}

void SpatiocyteStepper::point2global(Point aPoint, 
                                     unsigned& aGlobalRow,
                                     unsigned& aGlobalLayer,
                                     unsigned& aGlobalCol)
{
  double row(0);
  double layer(0);
  double col(0);
      //std::cout<<"inSS1"<<std::endl;
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
  //std::cout<<"inrow "<<aGlobalRow<<std::endl;
  //std::cout<<"inlayer "<<aGlobalLayer<<std::endl;
 // std::cout<<"incol "<<aGlobalCol<<std::endl;
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
  //std::cout<<"row "<<row<<std::endl;
  //std::cout<<"layer "<<layer<<std::endl;
  //std::cout<<"col "<<col<<std::endl;
  aGlobalRow = (unsigned)row;
  aGlobalLayer = (unsigned)layer;
  aGlobalCol = (unsigned)col;
  //std::cout<<"inSS3"<<std::endl;
  if(aGlobalCol >= theColSize)
    {
      aGlobalCol = theColSize-1;
    }
  if(aGlobalRow >= theRowSize)
    {
      aGlobalRow = theRowSize-1;
    }
  if(aGlobalLayer >= theLayerSize)
    {
      aGlobalLayer = theLayerSize-1;
    }
  //std::cout<<"inSS4"<<std::endl;
}

void SpatiocyteStepper::coord2global(unsigned aCoord,
                                     unsigned& aGlobalRow,
                                     unsigned& aGlobalLayer,
                                     unsigned& aGlobalCol) 
{
  aGlobalCol = (aCoord)/(theRowSize*theLayerSize);
  aGlobalLayer = ((aCoord)%(theRowSize*theLayerSize))/theRowSize;
  aGlobalRow = ((aCoord)%(theRowSize*theLayerSize))%theRowSize;
}

void SpatiocyteStepper::concatenateLayers(Voxel& aVoxel,
                                          unsigned a,
                                          unsigned aRow,
                                          unsigned aLayer,
                                          unsigned aCol)
{
  unsigned b(aRow+
                 theRowSize*aLayer+
                 theRowSize*theLayerSize*aCol);
  Voxel& anAdjoiningVoxel(theLattice[b]);
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      if((aLayer+1)%2+(aCol)%2 == 1)
        {
          aVoxel.adjoiningCoords[VENTRALN] = b;
          anAdjoiningVoxel.adjoiningCoords[DORSALS] = a;
          if(aRow < theRowSize-1)
            {
              unsigned c(aRow+1+ 
                             theRowSize*aLayer+
                             theRowSize*theLayerSize*aCol);
              Voxel& anAdjoiningVoxel(theLattice[c]);
              aVoxel.adjoiningCoords[VENTRALS] = c;
              anAdjoiningVoxel.adjoiningCoords[DORSALN] = a;
            }
        }
      else
        {
          aVoxel.adjoiningCoords[VENTRALS] = b;
          anAdjoiningVoxel.adjoiningCoords[DORSALN] = a;
          if(aRow > 0)
            {
              unsigned c(aRow-1+ 
                             theRowSize*aLayer+
                             theRowSize*theLayerSize*aCol);
              Voxel& anAdjoiningVoxel(theLattice[c]);
              aVoxel.adjoiningCoords[VENTRALN] = c;
              anAdjoiningVoxel.adjoiningCoords[DORSALS] = a;
            }
        }
      break;
    case CUBIC_LATTICE: 
      aVoxel.adjoiningCoords[VENTRAL] = b;
      anAdjoiningVoxel.adjoiningCoords[DORSAL] = a;
      break;
    }
}


void SpatiocyteStepper::concatenateCols(Voxel& aVoxel,
                                        unsigned a,
                                        unsigned aRow,
                                        unsigned aLayer,
                                        unsigned aCol)
{
  unsigned b(aRow+
                 theRowSize*aLayer+
                 theRowSize*theLayerSize*aCol);
  Voxel& anAdjoiningVoxel(theLattice[b]);
  switch(LatticeType)
    {
    case HCP_LATTICE: 
      if(aLayer%2 == 0)
        {
          if((aCol+1)%2 == 1)
            {
              aVoxel.adjoiningCoords[NW] = b;
              anAdjoiningVoxel.adjoiningCoords[SE] = a;
              if(aRow < theRowSize - 1)
                {
                  unsigned c(aRow+1+ 
                                 theRowSize*aLayer+
                                 theRowSize*theLayerSize*aCol);
                  Voxel& anAdjoiningVoxel(theLattice[c]);
                  aVoxel.adjoiningCoords[SW] = c;
                  anAdjoiningVoxel.adjoiningCoords[NE] = a;
                }
              if(aLayer < theLayerSize-1)
                {
                  unsigned c(aRow+ 
                                 theRowSize*(aLayer+1)+
                                 theRowSize*theLayerSize*aCol);
                  Voxel& anAdjoiningVoxel(theLattice[c]);
                  aVoxel.adjoiningCoords[WEST] = c;
                  anAdjoiningVoxel.adjoiningCoords[EAST] = a;
                }
            }
          else
            {
              aVoxel.adjoiningCoords[SW] = b;
              anAdjoiningVoxel.adjoiningCoords[NE] = a;
              if(aRow > 0)
                {
                  unsigned c(aRow-1+ 
                                 theRowSize*aLayer+
                                 theRowSize*theLayerSize*aCol);
                  Voxel& anAdjoiningVoxel(theLattice[c]);
                  aVoxel.adjoiningCoords[NW] = c;
                  anAdjoiningVoxel.adjoiningCoords[SE] = a;
                }
              if(aLayer > 0)
                {
                  unsigned c(aRow+ 
                                 theRowSize*(aLayer-1)+
                                 theRowSize*theLayerSize*aCol);
                  Voxel& anAdjoiningVoxel(theLattice[c]);
                  aVoxel.adjoiningCoords[WEST] = c;
                  anAdjoiningVoxel.adjoiningCoords[EAST] = a;
                }
            }
        }
      else
        {
          if((aCol+1)%2 == 1)
            {
              aVoxel.adjoiningCoords[SW] = b;
              anAdjoiningVoxel.adjoiningCoords[NE] = a;
              if(aRow > 0)
                {
                  unsigned c(aRow-1+ 
                                 theRowSize*aLayer+
                                 theRowSize*theLayerSize*aCol);
                  Voxel& anAdjoiningVoxel(theLattice[c]);
                  aVoxel.adjoiningCoords[NW] = c;
                  anAdjoiningVoxel.adjoiningCoords[SE] = a;
                }
              if(aLayer < theLayerSize-1)
                {
                  unsigned c(aRow+ 
                                 theRowSize*(aLayer+1)+
                                 theRowSize*theLayerSize*aCol);
                  Voxel& anAdjoiningVoxel(theLattice[c]);
                  aVoxel.adjoiningCoords[WEST] = c;
                  anAdjoiningVoxel.adjoiningCoords[EAST] = a;
                }
            }
          else
            {
              aVoxel.adjoiningCoords[NW] = b;
              anAdjoiningVoxel.adjoiningCoords[SE] = a;
              if(aRow < theRowSize - 1)
                {
                  unsigned c(aRow+1+ 
                                 theRowSize*aLayer+
                                 theRowSize*theLayerSize*aCol);
                  Voxel& anAdjoiningVoxel(theLattice[c]);
                  aVoxel.adjoiningCoords[SW] = c;
                  anAdjoiningVoxel.adjoiningCoords[NE] = a;
                }
              if(aLayer > 0)
                {
                  unsigned c(aRow+
                                 theRowSize*(aLayer-1)+
                                 theRowSize*theLayerSize*aCol);
                  Voxel& anAdjoiningVoxel(theLattice[c]);
                  aVoxel.adjoiningCoords[WEST] = c;
                  anAdjoiningVoxel.adjoiningCoords[EAST] = a;
                }
            }
        }
      break;
    case CUBIC_LATTICE:
      aVoxel.adjoiningCoords[WEST] = b;
      anAdjoiningVoxel.adjoiningCoords[EAST] = a;
      break;
    }
}

void SpatiocyteStepper::concatenateRows(Voxel& aVoxel,
                                        unsigned a,
                                        unsigned aRow,
                                        unsigned aLayer,
                                        unsigned aCol)
{
  unsigned b(aRow+ 
                 theRowSize*aLayer+ 
                 theRowSize*theLayerSize*aCol);
  Voxel& anAdjoiningVoxel(theLattice[b]);
  aVoxel.adjoiningCoords[NORTH] = b;
  anAdjoiningVoxel.adjoiningCoords[SOUTH] = a;
}


void SpatiocyteStepper::concatenatePeriodicSurfaces()
{
  Comp* aRootComp(theComps[0]);
  for(unsigned i(0); i<=theRowSize*theLayerSize*theColSize-theRowSize;
      i+=theRowSize)
    {
      unsigned j(i+theRowSize-1);
      unsigned srcCoord(coord2row(i)+ 
                            theRowSize*coord2layer(i)+ 
                            theRowSize*theLayerSize*coord2col(i)); 
      unsigned destCoord(coord2row(j)+
                             theRowSize*coord2layer(i)+
                             theRowSize*theLayerSize*coord2col(i)); 
      if(aRootComp->xyPlane == UNIPERIODIC)
        { 
          replaceUniVoxel(srcCoord, destCoord);
        }
      else if(aRootComp->xyPlane == PERIODIC)
        { 
          replaceVoxel(srcCoord, destCoord);
        }
      else if(!isPeriodicEdge)
        {
          //We cannot have valid voxels pointing to itself it it is not periodic
          //to avoid incorrect homodimerization reaction. So we set such
          //molecules to null ID.
          theLattice[srcCoord].idx = theNullID*theStride;
          theLattice[destCoord].idx = theNullID*theStride;
        }
    }
  for(unsigned i(0); i<=theRowSize*theLayerSize*(theColSize-1)+theRowSize;)
    {
      unsigned j(theRowSize*(theLayerSize-1)+i);
      unsigned srcCoord(coord2row(i)+
                            theRowSize*coord2layer(i)+ 
                            theRowSize*theLayerSize*coord2col(i)); 
      unsigned destCoord(coord2row(i)+
                             theRowSize*coord2layer(j)+ 
                             theRowSize*theLayerSize*coord2col(i)); 
      if(aRootComp->xzPlane == UNIPERIODIC)
        { 
          replaceUniVoxel(srcCoord, destCoord);
        }
      else if(aRootComp->xzPlane == PERIODIC)
        { 
          replaceVoxel(srcCoord, destCoord);
        }
      else if(!isPeriodicEdge)
        {
          theLattice[srcCoord].idx = theNullID*theStride;
          theLattice[destCoord].idx = theNullID*theStride;
        }
      ++i;
      if(coord2layer(i) != 0)
        {
          i += (theLayerSize-1)*theRowSize;
        }
    }
  for(unsigned i(0); i!=theRowSize*theLayerSize; ++i)
    {
      unsigned srcCoord(coord2row(i)+
                            theRowSize*coord2layer(i)+ 
                            theRowSize*theLayerSize*0); 
      unsigned destCoord(coord2row(i)+
                             theRowSize*coord2layer(i)+ 
                             theRowSize*theLayerSize*(theColSize-1)); 
      if(aRootComp->yzPlane == UNIPERIODIC)
        { 
          replaceUniVoxel(srcCoord, destCoord);
        }
      else if(aRootComp->yzPlane == PERIODIC)
        { 
          replaceVoxel(srcCoord, destCoord);
        }
      else if(!isPeriodicEdge)
        {
          theLattice[srcCoord].idx = theNullID*theStride;
          theLattice[destCoord].idx = theNullID*theStride;
        }
    }
}

unsigned SpatiocyteStepper::coord2row(unsigned aCoord)
{
  return (aCoord%(theRowSize*theLayerSize))%theRowSize;
}

unsigned SpatiocyteStepper::coord2layer(unsigned aCoord)
{
  return (aCoord%(theRowSize*theLayerSize))/theRowSize;
}

unsigned SpatiocyteStepper::coord2col(unsigned aCoord)
{
  return aCoord/(theRowSize*theLayerSize);
}

void SpatiocyteStepper::replaceVoxel(unsigned src, unsigned dest)
{
  Voxel& aSrcVoxel(theLattice[src]);
  Voxel& aDestVoxel(theLattice[dest]);
  if(aSrcVoxel.idx != theNullID*theStride && 
     aDestVoxel.idx != theNullID*theStride)
    {
      for(unsigned j(0); j!=theAdjoiningCoordSize; ++j)
        {
          if(aSrcVoxel.adjoiningCoords[j] == src &&
              aDestVoxel.adjoiningCoords[j] != dest)
            {
              aSrcVoxel.adjoiningCoords[j] = aDestVoxel.adjoiningCoords[j];
              for(unsigned k(0); k!=theAdjoiningCoordSize; ++k)
                {
                  Voxel& aVoxel(theLattice[aDestVoxel.adjoiningCoords[j]]);
                  if(aVoxel.adjoiningCoords[k] == dest)
                    {
                      aVoxel.adjoiningCoords[k] = src;
                    }
                }
            }
        }
      aDestVoxel.idx = theNullID*theStride;
    }
}

void SpatiocyteStepper::replaceUniVoxel(unsigned src, unsigned dest)
{
  Voxel& aSrcVoxel(theLattice[src]);
  Voxel& aDestVoxel(theLattice[dest]);
  if(aSrcVoxel.idx*theStride != theNullID &&
     aDestVoxel.idx != theNullID*theStride)
    {
      for(unsigned j(0); j!=theAdjoiningCoordSize; ++j)
        {
          if(aSrcVoxel.adjoiningCoords[j] == src)
            {
              for(unsigned k(0); k!=theAdjoiningCoordSize; ++k)
                {
                  Voxel& aVoxel(theLattice[aDestVoxel.adjoiningCoords[j]]);
                  if(aVoxel.adjoiningCoords[k] == dest)
                    {
                      aVoxel.adjoiningCoords[k] = src;
                    }
                }
            }
        }
      aDestVoxel.idx = theNullID*theStride;
    }
}

void SpatiocyteStepper::shuffleAdjoiningCoords()
{
  for(std::vector<Voxel>::iterator i(theLattice.begin());
      i != theLattice.end(); ++i)
    {
      if((*i).idx != theNullID*theStride)
        { 
          gsl_ran_shuffle(getRng(), (*i).adjoiningCoords, theAdjoiningCoordSize,
                          sizeof(unsigned));
        }
    }
}

void SpatiocyteStepper::setCompVoxelProperties()
{
  for(std::vector<Comp*>::iterator i(theComps.begin());
      i != theComps.end(); ++i)
    {
      if((*i)->isProcessComp)
        {
          continue;
        }
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
  if(aComp->geometry == CUBOID)
    {
      if(isPeriodicEdge)
        {
          aComp->vacantSpecies->removePeriodicEdgeVoxels();
        }
      aComp->vacantSpecies->removeSurfacePlanes();
    }
  setDiffusiveComp(aComp);
}

void SpatiocyteStepper::setVolumeCompProperties(Comp* aComp)
{
  setDiffusiveComp(aComp);
}

void SpatiocyteStepper::setDiffusiveComp(Comp* aComp)
{
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
          unsigned aCoord(aVacantSpecies->getCoord(i));
          aComp->diffusiveComp->vacantSpecies->addCompVoxel(aCoord);
        }
      aVacantSpecies->clearMolecules();
    }
}

void SpatiocyteStepper::setSurfaceVoxelProperties(Comp* aComp)
{
  if(!aComp->diffusiveComp)
    {
      Comp* aSuperComp(system2Comp(aComp->system->getSuperSystem())); 
      const unsigned intID(aSuperComp->vacantSpecies->getID());
      Comp* aSuperSuperComp(system2Comp(aSuperComp->system->getSuperSystem())); 
      unsigned extID(aSuperSuperComp->vacantSpecies->getID());
      if(aSuperSuperComp == aSuperComp)
        {
          extID = theNullID;
        }
      Species* aVacantSpecies(aComp->vacantSpecies);
      const unsigned surfaceID(aVacantSpecies->getID());
      //std::cout << "before:" << aVacantSpecies->size() << std::endl;
      for(unsigned i(0); i != aVacantSpecies->size(); )
        {
          const unsigned aCoord(aVacantSpecies->getCoord(i));
          if(setSurfaceDiffuseSize(aCoord, surfaceID, intID, extID,
                                   aSuperComp->system->isRootSystem(), i))
            {
              ++i;
            }
        }
      /*
      std::cout << "after:" << aVacantSpecies->size() << std::endl;
      for(unsigned i(0); i != aVacantSpecies->size(); )
        {
          const unsigned aCoord(aVacantSpecies->getCoord(i));
          if(setSurfaceDiffuseSize(aCoord, surfaceID, intID, extID,
                                   aSuperComp->system->isRootSystem(), i))
            {
              ++i;
            }
        }
      std::cout << "final:" << aVacantSpecies->size() << std::endl;
      */
      if(RemoveSurfaceBias)
        {
          for(unsigned i(0); i != aVacantSpecies->size(); ++i)
            {
              unsigned aCoord(aVacantSpecies->getCoord(i));
              addSurfaceAdjoins(aCoord, aComp);
            }
        }
    }
}


bool SpatiocyteStepper::setSurfaceDiffuseSize(const unsigned aCoord,
                                              const unsigned surfaceID,
                                              const unsigned intID,
                                              const unsigned extID,
                                              const bool isRemoveOrphan,
                                              const unsigned anIndex)
{
  Voxel& aVoxel(theLattice[aCoord]);
  unsigned* adjoins(aVoxel.adjoiningCoords);
  std::vector<unsigned> adjoinsCopy;
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      adjoinsCopy.push_back(adjoins[i]);
    }
  //Separate adjoining surface voxels and adjoining volume voxels.
  //Put the adjoining surface voxels at the beginning of the
  //adjoiningCoords list while the volume voxels are put at the end:
  unsigned f(0);
  unsigned b(theAdjoiningCoordSize-1);
  unsigned intCnt(0);
  unsigned extCnt(0);
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      const unsigned adjCoord(adjoinsCopy[i]);
      const unsigned anID(getID(theLattice[adjCoord]));
      if(anID == surfaceID)
        {
          adjoins[f++] = adjCoord; 
        }
      else
        {
          adjoins[b--] = adjCoord;
          if(anID == intID)
            {
              ++intCnt;
            }
          else if(anID == extID)
            {
              ++extCnt;
            }
        }
    }
  //Some surface voxels are not connected to the volume voxels. Remove these
  //orphaned surface voxels and update its neighbours.
  if(isRemoveOrphan && (!intCnt || !extCnt))
    {
      theSpecies[surfaceID]->removeCompVoxel(anIndex);
      if(!intCnt)
        {
          if(extID == theNullID)
            {
              theLattice[aCoord].idx = theNullID*theStride;
            }
          else
            {
              theSpecies[extID]->addCompVoxel(aCoord);
            }
        }
      else
        {
          theSpecies[intID]->addCompVoxel(aCoord);
        }
      for(unsigned i(0); i != f; ++i)
        {
          const unsigned adjCoord(adjoins[i]);
          setSurfaceDiffuseSize(adjCoord, surfaceID, intID, extID, false, 0);
        }
      return false;
    }
  aVoxel.diffuseSize = f;
  return true;
}

void SpatiocyteStepper::addSurfaceAdjoins(const unsigned aCoord,
                                          const Comp* aComp)
{
  Voxel& aVoxel(theLattice[aCoord]);
  if(aVoxel.diffuseSize > 5)
    {
      return;
    }
  std::vector<unsigned> extCoords;
  std::vector<double> extDists;
  const Point aPoint(coord2point(aCoord));
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      Voxel& adjoin(theLattice[aVoxel.adjoiningCoords[i]]);
      for(unsigned j(0); j != theAdjoiningCoordSize; ++j)
        { 
          const unsigned coord(adjoin.adjoiningCoords[j]);
          Voxel& extAdjoin(theLattice[coord]);
          if(extAdjoin.diffuseSize < 6 && coord != aCoord && 
             extAdjoin.idx != theNullID*theStride &&
             id2Comp(getID(extAdjoin))->dimension <= aComp->dimension)
            {
              if(std::find(extCoords.begin(), extCoords.end(),
                           coord) == extCoords.end())
                {
                  bool isAdjoined(false);
                  for(unsigned k(0); k != aVoxel.diffuseSize; ++k)
                    { 
                      if(aVoxel.adjoiningCoords[k] == coord)
                        {
                          isAdjoined = true;
                          break;
                        }
                    }
                  if(!isAdjoined)
                    {
                      extCoords.push_back(coord);
                      const Point extPoint(coord2point(coord));
                      extDists.push_back(distance(aPoint, extPoint));
                    }
                }
            }
        }
    }
  const unsigned addSize(6-aVoxel.diffuseSize);
  std::vector<double> tmpDists(extDists);
  std::sort(extDists.begin(), extDists.end());
  /*
  std::cout << "addSize:" << addSize << std::endl;
  for(unsigned j(0); j != extCoords.size(); ++j)
    {
      std::cout << "dist:" << extDists[j] << std::endl;
    }
    */
  for(unsigned i(0); i != addSize; ++i)
    {
      for(unsigned j(0); j != extCoords.size(); ++j)
        {
          if(tmpDists[j] == extDists[i])
            {
              Voxel& extAdjoin(theLattice[extCoords[j]]);
              aVoxel.adjoiningCoords[aVoxel.diffuseSize++] = extCoords[j];
              extAdjoin.adjoiningCoords[extAdjoin.diffuseSize++] = aCoord;
              extCoords[j] = extCoords.back();
              extCoords.pop_back();
              tmpDists[j] = tmpDists.back();
              tmpDists.pop_back();
              break;
            }
        }
    }
  /*
  if(aVoxel.diffuseSize != 6)
    {
      std::cout << "diffuse size is larger:" << aVoxel.diffuseSize << std::endl;
    }
    */
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
      if((*i)->isProcessComp)
        {
          continue;
        }
      //Only proceed if the volume is not enclosed:
      if(!(*i)->system->isRootSystem() && (*i)->dimension == 3)
        {
          Point a((*i)->centerPoint);
          Point u;
          u.x = (*i)->lengthX/2;
          u.y = (*i)->lengthY/2;
          u.z = (*i)->lengthZ/2;
          rotateX((*i)->rotateX, &u);
          rotateY((*i)->rotateY, &u);
          rotateZ((*i)->rotateZ, &u);
          for(std::vector<Comp*>::reverse_iterator j(theComps.rbegin());
              j != theComps.rend(); ++j)
            {
              if((*j)->isProcessComp)
                {
                  continue;
                }
              //Only proceed if the volume of the peer is enclosed:
              if((*i) != (*j) && !(*j)->system->isRootSystem() &&
                 (*j)->dimension == 3)
                {
                  //If i and j are peer volume compartments:
                  if((*i)->system->getSuperSystem() == 
                     (*j)->system->getSuperSystem())
                    {
                      Point b((*j)->centerPoint);
                      Point v;
                      v.x = (*j)->lengthX/2;
                      v.y = (*j)->lengthY/2;
                      v.z = (*j)->lengthZ/2;
                      rotateX((*j)->rotateX, &v);
                      rotateY((*j)->rotateY, &v);
                      rotateZ((*j)->rotateZ, &v);
                      if(std::min(a.x+u.x, a.x-u.x) <
                         std::max(b.x+v.x, b.x-v.x) &&
                         std::max(a.x+u.x, a.x-u.x) >
                         std::min(b.x+v.x, b.x-v.x) &&
                         std::min(a.y+u.y, a.y-u.y) <
                         std::max(b.y+v.y, b.y-v.y) &&
                         std::max(a.y+u.y, a.y-u.y) >
                         std::min(b.y+v.y, b.y-v.y) &&
                         std::min(a.z+u.z, a.z-u.z) <
                         std::max(b.z+v.z, b.z-v.z) &&
                         std::max(a.z+u.z, a.z-u.z) >
                         std::min(b.z+v.z, b.z-v.z))
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
      if((*i)->isProcessComp)
        {
          continue;
        }
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
  for(unsigned i(0); i != theLattice.size(); ++i)
    {
      if(theLattice[i].idx != theNullID*theStride)
        { 
          compartmentalizeVoxel(i, theComps[0]);
        }
    }
}

bool SpatiocyteStepper::compartmentalizeVoxel(unsigned aCoord, Comp* aComp)
{
  Voxel& aVoxel(theLattice[aCoord]);
  if(aComp->dimension == 3)
    {
      if(aComp->system->isRootSystem() ||
         isInsideCoord(aCoord, aComp, 0))
        {
          //Check if the voxel belongs to an intersecting peer compartment:
          if(!aComp->intersectPeers.empty())
            {
              //The function isPeerCoord also checks if the voxel is a future
              //surface voxel of the peer compartment: 
              if(isPeerCoord(aCoord, aComp))
                { 
                  return false;
                }
              if(aComp->surfaceSub && aComp->surfaceSub->enclosed)
                {
                  //Check if the voxel is neighbor of a peer voxel (either
                  //a future surface voxel)
                  if(isEnclosedSurfaceVoxel(aVoxel, aCoord, aComp))
                    {
                      aComp->surfaceSub->vacantSpecies->addCompVoxel(aCoord);
                      setMinMaxSurfaceDimensions(aCoord, aComp);
                      return true;
                    }
                }
            }
          if(aComp->isIntersectRoot)
            {
              Comp* aRootComp(system2Comp(aComp->system->getSuperSystem())); 
              if(aRootComp->surfaceSub && 
                 aRootComp->surfaceSub->enclosed && 
                 isRootSurfaceVoxel(aVoxel, aCoord, aRootComp))
                {
                  return false;
                }
              if(aComp->surfaceSub && 
                 isEnclosedRootSurfaceVoxel(aVoxel, aCoord, aComp, aRootComp))
                {
                  aComp->surfaceSub->vacantSpecies->addCompVoxel(aCoord);
                  setMinMaxSurfaceDimensions(aCoord, aComp);
                  return true;
                }
            }
          else if(aComp->isIntersectParent)
            {
              Comp* aParentComp(system2Comp(aComp->system->getSuperSystem())); 
              if(aComp->surfaceSub && aComp->surfaceSub->enclosed &&
                 isParentSurfaceVoxel(aVoxel, aCoord, aParentComp))
                {
                  aComp->surfaceSub->vacantSpecies->addCompVoxel(aCoord);
                  setMinMaxSurfaceDimensions(aCoord, aComp);
                  return true;
                }
            }
          for(unsigned i(0); i != aComp->immediateSubs.size(); ++i)
            {
              if(compartmentalizeVoxel(aCoord, aComp->immediateSubs[i]))
                {
                  setMinMaxSurfaceDimensions(aCoord, aComp);
                  return true;
                }
            }
          aComp->vacantSpecies->addCompVoxel(aCoord);
          setMinMaxSurfaceDimensions(aCoord, aComp);
          return true;
        }
      if(aComp->surfaceSub)
        { 
          if(isInsideCoord(aCoord, aComp, 4) &&
             isSurfaceVoxel(aVoxel, aCoord, aComp))
            {
              aComp->surfaceSub->vacantSpecies->addCompVoxel(aCoord);
              setMinMaxSurfaceDimensions(aCoord, aComp);
              return true;
            }
        }
    }
  //aComp dimension == 2 or 1
  else if(aComp->system->getSuperSystem()->isRootSystem())
    {
      Comp* aRootComp(system2Comp(aComp->system->getSuperSystem())); 
      if(!isInsideCoord(aCoord, aRootComp, -4) &&
         isRootSurfaceVoxel(aVoxel, aCoord, aRootComp))
        {
          aComp->vacantSpecies->addCompVoxel(aCoord);
          setMinMaxSurfaceDimensions(aCoord, aRootComp);
          return true;
        }
    }
  return false;
}

/*
 * This causes surface voxel bug in amoeba spiral
bool SpatiocyteStepper::isRootSurfaceVoxel(Voxel& aVoxel, unsigned aCoord,
                                           Comp* aComp)
{
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      const unsigned adjoinCoord(aVoxel.adjoiningCoords[i]);
      Voxel& adjoin(theLattice[adjoinCoord]);
      if((adjoin.idx == theNullID*theStride || adjoinCoord == aCoord) &&
         !isInsideCoord(adjoinCoord, aComp, 0))
        {
          return true;
        }
    }
  return false;
}
*/

bool SpatiocyteStepper::isRootSurfaceVoxel(Voxel& aVoxel, unsigned aCoord,
                                           Comp* aComp)
{
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      const unsigned adjoinCoord(aVoxel.adjoiningCoords[i]);
      const Voxel& adjoin(theLattice[adjoinCoord]);
      if(adjoin.idx == theNullID*theStride || adjoinCoord == aCoord)
        {
          return true;
        }
    }
  return false;
}


bool SpatiocyteStepper::isParentSurfaceVoxel(Voxel& aVoxel, unsigned aCoord,
                                             Comp* aComp)
{
  if(!isInsideCoord(aCoord, aComp, -4))
    {
      for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
        {
          if(!isInsideCoord(aVoxel.adjoiningCoords[i], aComp, 0))
            {
              return true;
            }
        }
    }
  return false;
}

bool SpatiocyteStepper::isEnclosedRootSurfaceVoxel(Voxel& aVoxel, 
                                                   unsigned aCoord,
                                                   Comp* aComp,
                                                   Comp* aRootComp)
{ 
  if(aComp->surfaceSub->enclosed || aRootComp->enclosed <= aComp->enclosed)
    { 
      if(!isInsideCoord(aCoord, aRootComp, -4))
        {
          if(!aRootComp->surfaceSub || 
             (aRootComp->surfaceSub && !aRootComp->surfaceSub->enclosed))
            {
              if(isRootSurfaceVoxel(aVoxel, aCoord, aComp))
                {
                  return true;
                }
            }
          else
            {
              for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
                {
                  unsigned adjoinCoord(aVoxel.adjoiningCoords[i]);
                  Voxel& adjoin(theLattice[adjoinCoord]);
                  if(isRootSurfaceVoxel(adjoin, adjoinCoord, aComp))
                    {
                      return true;
                    }
                }
            }
        }
    }
  return false;
}

bool SpatiocyteStepper::isSurfaceVoxel(Voxel& aVoxel, unsigned aCoord,
                                       Comp* aComp)
{
  if(isPeerCoord(aCoord, aComp))
    { 
      return false;
    }
  if(aComp->surfaceSub && !aComp->surfaceSub->enclosed &&
     isLowerPeerCoord(aCoord, aComp))
    {
      return false;
    }
  if(aComp->isIntersectRoot)
    {
      Comp* aRootComp(system2Comp(aComp->system->getSuperSystem())); 
      if(aRootComp->surfaceSub && aRootComp->surfaceSub->enclosed && 
         isRootSurfaceVoxel(aVoxel, aCoord, aRootComp))
        {
          return false;
        }
    }
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      if(isInsideCoord(aVoxel.adjoiningCoords[i], aComp, 0))
        {
          return true;
        }
    }
  return false;
}

bool SpatiocyteStepper::isEnclosedSurfaceVoxel(Voxel& aVoxel, 
                                               unsigned aCoord,
                                               Comp* aComp)
{
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      if(isPeerCoord(aVoxel.adjoiningCoords[i], aComp))
        {
          return true;
        }
    }
  return false;
}

bool SpatiocyteStepper::isLowerPeerCoord(unsigned aCoord, Comp* aComp)
{
  for(std::vector<Comp*>::iterator i(aComp->intersectLowerPeers.begin());
      i != aComp->intersectLowerPeers.end(); ++i)
    {
      if(isInsideCoord(aCoord, *i, 0))
        {
          return true;
        }
    }
  return false;
}

bool SpatiocyteStepper::isPeerCoord(unsigned aCoord, Comp* aComp)
{
  for(std::vector<Comp*>::iterator i(aComp->intersectPeers.begin());
      i != aComp->intersectPeers.end(); ++i)
    {
      if(isInsideCoord(aCoord, *i, 0))
        {
          return true;
        }
    }
  //The voxel is not inside one of the intersecting peer compartments,
  //we are now checking if the one of voxel's neighbor belongs to the
  //peer comp, to determine if aVoxel is a surface voxel of the peer:
  for(std::vector<Comp*>::iterator i(aComp->intersectPeers.begin());
      i != aComp->intersectPeers.end(); ++i)
    {
      if((*i)->surfaceSub && (*i)->surfaceSub->enclosed)
        { 
          if(isInsideCoord(aCoord, *i, 4))
            {
              for(unsigned j(0); j != theAdjoiningCoordSize; ++j)
                {
                  if(isInsideCoord(theLattice[aCoord].adjoiningCoords[j],
                                   *i, 0))
                    {
                      return true;
                    }
                }
            }
        }
    }
  return false;
}

bool SpatiocyteStepper::isLineVoxel(Voxel& aVoxel, unsigned aCoord,
                                    Comp* aComp)
{
    const double safety(2.0);
    const Point aPoint(coord2point(aCoord));

    double dist(aPoint.x - aComp->centerPoint.x);
    if (-safety < dist && dist <= 0)
    {
        // This is not efficient because we don't need to check volume voxels.
        // However, at this time, aVoxel.adjoinigVoxels is not properly 
        // aligned yet.
        for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
        {
            const Point adjoiningPoint(coord2point(aVoxel.adjoiningCoords[i]));
            const double dist_i(adjoiningPoint.x - aComp->centerPoint.x);
            if (dist_i > 0)
            {
                return true;
            }
        }
    }
                       
    dist = aPoint.y - aComp->centerPoint.y;
    if (-safety < dist && dist <= 0)
    {
        // This is not efficient because we don't need to check volume voxels.
        // However, at this time, aVoxel.adjoinigVoxels is not properly 
        // aligned yet.
        for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
        {
            const Point adjoiningPoint(coord2point(aVoxel.adjoiningCoords[i]));
            const double dist_i(adjoiningPoint.y - aComp->centerPoint.y);
            if (dist_i > 0)
            {
                return true;
            }
        }
    }
    return false;
}

unsigned SpatiocyteStepper::getStartCoord()
{
  return 0;
}


void SpatiocyteStepper::setMinMaxSurfaceDimensions(unsigned aCoord, 
                                                   Comp* aComp)
{
  //maxCoord, minCoord, maxPoint and minPoint are only valid for unrotated
  //cuboid compartments
  if(aComp->geometry != CUBOID || aComp->rotateX != 0 || aComp->rotateY != 0 || 
     aComp->rotateZ != 0)
    {
      return;
    }
  unsigned aRow;
  unsigned aLayer;
  unsigned aCol;
  coord2global(aCoord, aRow, aLayer, aCol);
  if(aRow < aComp->minCoord.row)
    {
      aComp->minCoord.row = aRow;
      if(aComp->surfaceSub)
        {
          aComp->surfaceSub->minCoord.row = aRow;
        }
    }
  else if(aRow > aComp->maxCoord.row)
    {
      aComp->maxCoord.row = aRow;
      if(aComp->surfaceSub)
        {
          aComp->surfaceSub->maxCoord.row = aRow;
        }
    }
  if(aCol < aComp->minCoord.col)
    {
      aComp->minCoord.col = aCol;
      if(aComp->surfaceSub)
        {
          aComp->surfaceSub->minCoord.col = aCol;
        }
    }
  else if(aCol > aComp->maxCoord.col)
    {
      aComp->maxCoord.col = aCol;
      if(aComp->surfaceSub)
        {
          aComp->surfaceSub->maxCoord.col = aCol;
        }
    }
  if(aLayer < aComp->minCoord.layer)
    {
      aComp->minCoord.layer = aLayer;
      if(aComp->surfaceSub)
        {
          aComp->surfaceSub->minCoord.layer = aLayer;
        }
    }
  else if(aLayer > aComp->maxCoord.layer)
    {
      aComp->maxCoord.layer = aLayer;
      if(aComp->surfaceSub)
        {
          aComp->surfaceSub->maxCoord.layer = aLayer;
        }
    }
  Point aPoint(coord2point(aCoord));
  if(aPoint.x < aComp->minPoint.x)
    {
      aComp->minPoint.x = aPoint.x;
    }
  else if(aPoint.x > aComp->maxPoint.x)
    {
      aComp->maxPoint.x = aPoint.x;
    }
  if(aPoint.y < aComp->minPoint.y)
    {
      aComp->minPoint.y = aPoint.y;
    }
  else if(aPoint.y > aComp->maxPoint.y)
    {
      aComp->maxPoint.y = aPoint.y;
    }
  if(aPoint.z < aComp->minPoint.z)
    {
      aComp->minPoint.z = aPoint.z;
    }
  else if(aPoint.z > aComp->maxPoint.z)
    {
      aComp->maxPoint.z = aPoint.z;
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

bool SpatiocyteStepper::isInsideCoord(unsigned aCoord,
                                      Comp* aComp, double delta)
{
  Point aPoint(coord2point(aCoord));
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
      if((aPoint.x >= aWestPoint.x-delta && aPoint.x <= anEastPoint.x+delta &&
          distance(aPoint, aCenterPoint) <= aRadius+delta))
        { 
          return true;
        }
      break;
    case OPEN_CYLINDER: 
      //The axial point of the cylindrical portion of the rod:
      aCenterPoint.x = aPoint.x; 
      aRadius = aComp->lengthY/2+nVoxelRadius;
      //If the distance between the voxel and the center point is less than 
      //or equal to the radius, then the voxel must be inside the Comp:
      if(distance(aPoint, aCenterPoint) <= aRadius+delta)
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
      if((aPoint.x >= aWestPoint.x-delta && aPoint.x <= anEastPoint.x+delta &&
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
  //number of molecules that needs to be populated on a vacant species:
  std::vector<unsigned> vacantPopulations;
  std::vector<Species*> prioritySpecies;
  std::vector<Species*> multiscaleSpecies;
  std::vector<Species*> diffusiveSpecies;
  std::vector<Species*> normalSpecies;
  vacantPopulations.resize(theSpecies.size());
  for(unsigned i(0); i != theSpecies.size(); ++i)
    {
      vacantPopulations[i] = 0;
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
          vacantPopulations[(*i)->getVacantSpecies()->getID()] += 
            (*i)->getPopulateCoordSize();
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
  for(unsigned i(0); i != prioritySpecies.size(); ++i)
    {
      //First get the highest priority species and get its vacant species.
      //Then populate all the species that are supposed to be populated on
      //the vacant species according to priority:
      if(!prioritySpecies[i]->getIsPopulated())
        {
          unsigned aVacantID(prioritySpecies[i]->getVacantSpecies()->getID());
          if(vacantPopulations[aVacantID])
            {
              unsigned available(theSpecies[aVacantID]->getPopulatableSize());
              if(vacantPopulations[aVacantID] > available)
                {
                  THROW_EXCEPTION(ValueError, String(
                          getPropertyInterface().getClassName()) +
                          ": There are " +
                          int2str(vacantPopulations[aVacantID]) +
                          " total molecules that must be uniformly " +
                          "populated,\nbut there are only "
                          + int2str(available) + " vacant voxels of [" + 
                          theSpecies[aVacantID]->getIDString() +
                          "] that can be populated on.");
                } 
              if(double(vacantPopulations[aVacantID])/available > 0.2)
                { 
                  populateSpeciesDense(prioritySpecies, theSpecies[aVacantID],
                                       vacantPopulations[aVacantID], available);
                }
              else
                {
                  populateSpeciesSparse(prioritySpecies, theSpecies[aVacantID],
                                        available);
                }
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

unsigned SpatiocyteStepper::getID(const Voxel* aVoxel) const
{
  return aVoxel->idx/theStride;
}

unsigned SpatiocyteStepper::getID(const Voxel& aVoxel) const
{
  return aVoxel.idx/theStride;
}

void SpatiocyteStepper::populateSpeciesDense(std::vector<Species*>&
                                             aSpeciesList,
                                             Species* aVacantSpecies,
                                             unsigned aSize,
                                             unsigned availableVoxelSize)
{
  cout << "    Populating densely on " << aVacantSpecies->getIDString() <<
    ", available size:" << availableVoxelSize << std::endl;
  unsigned count(0);
  unsigned* populateVoxels(new unsigned[aSize]);
  unsigned* availableVoxels(new unsigned [availableVoxelSize]); 
  unsigned k(0);
  for(unsigned i(0); i != aVacantSpecies->size(); ++i)
    {
      if(aVacantSpecies->getIsCompVacant())
        {
          //Some of the voxels of the comp vacant species could have been
          //occupied by interface species, so we need to only include indices
          //that are still vacant:
          if(getID(aVacantSpecies->getMolecule(i)) == aVacantSpecies->getID())
            {
              availableVoxels[k++] = i;
            }
        }
      else
        {
          availableVoxels[i] = i;
        }
    }
  gsl_ran_choose(getRng(), populateVoxels, aSize, availableVoxels,
                 availableVoxelSize, sizeof(unsigned));
  //gsl_ran_choose arranges the position ascending, so we need
  //to shuffle the order of voxel positions:
  gsl_ran_shuffle(getRng(), populateVoxels, aSize, sizeof(unsigned)); 
  for(std::vector<Species*>::const_iterator i(aSpeciesList.begin());
      i != aSpeciesList.end(); ++i)
    {
      if((*i)->getVacantSpecies()->getID() == aVacantSpecies->getID())
        {
          (*i)->populateCompUniformDense(populateVoxels, &count);
        }
    }
  delete[] populateVoxels;
  delete[] availableVoxels;
}

void SpatiocyteStepper::populateSpeciesSparse(std::vector<Species*>&
                                              aSpeciesList,
                                              Species* aVacantSpecies,
                                              unsigned availableVoxelSize)
{
  cout << "    Populating sparsely on " << aVacantSpecies->getIDString() <<
    ", available size:" << availableVoxelSize << std::endl;
  for(std::vector<Species*>::const_iterator i(aSpeciesList.begin());
      i != aSpeciesList.end(); ++i)
    {
      if((*i)->getVacantSpecies()->getID() == aVacantSpecies->getID())
        {
          (*i)->populateCompUniformSparse();
        }
    }
}


void SpatiocyteStepper::clearComp(Comp* aComp)
{
  for(std::vector<Species*>::const_iterator i(aComp->species.begin());
      i != aComp->species.end(); ++i)
    {
      if(!(*i)->getIsCompVacant())
        {
          (*i)->removeMolecules();
          (*i)->updateMolecules();
        }
    }
  //A Comp may have more than one CompVacant species. Some CompVacant species
  //are isReactiveVacant or isDiffusiveVacant, so their numbers must be updated
  //after species occupying on them removed first (above):
  for(std::vector<Species*>::const_iterator i(aComp->species.begin());
      i != aComp->species.end(); ++i)
    {
      if((*i)->getIsCompVacant())
        {
          (*i)->updateMolecules();
        }
    }
}


std::vector<Comp*> const& SpatiocyteStepper::getComps() const
{
  return theComps;
}

