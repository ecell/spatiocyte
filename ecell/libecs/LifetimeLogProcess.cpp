//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 2006-2009 Keio University
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

#include <libecs/LifetimeLogProcess.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(LifetimeLogProcess, Process); 

void LifetimeLogProcess::initialize()
{
  if(isInitialized)
    {
      return;
    }
  SpatiocyteProcess::initialize();
  isPriorityQueued = true;
  theTotalIterations = Iterations;
}

void LifetimeLogProcess::initializeFirst()
{
  SpatiocyteProcess::initializeFirst();
  isBindingSiteReaction.resize(getStepper()->getProcessVector().size(), false);
  isAddDimerReaction.resize(getStepper()->getProcessVector().size(), false);
  isRemoveDimerReaction.resize(getStepper()->getProcessVector().size(), false);
  isDimerizationReaction.resize(getStepper()->getProcessVector().size(), false);
  isDedimerizationReaction.resize(getStepper()->getProcessVector().size(),
                                  false);
  isTrackedSpecies.resize(theSpecies.size(), false);
  trackedSpeciesSize.resize(theSpecies.size(), 0);
  isTrackedDimerSpecies.resize(theSpecies.size(), false);
  unsigned cntTracked(0);
  for(VariableReferenceVector::iterator
      i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      Species* aSpecies(theSpatiocyteStepper->variable2species(
                                     (*i).getVariable())); 
      isTrackedSpecies[aSpecies->getID()] = true;
      aSpecies->setIsTagged();
      ++cntTracked;
    }
  setTrackedDimerSpecies();
  //This process is only for tracking the lifetime of one particular species
  //that goes through several intermediate species before finally untracked
  //when it becomes an untracked species. For example if we want to track
  //the lifetime of a MinDatp_m molecule after it is bound to the membrane and
  //it goes through another state such as MinDadp_m before dissociating and
  //becoming a Vacant species:
  //trackedSpecies = MinDatp_m, MinDadp_m (coefficient -1)
  //untrackedSpecies = Vacant (coefficient 1)
  if(!cntTracked)
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + 
                      "]: Lifetime logging requires at least one " +
                      "nonHD variable as the tracked species, " +
                      "but none is given."); 
    }
}

void LifetimeLogProcess::setTrackedDimerSpecies()
{
  std::vector<Process*> aProcessVector(getStepper()->getProcessVector());
  for(std::vector<Process*>::const_iterator i(aProcessVector.begin());
      i != aProcessVector.end(); ++i)
    {
      const ReactionProcess* aProcess(dynamic_cast<ReactionProcess*>(*i));
      if(aProcess)
        {
          Species* A(aProcess->getA());
          Species* B(aProcess->getB());
          Species* C(aProcess->getC());
          Species* D(aProcess->getD());
          //Dimer -> Monomer + Monomer
          //See if a tracked substrate species is a dimer by
          //checking if the products are made up of two tracked species:
          if(((A && isTrackedSpecies[A->getID()] && !B) ||
              (B && isTrackedSpecies[B->getID()] && !A)) && C && D)
            {
              if(isTrackedSpecies[C->getID()] && isTrackedSpecies[D->getID()])
                {
                  if(A)
                    {
                      isTrackedDimerSpecies[A->getID()] = true;
                    }
                  else
                    {
                      isTrackedDimerSpecies[B->getID()] = true;
                    }
                  isDedimerizationReaction[aProcess->getID()] = true;
                }
            }
          //Monomer + Monomer -> Dimer
          //See if a tracked product species is a dimer by
          //checking if the two substrates are tracked species:
          else if(A && isTrackedSpecies[A->getID()] &&
                  B && isTrackedSpecies[B->getID()] && C && !D)
            {
              isTrackedDimerSpecies[C->getID()] = true;
              isDimerizationReaction[aProcess->getID()] = true;
            }
        }
    }
}

void LifetimeLogProcess::initializeSecond()
{
  availableTagIDs.resize(0);
  theTagTimes.resize(0);
}

void LifetimeLogProcess::initializeLastOnce()
{
  if(LogEnd == libecs::INF)
    {
      THROW_EXCEPTION(ValueError, String(
                  getPropertyInterface().getClassName()) +
                  "[" + getFullID().asString() + 
                  "]: LogEnd must be specified for " + 
                  "LifetimeLogProcess when LogInterval > 0."); 
    }
  if(LogEnd <= LogStart)
    {
      THROW_EXCEPTION(ValueError, String(
                  getPropertyInterface().getClassName()) +
                  "[" + getFullID().asString() + 
                  "]: LogEnd must be larger than LogStart " + 
                  "LifetimeLogProcess."); 
    }
  timePoints = (unsigned)ceil((LogEnd-LogStart)/theInterval)+1;
  theLogValues.resize(timePoints);
  unsigned aDataSize(1); //Only average diffusion now
  for(unsigned i(0); i != timePoints; ++i)
    {
      theLogValues[i].resize(aDataSize, 0);
      for(unsigned j(0); j != aDataSize; ++j)
        {
          theLogValues[i][j] = 0;
        }
    }
}

void LifetimeLogProcess::initializeFourth()
{
  reset();
}

void LifetimeLogProcess::reset()
{
  konCnt = 0;
  logCnt = 0;
  totalSize = 0;
  completedSquaredDisplacement = 0;
  totalDuration = 0;
  timePointCnt = 0;
  //The following is done after populating all molecules:
  for(unsigned i(0); i != isTrackedSpecies.size(); ++i)
    {
      if(isTrackedSpecies[i])
        {
          Species* aSpecies(theSpecies[i]);
          for(unsigned j(0); j != aSpecies->size(); ++j)
            {
              initTrackedMolecule(aSpecies, j);
            }
        }
    }
  for(unsigned i(0); i != isTrackedDimerSpecies.size(); ++i)
    {
      if(isTrackedDimerSpecies[i])
        {
          Species* aSpecies(theSpecies[i]);
          for(unsigned j(0); j != aSpecies->size(); ++j)
            {
              initTrackedDimer(aSpecies, j);
            }
        }
    }
}

void LifetimeLogProcess::initializeFifth()
{
  //Only fire at fixed intervals if LogInterval is given:
  if(LogInterval > 0)
    {
      theInterval = LogInterval;
    }
  if(!LogStart)
    {
      LogStart = 2*getStepper()->getMinStepInterval();
      theTime = getStepper()->getMinStepInterval();
    }
  else
    {
      theTime = std::max(LogStart-(theInterval/2), 1e-10);
    }
  thePriorityQueue->move(theQueueID);
}
void LifetimeLogProcess::fire()
{
  if(theTime < LogStart)
    {
      reset();
      doPreLog();
      if(LogInterval > 0)
        {
          theTime = LogStart;
        }
      else
        {
          theTime = LogEnd;
        }
      thePriorityQueue->moveTop();
      return;
    }
  if(theTime <= LogEnd)
    {
      logValues();
      ++timePointCnt;
    }
  if(theTotalIterations == 1)
    {
      saveATimePoint(theLogFile, theTime, 1, timePointCnt-1);
      if(theTime >= LogEnd)
        {
          theLogFile.close();
        }
    }
  else if(theTime >= LogEnd && Iterations > 0)
    {
      --Iterations;
      saveFile();
      if(Iterations)
        {
          theSpatiocyteStepper->reset(Iterations);
          return;
        }
      else
        {
          theInterval = libecs::INF;
        }
    }
  theTime += theInterval;
  thePriorityQueue->moveTop();
}

void LifetimeLogProcess::logValues()
{
  theLogValues[timePointCnt][0] += getAverageDiffusion();
}

void LifetimeLogProcess::saveATimePoint(std::ofstream& aFile,
                                         const double aTime,
                                         const unsigned anIteration,
                                         const unsigned aCnt)
{
  logFile(aTime, 0, 0, 0, theLogValues[aCnt][0]/anIteration);
}

void LifetimeLogProcess::saveFile()
{
  unsigned aCompletedIterations(1);
  if(!Iterations)
    {
      cout << "Saving data in: " << FileName.c_str() << std::endl;
      theLogFile.open(FileName.c_str(), std::ios::trunc);
      aCompletedIterations = theTotalIterations;
    }
  else
    {
      return;
    }
  saveFileHeader(theLogFile);
  saveFileData(theLogFile, aCompletedIterations);
}

void LifetimeLogProcess::doPreLog()
{
  cout << "Iterations left:" << Iterations << " of " << theTotalIterations
    << std::endl;
  if(theTotalIterations == 1)
    {
      cout << "Saving single iteration data in: " << FileName.c_str()
        << std::endl;
      theLogFile.open(FileName.c_str(), std::ios::trunc);
      saveFileHeader(theLogFile);
    }
}

void LifetimeLogProcess::interruptedPre(ReactionProcess* aProcess)
{
  if(isDimerizationReaction[aProcess->getID()])
    {
      saveDimerizingMonomerTag(aProcess);
    }
  else if(isRemoveDimerReaction[aProcess->getID()])
    {
      logTrackedDimer(aProcess->getA(), aProcess->getMoleculeA());
    }
  else if(aProcess->getA() && isTrackedSpecies[aProcess->getA()->getID()])
    {
      logTrackedMolecule(aProcess, aProcess->getA(), aProcess->getMoleculeA());
    }
  else if(aProcess->getB() && isTrackedSpecies[aProcess->getB()->getID()])
    {
      logTrackedMolecule(aProcess, aProcess->getB(), aProcess->getMoleculeB());
    }
}

void LifetimeLogProcess::interruptedPost(ReactionProcess* aProcess)
{
  if(isDedimerizationReaction[aProcess->getID()])
    {
      initDedimerizingMonomerTag(aProcess);
    }
  else if(isAddDimerReaction[aProcess->getID()])
    {
      initTrackedDimer(aProcess->getC(), aProcess->getC()->size()-1);
    }
  else if(aProcess->getC() && isTrackedSpecies[aProcess->getC()->getID()])
    {
      initTrackedMolecule(aProcess->getC(), aProcess->getC()->size()-1);
    }
  else if(aProcess->getD() && isTrackedSpecies[aProcess->getD()->getID()])
    {
      initTrackedMolecule(aProcess->getD(), aProcess->getD()->size()-1);
    }
}

void LifetimeLogProcess::logTrackedMolecule(ReactionProcess* aProcess,
                                            Species* aSpecies,
                                            const Voxel* aMolecule)
{
  if(isBindingSiteReaction[aProcess->getID()])
    {
      if(aProcess->getMoleculeC() &&
         isTrackedSpecies[aProcess->getC()->getID()])
        {
          return;
        }
      else if(aProcess->getMoleculeD() &&
              isTrackedSpecies[aProcess->getD()->getID()])
        {
          return;
        }
    }
  const unsigned anIndex(aSpecies->getIndex(aMolecule));
  logTag(aSpecies, aSpecies->getTag(anIndex), anIndex);
}

void LifetimeLogProcess::logTrackedDimer(Species* aSpecies,
                                         const Voxel* aMolecule)
{
  const unsigned anIndex(aSpecies->getIndex(aMolecule));
  Tag& aTag(aSpecies->getTag(anIndex));
  logTag(aSpecies, aTag, anIndex);
  logTag(aSpecies, theDimerizingMonomerTags[aTag.molID], anIndex);
}

double LifetimeLogProcess::getAverageDiffusion()
{
  double activeSquaredDisplacement(0);
  double activeTotalTime(0);
  const double now(getStepper()->getCurrentTime());
  double aCoeff(0);
  for(unsigned i(0); i != theSpecies.size(); ++i)
    {
      if(isTrackedSpecies[i])
        {
          Species* aSpecies(theSpecies[i]);
          aCoeff = aSpecies->getDimension()*2;
          /*
          double disp(0);
          double time(0);
          */
          for(unsigned j(0); j != aSpecies->size(); ++j)
            {
              /*
              disp += aSpecies->getSquaredDisplacement(j);
              time += now-theTagTimes[aSpecies->getTag(j).molID];
              */
              activeSquaredDisplacement += aSpecies->getSquaredDisplacement(j);
              activeTotalTime += now-theTagTimes[aSpecies->getTag(j).molID];
            } 
          /*
          std::cout << now << " species:" << aSpecies->getIDString() << 
          " dif:" << disp/(time*aCoeff) << " number:" << aSpecies->size() <<
          std::endl;
          */
        }
    }
  return (completedSquaredDisplacement+activeSquaredDisplacement)/
    ((totalDuration+activeTotalTime)*aCoeff);
}

void LifetimeLogProcess::logTag(Species* aSpecies, Tag& aTag,
                                const unsigned anIndex)
{
  availableTagIDs.push_back(aTag.molID);
  double startTime(theTagTimes[aTag.molID]);
  double endTime(getStepper()->getCurrentTime());
  double duration(endTime-startTime);
  double averageDiffusion(getAverageDiffusion());
  double squaredDisplacement(aSpecies->getSquaredDisplacement(anIndex));
  totalDuration += duration;
  completedSquaredDisplacement += squaredDisplacement;
  ++logCnt;
  logFile(endTime, duration, squaredDisplacement, startTime, averageDiffusion);
  if(Verbose)
    {
      std::cout << "mean kon:" << konCnt/endTime << 
        " mean koff:" << 1/(totalDuration/logCnt);
      for(unsigned i(0); i != theSpecies.size(); ++i)
        {
          if(isTrackedSpecies[i])
            {
              trackedSpeciesSize[i] += theSpecies[i]->size();
              std::cout << " " << theSpecies[i]->getIDString() <<
                "[" << trackedSpeciesSize[i]/logCnt << "]";
              totalSize += theSpecies[i]->size();
            }
        }
      std::cout << " all[" << totalSize/logCnt << "]" << std::endl;
    }
}

void LifetimeLogProcess::logFile(const double time, const double duration,
                                 const double squaredDisplacement,
                                 const double startTime,
                                 const double averageDiffusion)
{
  double averageKoff(libecs::INF);
  if(totalDuration > 0)
    {
      averageKoff = logCnt/totalDuration;
    }
  theLogFile << time << "," << duration << "," <<
    squaredDisplacement << "," << startTime << "," << 
    konCnt << "," << konCnt/time << "," << averageKoff << "," <<
    averageDiffusion << std::endl;
}

void LifetimeLogProcess::saveFileHeader(std::ofstream& aFile)
{
  aFile << "Time, Duration of";
  for(unsigned i(0); i != isTrackedSpecies.size(); ++i)
    {
      if(isTrackedSpecies[i])
        {
          aFile << ":" << theSpecies[i]->getIDString();
        }
    }
  aFile << "(s), Squared Displacement(m^2), Start time(s), " <<
    "kon counts [#kon] (some still haven't koff), Average kon so far " <<
    "[#kon/time](1/s), Average koff so far(1/s), Average diffusion coefficient"
    << " so far(m^2/s)" <<
    std::endl; 
}

void LifetimeLogProcess::initDedimerizingMonomerTag(ReactionProcess* aProcess)
{
  //Monomer C's tag will be the dimers tag, we need to set
  //monomer D's tag to the previously saved tag:
  Species* C(aProcess->getC());
  Species* D(aProcess->getD());
  //Since it is first order dedimerization reaction, the reaction will be
  //executed by SNRP and STPL both of which have valid moleculeC and moleculeD.
  //We only need to access moleculeC and moleculeD as C->size()-1 index
  //when the reaction is DIRP since moleculeC and moleculeD are not valid
  //for DIRP.
  const Voxel* molC(aProcess->getMoleculeC());
  const Voxel* molD(aProcess->getMoleculeD());
  const unsigned indexC(C->getIndex(molC));
  const unsigned indexD(D->getIndex(molD));
  Tag& aTag(D->getTag(indexD));
  aTag = theDimerizingMonomerTags[C->getTag(indexC).molID];
}

void LifetimeLogProcess::initTrackedDimer(Species* aSpecies,
                                          const unsigned anIndex)
{
  initTrackedMolecule(aSpecies, anIndex);
  Tag& tagA(aSpecies->getTag(anIndex));
  Tag tagB(tagA);
  addTagTime(tagB);
  theDimerizingMonomerTags.resize(theTagTimes.size());
  theDimerizingMonomerTags[tagA.molID] = tagB;
  ++konCnt;
}

void LifetimeLogProcess::initTrackedMolecule(Species* aSpecies,
                                             const unsigned anIndex)
{
  aSpecies->resetTagOrigin(anIndex);
  Tag& aTag(aSpecies->getTag(anIndex));
  addTagTime(aTag);
  ++konCnt;
}


void LifetimeLogProcess::addTagTime(Tag& aTag)
{
  if(availableTagIDs.size())
    {
      aTag.molID = availableTagIDs.back();
      theTagTimes[availableTagIDs.back()] = getStepper()->getCurrentTime();
      availableTagIDs.pop_back();
    }
  else
    {
      aTag.molID = theTagTimes.size();
      theTagTimes.push_back(getStepper()->getCurrentTime());
    }
}

void LifetimeLogProcess::saveDimerizingMonomerTag(ReactionProcess* aProcess)
{
  //Monomer A's tag will be taken up by the product dimer, so save the 
  //monomer B tag before it is destroyed after the dimerization
  //reaction:
  Species* A(aProcess->getA());
  Species* B(aProcess->getB());
  const Voxel* molA(aProcess->getMoleculeA());
  const Voxel* molB(aProcess->getMoleculeB());
  const unsigned indexA(A->getIndex(molA));
  const unsigned indexB(B->getIndex(molB));
  theDimerizingMonomerTags.resize(theTagTimes.size());
  theDimerizingMonomerTags[A->getTag(indexA).molID] = B->getTag(indexB);
}


bool LifetimeLogProcess::isDependentOnPre(const Process* aProcess)
{
  const ReactionProcess* aReactionProcess(dynamic_cast<
                                          const ReactionProcess*>(aProcess));
  if(!aReactionProcess)
    {
      return false;
    }
  if(isDimerizationReaction[aReactionProcess->getID()])
    {
      return true;
    }
  Species* A(aReactionProcess->getA());
  Species* B(aReactionProcess->getB());
  Species* C(aReactionProcess->getC());
  Species* D(aReactionProcess->getD());
  // Dimer -> 0
  if(A && isTrackedDimerSpecies[A->getID()] && 
     getVariableNetCoefficient(aProcess, A->getVariable()) &&
     (!C || (C && !isTrackedSpecies[C->getID()])))
    {
      isRemoveDimerReaction[aReactionProcess->getID()] = true;
      return true;
    }
  if((A && isTrackedSpecies[A->getID()]) && 
     getVariableNetCoefficient(aProcess, A->getVariable()) ||
     (B && isTrackedSpecies[B->getID()] && 
      getVariableNetCoefficient(aProcess, B->getVariable())))
    {
      if((C && !isTrackedSpecies[C->getID()]) ||
         (D && !isTrackedSpecies[D->getID()]))
        {
          if((C && isTrackedSpecies[C->getID()]) ||
             (D && isTrackedSpecies[D->getID()]))
            {
              isBindingSiteReaction[aReactionProcess->getID()] = true;
            }
          return true;
        }
    }
  return false;
}

bool LifetimeLogProcess::isDependentOnPost(const Process* aProcess)
{
  const ReactionProcess* aReactionProcess(dynamic_cast<
                                          const ReactionProcess*>(aProcess));
  if(!aReactionProcess)
    {
      return false;
    }
  if(isDedimerizationReaction[aReactionProcess->getID()])
    {
      return true;
    }
  Species* A(aReactionProcess->getA());
  Species* B(aReactionProcess->getB());
  Species* C(aReactionProcess->getC());
  Species* D(aReactionProcess->getD());
  // 0 -> Dimer
  if(!A && C && isTrackedDimerSpecies[C->getID()])
    {
      isAddDimerReaction[aReactionProcess->getID()] = true;
      return true;
    }
  if((C && isTrackedSpecies[C->getID()]) ||
     (D && isTrackedSpecies[D->getID()]))
    {
      if((A && isTrackedSpecies[A->getID()]) ||
         (B && isTrackedSpecies[B->getID()]))
        {
          return false;
        }
      return true;
    }
  return false;
}

}
