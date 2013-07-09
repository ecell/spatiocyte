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
  isTrackedDimerSpecies.resize(theSpecies.size(), false);
  isUntrackedSpecies.resize(theSpecies.size(), false);
  unsigned cntTracked(0);
  unsigned cntUntracked(0);
  for(VariableReferenceVector::iterator
      i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      Species* aSpecies(theSpatiocyteStepper->variable2species(
                                     (*i).getVariable())); 
      if((*i).getCoefficient() == -1)
        {
          isTrackedSpecies[aSpecies->getID()] = true;
          aSpecies->setIsTagged();
          ++cntTracked;
        }
      else if((*i).getCoefficient() == 1)
        {
          isUntrackedSpecies[aSpecies->getID()] = true;
          ++cntUntracked;
        }
    }
  setTrackedDimerSpecies();
  //This process is only for tracking the lifetime of one particular species
  //that goes through several intermediate species before finally untracked
  //when it becomes a final untracked species. For example if we want to track
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
                      "nonHD variable reference with -1 " +
                      "coefficient as the tracked species, " +
                      "but none is given."); 
    }
  if(!cntUntracked)
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + 
                      "]: Lifetime logging requires at least one " +
                      "nonHD variable reference with 1 " +
                      "coefficient as the untracked species, " +
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
  theLogFile.open(FileName.c_str(), std::ios::trunc);
}

void LifetimeLogProcess::initializeFifth()
{
  theInterval = LogEnd;
  if(!LogStart)
    {
      LogStart = theInterval;
    }
  theTime = std::max(LogStart-theInterval, getStepper()->getMinStepInterval());
  thePriorityQueue->move(theQueueID);
}

void LifetimeLogProcess::fire()
{
  if(theTime < LogStart)
    {
      cout << "Iterations left:" << Iterations << " of " <<
        theTotalIterations << std::endl;    
      //Since there is only one file even when there are more than one
      //iterations, save header only once:
      if(Iterations == theTotalIterations)
        {
          saveFileHeader(theLogFile);
        }
    }
  if(theTime >= LogEnd && Iterations > 0)
    {
      --Iterations;
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

void LifetimeLogProcess::interruptedPre(ReactionProcess* aProcess)
{
  //std::cout << "me:" << getIDString() << " interruptedPre:" << aProcess->getIDString() << std::endl;
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
  //std::cout << "me:" << getIDString() << " interruptedPost:" << aProcess->getIDString() << std::endl;
  if(isDedimerizationReaction[aProcess->getID()])
    {
      initDedimerizingMonomerTag(aProcess);
    }
  else if(isAddDimerReaction[aProcess->getID()])
    {
      initDimerNewTags(aProcess->getC());
    }
  else if(aProcess->getC() && isTrackedSpecies[aProcess->getC()->getID()])
    {
      initTrackedMolecule(aProcess->getC());
    }
  else if(aProcess->getD() && isTrackedSpecies[aProcess->getD()->getID()])
    {
      initTrackedMolecule(aProcess->getD());
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
  logTag(aSpecies, theDimerizingMonomerTags[aTag.id], anIndex);
}

void LifetimeLogProcess::logTag(Species* aSpecies, Tag& aTag,
                                const unsigned anIndex)
{
  const Point aPoint(aSpecies->getPoint(anIndex));
  const Point anOrigin(aSpecies->coord2point(aTag.origin));
  availableTagIDs.push_back(aTag.id);
  double aTime(getStepper()->getCurrentTime());
  double duration(aTime-theTagTimes[aTag.id]);
  totalDuration += duration;
  ++logCnt;
  theLogFile << std::setprecision(15) << duration << "," <<
    distance(aPoint, anOrigin)*2*theSpatiocyteStepper->getVoxelRadius() << ","
    << theTagTimes[aTag.id] << "," << aTime << "," << konCnt << "," <<
    1/(totalDuration/logCnt) << std::endl;
  std::cout << "average koff:" << 1/(totalDuration/logCnt) << std::endl;
}

void LifetimeLogProcess::saveFileHeader(std::ofstream& aFile)
{
  aFile << "Duration of";
  for(unsigned i(0); i != isTrackedSpecies.size(); ++i)
    {
      if(isTrackedSpecies[i])
        {
          aFile << ":" << theSpecies[i]->getIDString();
        }
    }
  aFile << "(s), Distance(m), Start time(s), End time(s), " <<
    "kon counts (some still haven't koff), Average koff so far(1/s)" <<
    std::endl; 
}

void LifetimeLogProcess::initDedimerizingMonomerTag(ReactionProcess* aProcess)
{
  //Monomer C's tag will be the dimers tag, we need to set
  //monomer D's tag to the previously saved tag:
  //std::cout << "1" << std::endl;
  Species* C(aProcess->getC());
  //std::cout << "2" << std::endl;
  Species* D(aProcess->getD());
  //std::cout << "3" << std::endl;
  //Since it is first order dedimerization reaction, the reaction will be
  //executed by SNRP and STPL both of which have valid moleculeC and moleculeD.
  //We only need to access moleculeC and moleculeD as C->size()-1 index
  //when the reaction is DIRP since moleculeC and moleculeD are not valid
  //for DIRP.
  const Voxel* molC(aProcess->getMoleculeC());
  //std::cout << "4" << std::endl;
  const Voxel* molD(aProcess->getMoleculeD());
  //std::cout << "5" << std::endl;
  const unsigned indexC(C->getIndex(molC));
  //std::cout << "6" << std::endl;
  const unsigned indexD(D->getIndex(molD));
  //std::cout << "7" << std::endl;
  Tag& aTag(D->getTag(indexD));
  //std::cout << "8 id:" << C->getTag(indexC).id << " size:" << theDimerizingMonomerTags.size() << std::endl;
  aTag = theDimerizingMonomerTags[C->getTag(indexC).id];
  //std::cout << "9" << std::endl;
}

void LifetimeLogProcess::initDimerNewTags(Species* aDimer)
{
  const unsigned anIndex(aDimer->size()-1);
  Tag& tagA(aDimer->getTag(anIndex));
  tagA.origin = aDimer->getCoord(anIndex);
  addTagTime(tagA);
  Tag tagB(tagA);
  addTagTime(tagB);
  theDimerizingMonomerTags.resize(theTagTimes.size());
  theDimerizingMonomerTags[tagA.id] = tagB;
  konCnt += 2;
}

void LifetimeLogProcess::addTagTime(Tag& aTag)
{
  if(availableTagIDs.size())
    {
      aTag.id = availableTagIDs.back();
      theTagTimes[availableTagIDs.back()] = getStepper()->getCurrentTime();
      availableTagIDs.pop_back();
    }
  else
    {
      aTag.id = theTagTimes.size();
      theTagTimes.push_back(getStepper()->getCurrentTime());
    }
}

void LifetimeLogProcess::saveDimerizingMonomerTag(ReactionProcess* aProcess)
{
  //std::cout << "saving" << std::endl;
  //Monomer A's tag will be taken up by the product dimer, so save the 
  //monomer B tag before it is destroyed after the dimerization
  //reaction:
  Species* A(aProcess->getA());
  Species* B(aProcess->getB());
  const Voxel* molA(aProcess->getMoleculeA());
  const Voxel* molB(aProcess->getMoleculeB());
  const unsigned indexA(A->getIndex(molA));
  const unsigned indexB(B->getIndex(molB));
  //std::cout << "tag time size:" << theTagTimes.size() << std::endl;
  theDimerizingMonomerTags.resize(theTagTimes.size());
  theDimerizingMonomerTags[A->getTag(indexA).id] = B->getTag(indexB);
}


void LifetimeLogProcess::initTrackedMolecule(Species* aSpecies)
{
  const unsigned anIndex(aSpecies->size()-1);
  Tag& aTag(aSpecies->getTag(anIndex));
  aTag.origin = aSpecies->getCoord(anIndex);
  addTagTime(aTag);
  ++konCnt;
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
     (!C || (C && !isTrackedSpecies[C->getID()])))
    {
      isRemoveDimerReaction[aReactionProcess->getID()] = true;
      return true;
    }
  if((A && isTrackedSpecies[A->getID()]) ||
     (B && isTrackedSpecies[B->getID()]))
    {
      if((C && isUntrackedSpecies[C->getID()]) ||
         (D && isUntrackedSpecies[D->getID()]))
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
