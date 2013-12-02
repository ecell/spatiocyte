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


#ifndef __LifetimeLogProcess_hpp
#define __LifetimeLogProcess_hpp

#include <libecs/IteratingLogProcess.hpp>
#include <libecs/ReactionProcess.hpp>

namespace libecs
{

LIBECS_DM_CLASS(LifetimeLogProcess, IteratingLogProcess)
{ 
public:
  LIBECS_DM_OBJECT(LifetimeLogProcess, Process)
    {
      INHERIT_PROPERTIES(IteratingLogProcess);
    }
  LifetimeLogProcess():
    konCnt(0),
    logCnt(0),
    totalSize(0),
    totalDuration(0)
  {
    FileName = "LifetimeLog.csv";
  }
  virtual ~LifetimeLogProcess() {}
  virtual void initialize();
  virtual void initializeFirst();
  virtual void initializeSecond();
  virtual void initializeFourth();
  virtual void initializeFifth();
  virtual void initializeLastOnce();
  virtual void fire();
  virtual void interruptedPre(ReactionProcess*);
  virtual void interruptedPost(ReactionProcess*);
  virtual bool isDependentOnPre(const Process*);
  virtual bool isDependentOnPost(const Process*);
  virtual void saveFileHeader(std::ofstream&);
private:
  void logTrackedMolecule(ReactionProcess*, Species*, const Voxel*);
  void logTrackedDimer(Species*, const Voxel*);
  void initTrackedMolecule(Species*, const unsigned);
  void initTrackedDimer(Species*, const unsigned);
  void saveDimerizingMonomerTag(ReactionProcess*);
  void initDedimerizingMonomerTag(ReactionProcess*);
  void setTrackedDimerSpecies();
  void logTag(Species*, Tag&, const unsigned);
  void addTagTime(Tag&);
private:
  unsigned konCnt;
  unsigned logCnt;
  double totalSize;
  double totalDuration;
  std::vector<bool> isAddDimerReaction;
  std::vector<bool> isBindingSiteReaction;
  std::vector<bool> isDedimerizationReaction;
  std::vector<bool> isDimerizationReaction;
  std::vector<bool> isRemoveDimerReaction;
  std::vector<bool> isTrackedSpecies;
  std::vector<bool> isTrackedDimerSpecies;
  std::vector<unsigned> availableTagIDs;
  std::vector<Tag> theDimerizingMonomerTags;
  std::vector<double> theTagTimes;
  std::vector<double> trackedSpeciesSize;
};

}

#endif /* __LifetimeLogProcess_hpp */
