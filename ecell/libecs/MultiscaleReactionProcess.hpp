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


#ifndef __MultiscaleReactionProcess_hpp
#define __MultiscaleReactionProcess_hpp

#include <sstream>
#include <libecs/DiffusionInfluencedReactionProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>

namespace libecs
{

LIBECS_DM_CLASS(MultiscaleReactionProcess, DiffusionInfluencedReactionProcess)
{ 
  typedef bool (MultiscaleReactionProcess::*Method)(Voxel*, Voxel*,
                                              const unsigned, const unsigned);
public:
  LIBECS_DM_OBJECT(MultiscaleReactionProcess, Process)
    {
      INHERIT_PROPERTIES(DiffusionInfluencedReactionProcess);
    }
  MultiscaleReactionProcess() {}
  virtual ~MultiscaleReactionProcess() {}
  virtual void initializeThird();
  virtual void initializeMultiscaleWalkBindUnbind();
  virtual void initializeMultiscaleCompReaction();
  virtual void setReactMethod();
  virtual void bind(Voxel* aVoxel, const unsigned multiIdx)
    {
      const unsigned index(aVoxel->idx%theStride);
      M->addMoleculeInMulti(aVoxel, multiIdx, N->getTag(index).boundCnt);
      N->removeMoleculeBoundDirect(index);
    }
  virtual void unbind(Voxel* aVoxel)
    {
      const unsigned index(aVoxel->idx%theStride);
      N->addMoleculeExMulti(aVoxel, M->getTag(index).boundCnt);
      M->removeMoleculeBoundDirect(index);
    }
  virtual bool react(Voxel* molA, Voxel* molB, const unsigned indexA,
                     const unsigned indexB)
    {
      //We need an explicit declaration of the react method here
      //although it is identical to DIRP's react method because
      //this-> in DIRP would refer to its own method.
      moleculeA = molA;
      moleculeB = molB;
      if((this->*reactM)(molA, molB, indexA, indexB))
        {
          interruptProcessesPost();
          return true;
        }
      return false;
    }
protected:
  unsigned getIdx(Species*, Voxel*, const unsigned);
  bool reactMuAtoMuC(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactAllMuAtoMuC(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuBtoMuC(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactAtoC_MuBtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuAtoMuC_BtoD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactBtoC_MuAtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuBtoMuC_AtoD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactAeqC_MuBtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuAeqMuC_BtoD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuAeqMuC_MuBtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactBeqC_MuAtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuBeqMuC_AtoD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuBtoMuC_AeqD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactBtoC_MuAeqMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuAtoMuC_BeqD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuAtoMuC_MuBeqMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactAtoC_MuBeqMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactAtoC_Multi(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactBtoC_Multi(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuAtoMuC_MuBtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  bool reactMuBeqMuC_MuAtoMuD(Voxel*, Voxel*, const unsigned, const unsigned);
  void setReactVarC_D();
  void setReactVarD_C();
  void setReactD();
  void setReactC();
protected:
  Species* M;
  Species* N;
  Species* theMultiscale;
  Method reactM;
};

}

#endif /* __MultiscaleReactionProcess_hpp */
