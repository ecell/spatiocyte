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


#ifndef __FilamentProcess_hpp
#define __FilamentProcess_hpp

#include <sstream>
#include <libecs/MicrotubuleProcess.hpp>

namespace libecs
{

LIBECS_DM_CLASS(FilamentProcess, MicrotubuleProcess)
{ 
public:
  LIBECS_DM_OBJECT(FilamentProcess, Process)
    {
      INHERIT_PROPERTIES(MicrotubuleProcess);
    }
  FilamentProcess()
  {
    Filaments = 1;
  }
  virtual ~FilamentProcess() {}
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      SpatiocyteProcess::initialize();
      theInterfaceSpecies = theSpatiocyteStepper->addSpecies(
                                                       theInterfaceVariable);
      theInterfaceSpecies->setIsInterface();
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          Species* aSpecies(theSpatiocyteStepper->variable2species(
                                   (*i).getVariable())); 
          if((*i).getCoefficient())
            {
              if((*i).getCoefficient() == -1)
                {
                  if(theVacantSpecies)
                    {
                      THROW_EXCEPTION(ValueError, String(
                                      getPropertyInterface().getClassName()) +
                                      "[" + getFullID().asString() + 
                                      "]: A FilamentProcess requires only " +
                                      "one vacant variable reference with -1 " +
                                      "coefficient as the vacant species of " +
                                      "the Filament compartment, but " +
                                      getIDString(theVacantSpecies) + " and " +
                                      getIDString(aSpecies) + " are given."); 
                    }
                  theVacantSpecies = aSpecies;
                }
            }
          else
            {
              theFilamentSpecies.push_back(aSpecies);
            }
        }
      if(!theFilamentSpecies.size())
        {
          THROW_EXCEPTION(ValueError, String(
                          getPropertyInterface().getClassName()) +
                          "[" + getFullID().asString() + 
                          "]: A FilamentProcess requires at least one " +
                          "nonHD variable reference with zero coefficient " +
                          "as the filament species, but none is given."); 
        }
      if(!theVacantSpecies)
        {
          THROW_EXCEPTION(ValueError, String(
                          getPropertyInterface().getClassName()) +
                          "[" + getFullID().asString() + 
                          "]: A FilamentProcess requires one " +
                          "nonHD variable reference with negative " +
                          "coefficient as the vacant species, " +
                          "but none is given."); 
        }
      if(!DiffuseRadius)
        {
          if(SubunitRadius)
            {
              DiffuseRadius = SubunitRadius;
            }
          else
            {
              DiffuseRadius = theSpatiocyteStepper->getVoxelRadius();
            }
        }
      if(!SubunitRadius)
        {
          SubunitRadius = DiffuseRadius;
        }
      VoxelRadius = theSpatiocyteStepper->getVoxelRadius();
      //Normalized off-lattice voxel radius:
      nSubunitRadius = SubunitRadius/(VoxelRadius*2);
      nDiffuseRadius = DiffuseRadius/(VoxelRadius*2);
      Radius = SubunitRadius;
      nRadius = Radius/(VoxelRadius*2);
      nGridSize = 10*nDiffuseRadius;
    }
  virtual void initializeFirst()
    {
      CompartmentProcess::initializeFirst();
      for(unsigned i(0); i != theFilamentSpecies.size(); ++i)
        {
          theFilamentSpecies[i]->setIsOffLattice();
          theFilamentSpecies[i]->setDimension(1);
          theFilamentSpecies[i]->setVacantSpecies(theVacantSpecies);
          theFilamentSpecies[i]->setComp(theComp);
          theFilamentSpecies[i]->resetFixedAdjoins();
        }
    }
  virtual unsigned getLatticeResizeCoord(unsigned);
  virtual void setCompartmentDimension();
  virtual void initializeVectors();
  virtual void initializeFilaments(Point&, unsigned, unsigned, double, Species*,
                                   unsigned);
  virtual void initializeThird();
  virtual void setSubunitStart();
  virtual void connectFilaments(unsigned, unsigned, unsigned);
  virtual void elongateFilaments(Species*, unsigned, unsigned, unsigned,
                                 double);
  virtual void addPlaneIntersectInterfaceVoxel(Voxel&, Point&);
  virtual bool isInside(Point&);
  virtual bool isOnAboveSurface(Point&);
  void connectTrailTubulins(unsigned, unsigned, unsigned);
  void setTrailSize(unsigned, unsigned);
protected:
  std::vector<Species*> theFilamentSpecies;
};

}

#endif /* __FilamentProcess_hpp */




