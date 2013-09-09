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


#ifndef __SurfaceCompartmentProcess_hpp
#define __SurfaceCompartmentProcess_hpp

#include <sstream>
#include <libecs/CompartmentProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>

LIBECS_DM_CLASS(SurfaceCompartmentProcess, CompartmentProcess)
{ 
public:
  LIBECS_DM_OBJECT(SurfaceCompartmentProcess, Process)
    {
      PROPERTYSLOT_SET_GET(String, FileName);
    }
  SurfaceCompartmentProcess() {}
  virtual ~SurfaceCompartmentProcess() {}
  SIMPLE_SET_GET_METHOD(String, FileName);
  virtual void prepreinitialize()
    {
      SpatiocyteProcess::prepreinitialize();
      theInterfaceVariable = createVariable("Interface");
      theVacantVariable = createVariable("Vacant");
    }
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
      theVacantSpecies = theSpatiocyteStepper->addSpecies(theVacantVariable);
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          Species* aSpecies(theSpatiocyteStepper->variable2species(
                                   (*i).getVariable())); 
          if(!(*i).getCoefficient())
            {
              theVacantCompSpecies.push_back(aSpecies);
            }
        }
      VoxelRadius = theSpatiocyteStepper->getVoxelRadius();
      if(!SubunitRadius)
        {
          SubunitRadius = VoxelRadius;
        }
      //SubunitRadius is the actual radius of a molecule.
      nSubunitRadius = SubunitRadius/(VoxelRadius*2);
    }
  virtual void initializeFirst()
    {
      SpatiocyteProcess::initializeFirst();
      theComp = new Comp;
      theVacantSpecies->setIsCompVacant();
      theVacantSpecies->setIsOffLattice();
      theVacantSpecies->setComp(theComp);
      for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
        {
          theVacantCompSpecies[i]->setIsOffLattice();
          //setVacantSpecies must be declared here since it needs
          //to be overwritten by DiffusionProcess in initializeSecond:
          theVacantCompSpecies[i]->setVacantSpecies(theVacantSpecies);
          theVacantCompSpecies[i]->setComp(theComp);
        }
    }
  virtual unsigned getLatticeResizeCoord(unsigned);
  virtual void initializeThird();
  virtual void printParameters() {}
protected:
  String FileName;
};

#endif /* __SurfaceCompartmentProcess_hpp */




