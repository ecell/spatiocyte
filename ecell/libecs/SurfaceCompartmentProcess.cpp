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

#include <libecs/SurfaceCompartmentProcess.hpp>
#include <libecs/SpatiocyteVector.hpp>

LIBECS_DM_INIT_STATIC(SurfaceCompartmentProcess, Process); 

unsigned SurfaceCompartmentProcess::getLatticeResizeCoord(unsigned aStartCoord)
{
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  if(aComp->diffusiveComp)
    {
      Comp* aDiffusiveComp(aComp->diffusiveComp);
      if(aDiffusiveComp->baseInterfaceID == theSpecies.size())
        {
          aDiffusiveComp->baseInterfaceID = theInterfaceSpecies->getID();
        }
      else
        {
          theInterfaceSpecies = theSpecies[aDiffusiveComp->baseInterfaceID];
        }
    }
  aComp->baseInterfaceID = theInterfaceSpecies->getID();
  *theComp = *aComp;
  theVacantSpecies->resetFixedAdjoins();
  theVacantSpecies->setMoleculeRadius(SubunitRadius);
  Origin = aComp->centerPoint;
  theDimension = 2;
  theComp->dimension = theDimension;
  subStartCoord = aStartCoord;
  endCoord = aStartCoord+1000;
  return endCoord-aStartCoord;
}

void SurfaceCompartmentProcess::initializeThird()
{
}

