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

#include <libecs/MicrotubuleProcess.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(MicrotubuleProcess, Process); 

void MicrotubuleProcess::initialize() {
  FilamentProcess::initialize();
  nMonomerPitch = MonomerPitch/(VoxelRadius*2);
}

void MicrotubuleProcess::setSubunitStart() {
  Point R; //Initialize a random point on the plane attached at the minus end
  if(Minus.x != Plus.x)
    {
      R.y = 10;
      R.z = 30; 
      R.x = (Minus.x*lengthVector.x+Minus.y*lengthVector.y-R.y*lengthVector.y+
             Minus.z*lengthVector.z-R.z*lengthVector.z)/lengthVector.x;
    }
  else if(Minus.y != Plus.y)
    {
      R.x = 10; 
      R.z = 30;
      R.y = (Minus.x*lengthVector.x-R.x*lengthVector.x+Minus.y*lengthVector.y+
             Minus.z*lengthVector.z-R.z*lengthVector.z)/lengthVector.y;
    }
  else
    {
      R.x = 10; 
      R.y = 30;
      R.z = (Minus.x*lengthVector.x-R.x*lengthVector.x+Minus.y*lengthVector.y-
             R.y*lengthVector.y+Minus.z*lengthVector.z)/lengthVector.z;
    }
  //The direction vector from the minus end to the random point, R
  Point D(sub(R, Minus));
  norm_(D);
  subunitStart = disp(Minus, D, nRadius);
}

void MicrotubuleProcess::initializeThird() {
  if(!isCompartmentalized)
    {
      thePoints.resize(endCoord-subStartCoord);
      vacStartIndex = theVacantSpecies->size();
      intStartIndex = theInterfaceSpecies->size();
      initializeVectors();
      initializeFilaments(subunitStart, Filaments, Subunits, nMonomerPitch,
                          theVacantSpecies, subStartCoord);
      elongateFilaments(theVacantSpecies, subStartCoord, Filaments, Subunits,
                        nDiffuseRadius);
      connectFilaments(subStartCoord, Filaments, Subunits);
      setDiffuseSize(subStartCoord, lipStartCoord);
      connectTrailSubunits(subStartCoord, Filaments, Subunits);
      setTrailSize(subStartCoord, lipStartCoord);
      setGrid(theVacantSpecies, theVacGrid, subStartCoord);
      interfaceSubunits();
      isCompartmentalized = true;
    }
  theVacantSpecies->setIsPopulated();
  theInterfaceSpecies->setIsPopulated();
  theMinusSpecies->setIsPopulated();
  thePlusSpecies->setIsPopulated();
}

void MicrotubuleProcess::initializeFilaments(Point& aStartPoint, unsigned aRows,
                                             unsigned aCols, double aRadius,
                                             Species* aVacant,
                                             unsigned aStartCoord) {
  Voxel* aVoxel(addCompVoxel(0, 0, aStartPoint, aVacant, aStartCoord, aCols));
  theMinusSpecies->addMolecule(aVoxel);
  Point U(aStartPoint);
  for(unsigned i(1); i < aRows; ++i)
    {
      double angle(2*M_PI/aRows);
      rotatePointAlongVector(U, Minus, lengthVector, angle);
      disp_(U, lengthVector, aRadius/(aRows-1));
      aVoxel = addCompVoxel(i, 0, U, aVacant, aStartCoord, aCols);
      theMinusSpecies->addMolecule(aVoxel);
    }
}

}
