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

void MicrotubuleProcess::initializeVectors() {
  FilamentProcess::initializeVectors();
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

void MicrotubuleProcess::initializeCompartment() {
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
  /*
  theSpecies[0]->setIsPopulated();
  theSpecies[9]->setIsPopulated();
  theSpecies[10]->setIsPopulated();
  */
}

void MicrotubuleProcess::initializeFilaments(Point& aStartPoint, unsigned aRows,
                                             unsigned aCols, double aRadius,
                                             Species* aVacant,
                                             unsigned aStartCoord) {
  Voxel* aVoxel(addCompVoxel(0, 0, aStartPoint, aVacant, aStartCoord, aCols));
  theMinusSpecies->addMolecule(aVoxel);
  theMinusPoints.push_back(aStartPoint);
  Point U(aStartPoint);
  for(unsigned i(1); i < aRows; ++i)
    {
      double angle(2*M_PI/aRows);
      rotatePointAlongVector(U, Minus, lengthVector, angle);
      disp_(U, lengthVector, aRadius/(aRows-1));
      aVoxel = addCompVoxel(i, 0, U, aVacant, aStartCoord, aCols);
      theMinusSpecies->addMolecule(aVoxel);
      theMinusPoints.push_back(U);
    }
}

unsigned MicrotubuleProcess::getNearestFilament(Point& aPoint)
{
  double nearestDist(libecs::INF);
  unsigned nearestFilament(0);
  for(unsigned i(0); i != theMinusPoints.size(); ++i)
    {
      double aDist(point2lineDist(aPoint, lengthVector, theMinusPoints[i]));
      if(aDist < nearestDist)
        {
          nearestDist = aDist;
          nearestFilament = i;
        }
    }
  return nearestFilament;
}


//Is inside the parent compartment and confined by the length of the MT:
bool MicrotubuleProcess::isInside(Point& aPoint)
{
  if(nVoxelRadius >= nRadius)
    {
      return FilamentProcess::isInside(aPoint);
    }
  unsigned nearestFilament(getNearestFilament(aPoint));
  double aLengthDisplace(dot(lengthVector, theMinusPoints[nearestFilament]));
  double dispA(point2planeDisp(aPoint, lengthVector, aLengthDisplace));
  if(dispA >= -nDiffuseRadius/2 && dispA <= nLength+nDiffuseRadius/2)
    {
      return true;
    }
  return false;
}


//Overload FilamentProcess::extendInterfacesOverSurface():
void MicrotubuleProcess::extendInterfacesOverSurface()
{
  if(nVoxelRadius >= nRadius)
    {
      FilamentProcess::extendInterfacesOverSurface();
      //CompartmentProcess::extendInterfacesOverSurface();
    }
  else
    { 
      CompartmentProcess::extendInterfacesOverSurface();
    }
}

//nRadius is the radius of MT cylinder
double MicrotubuleProcess::getDisplacementToSurface(Point& aPoint)
{
  if(nVoxelRadius >= nRadius)
    {
      return FilamentProcess::getDisplacementToSurface(aPoint);
    }
  return point2lineDist(aPoint, lengthVector, lengthStart)-nRadius;
}

void MicrotubuleProcess::removeAdjoinsFromNonBindingSide(Voxel& interface)
{
  for(unsigned i(0); i != interface.diffuseSize; ++i)
    {
      unsigned coord(interface.adjoiningCoords[i]);
      Voxel& adjoin((*theLattice)[coord]);
      if(theSpecies[getID(adjoin)]->getIsCompVacant() && 
         !theSpecies[getID(adjoin)]->getIsInterface() &&
         !isBindingSide(adjoin.coord))
        {
          /*
          if(getID(adjoin) != theSpecies[10]->getID())
            {
              theSpecies[10]->addMolecule(&adjoin);
            }
            */
          adjoin.idx = theNullID*theStride;;
        }
    }
}


//When the nVoxelRadius >= nRadius
bool MicrotubuleProcess::isWithinMTDiameter(Point& A, Point& B)
{
  /*
  const double distA(point2lineDist(A, lengthVector, lengthStart));
  const double distB(point2lineDist(B, lengthVector, lengthStart));
  return (distA <= nRadius && distB > nRadius) || 
    (distB <= nRadius && distA > nRadius);
    */
  Point iA(point2lineIntersect(A, lengthVector, lengthStart));
  Point surfNormal(sub(A, iA));
  //Point surfPoint(disp(iA, surfNormal, nRadius));
  //double surfDisp(dot(surfNormal, surfPoint));
  double surfDisp(dot(surfNormal, iA));
  const double dispA(point2planeDisp(A, surfNormal, surfDisp));
  const double dispB(point2planeDisp(B, surfNormal, surfDisp));
  return (dispA < 0) != (dispB < 0);
  /*
  Point iA(point2lineIntersect(A, lengthVector, lengthStart));
  Point surfNormal(sub(A, iA));
  Point surfPoint(disp(iA, surfNormal, nRadius));
  double surfDisp(dot(surfNormal, surfPoint));
  const double dispA(point2planeDisp(A, surfNormal, surfDisp));
  const double dispB(point2planeDisp(B, surfNormal, surfDisp));
  return (dispA < 0) != (dispB < 0);
  */
}


void MicrotubuleProcess::addSurfaceIntersectInterfaceVoxel(Voxel& aVoxel,
                                                           Point& aPoint)
{
  //CompartmentProcess::addSurfaceIntersectInterfaceVoxel(aVoxel, aPoint);
  double nearestDist(libecs::INF);
  Voxel* nearestVoxel;
  //We have already checked that aVoxel is a compVacant
  const double delta(std::max(nDiffuseRadius, nVoxelRadius)/2);
  double dispA(getDisplacementToSurface(aPoint));
  double distA(std::fabs(dispA));
  bool isInsideA(isInside(aPoint));
  if(dispA == 0 && !theSpecies[getID(aVoxel)]->getIsInterface() && isInsideA)
    {
      addInterfaceVoxel(aVoxel);
    }
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      Voxel& adjoin((*theLattice)[aVoxel.adjoiningCoords[i]]);
      Point pointB(theSpatiocyteStepper->coord2point(adjoin.coord));
      double dispB(getDisplacementToSurface(pointB));
      bool isInsideB(isInside(pointB));
      if(dispB == 0 && theSpecies[getID(adjoin)]->getIsCompVacant() &&
         !theSpecies[getID(adjoin)]->getIsInterface() && isInsideB)
        {
          addInterfaceVoxel(adjoin);
        }
      const double distB(std::fabs(dispB));
      if(!theSpecies[getID(aVoxel)]->getIsInterface() &&
         !theSpecies[getID(adjoin)]->getIsInterface() &&
         ((dispA < 0) != (dispB < 0)))
        {
          //If aVoxel is nearer to the plane:
          if(distA <= distB)
            {
              if(isInsideA)
                {
                  if(distA < nearestDist)
                    {
                      nearestDist = distA;
                      nearestVoxel = &aVoxel;
                    }
                }
              else if(theSpecies[getID(adjoin)]->getIsCompVacant() &&
                      isInsideB)
                {
                  if(distB < nearestDist && distB <= delta)
                    {
                      nearestDist = distB;
                      nearestVoxel = &adjoin;
                    }
                }
            }
          else if(theSpecies[getID(adjoin)]->getIsCompVacant())
            {
              if(isInsideB)
                {
                  if(distB < nearestDist)
                    {
                      nearestDist = distB;
                      nearestVoxel = &adjoin;
                    }
                }
              else if(isInsideA)
                {
                  if(distA < nearestDist && distA <= delta)
                    {
                      nearestDist = distA;
                      nearestVoxel = &aVoxel;
                    }
                }
            }
        }
    }
  if(nearestDist != libecs::INF)
    { 
      addInterfaceVoxel(*nearestVoxel);
    }
}


}
