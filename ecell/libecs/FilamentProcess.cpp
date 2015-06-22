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

#include <libecs/FilamentProcess.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(FilamentProcess, Process); 

void FilamentProcess::prepreinitialize() {
  SpatiocyteProcess::prepreinitialize();
  theInterfaceVariable = createVariable("Interface");
}

void FilamentProcess::initialize() {
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
                                  "]: This compartment requires only " +
                                  "one vacant variable reference with -1 " +
                                  "coefficient as the vacant species of " +
                                  "the compartment, but " +
                                  getIDString(theVacantSpecies) + " and " +
                                  getIDString(aSpecies) + " are given."); 
                }
              theVacantSpecies = aSpecies;
            }
          else if((*i).getCoefficient() == -2)
            {
              if(theMinusSpecies)
                {
                  THROW_EXCEPTION(ValueError, String(
                                  getPropertyInterface().getClassName()) +
                                  "[" + getFullID().asString() + 
                                  "]: This compartment requires only " +
                                  "one variable reference with -2 " +
                                  "coefficient as the minus end species " +
                                  "of the compartment, but " +
                                  getIDString(theMinusSpecies) + " and " +
                                  getIDString(aSpecies) + " are given."); 
                }
              theMinusSpecies = aSpecies;
            }
          else if((*i).getCoefficient() == -3)
            {
              if(thePlusSpecies)
                {
                  THROW_EXCEPTION(ValueError, String(
                                  getPropertyInterface().getClassName()) +
                                  "[" + getFullID().asString() + 
                                  "]: This compartment requires only " +
                                  "one variable reference with -3 " +
                                  "coefficient as the plus end species " +
                                  "of the microtubule compartment, but " +
                                  getIDString(thePlusSpecies) + " and " +
                                  getIDString(aSpecies) + " are given."); 
                }
              thePlusSpecies = aSpecies;
            }
        }
      else
        {
          theBindingSpecies.push_back(aSpecies);
        }
    }
  if(!theBindingSpecies.size())
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + 
                      "]: This compartment requires at least one " +
                      "nonHD variable reference with zero coefficient " +
                      "as the binding species, but none is given."); 
    }
  if(!theVacantSpecies)
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + 
                      "]: This compartment requires one " +
                      "nonHD variable reference with negative " +
                      "coefficient as the vacant species, " +
                      "but none is given."); 
    }
  if(!theMinusSpecies)
    {
      theMinusSpecies = theVacantSpecies;
    }
  if(!thePlusSpecies)
    {
      thePlusSpecies = theVacantSpecies;
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
  nRadius = Radius/(VoxelRadius*2);
  nGridSize = 10*nDiffuseRadius;
}

void FilamentProcess::initializeFirst() {
  CompartmentProcess::initializeFirst();
  theMinusSpecies->setIsOffLattice();
  theMinusSpecies->setComp(theComp);
  thePlusSpecies->setIsOffLattice();
  thePlusSpecies->setComp(theComp);
  for(unsigned i(0); i != theBindingSpecies.size(); ++i)
    {
      theBindingSpecies[i]->setIsOffLattice();
      theBindingSpecies[i]->setDimension(1);
      theBindingSpecies[i]->setVacantSpecies(theVacantSpecies);
      theBindingSpecies[i]->setComp(theComp);
      theBindingSpecies[i]->resetFixedAdjoins();
    }
}

unsigned FilamentProcess::getLatticeResizeCoord(unsigned aStartCoord) {
  const unsigned aSize(CompartmentProcess::getLatticeResizeCoord(aStartCoord));
  theMinusSpecies->resetFixedAdjoins();
  theMinusSpecies->setMoleculeRadius(DiffuseRadius);
  thePlusSpecies->resetFixedAdjoins();
  thePlusSpecies->setMoleculeRadius(DiffuseRadius);
  for(unsigned i(0); i != theBindingSpecies.size(); ++i)
    {
      theBindingSpecies[i]->setMoleculeRadius(DiffuseRadius);
    }
  return aSize;
}

void FilamentProcess::setCompartmentDimension() {
  if(Length)
    {
      Subunits = (unsigned)rint(Length/(2*DiffuseRadius));
    }
  Length = Subunits*2*DiffuseRadius;
  Width = Radius*2;
  Height = Radius*2;
  theDimension = 1;
  /*
  Origin.x += OriginX*theComp->lengthX/2;
  Origin.y += OriginY*theComp->lengthY/2;
  Origin.z += OriginZ*theComp->lengthZ/2;
  */
  allocateGrid();
}

void FilamentProcess::initializeCompartment() {
  if(!isCompartmentalized)
    {
      thePoints.resize(endCoord-subStartCoord);
      vacStartIndex = theVacantSpecies->size();
      intStartIndex = theInterfaceSpecies->size();
      initializeVectors();
      initializeFilaments(subunitStart, Filaments, Subunits, 0,
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

void FilamentProcess::setTrailSize(unsigned start, unsigned end) {
  for(unsigned i(start); i != end; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
      subunit.trailSize = subunit.adjoiningSize;
    }
}

void FilamentProcess::initializeVectors() { 
  //Minus end
  Minus.x = -nLength/2;
  Minus.y = 0;
  Minus.z = 0;
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  Point tmpOrigin;
  tmpOrigin.x = OriginX*aComp->lengthX/2;
  tmpOrigin.y = OriginY*aComp->lengthY/2;
  tmpOrigin.z = OriginZ*aComp->lengthZ/2;
  //Rotated Minus end
  theSpatiocyteStepper->rotateX(theComp->rotateX, &Minus, -1);
  theSpatiocyteStepper->rotateY(theComp->rotateY, &Minus, -1);
  theSpatiocyteStepper->rotateZ(theComp->rotateZ, &Minus, -1);
  theSpatiocyteStepper->rotateX(RotateX, &Minus, 1);
  theSpatiocyteStepper->rotateY(RotateY, &Minus, 1);
  theSpatiocyteStepper->rotateZ(RotateZ, &Minus, 1);
  theSpatiocyteStepper->rotateX(theComp->rotateX, &tmpOrigin, -1);
  theSpatiocyteStepper->rotateY(theComp->rotateY, &tmpOrigin, -1);
  theSpatiocyteStepper->rotateZ(theComp->rotateZ, &tmpOrigin, -1);
  add_(tmpOrigin, Origin);
  add_(Minus, tmpOrigin);
  //Direction vector from the Minus end to center
  lengthVector = sub(tmpOrigin, Minus);
  //Make direction vector a unit vector
  norm_(lengthVector);
  //Rotated Plus end
  Plus = disp(Minus, lengthVector, nLength);
  setSubunitStart();
}

void FilamentProcess::setSubunitStart()
{
  subunitStart = Minus;
}

void FilamentProcess::initializeFilaments(Point& aStartPoint, unsigned aRows,
                                          unsigned aCols, double aRadius,
                                          Species* aVacant,
                                          unsigned aStartCoord) {
  Voxel* aVoxel(addCompVoxel(0, 0, aStartPoint, aVacant, aStartCoord, aCols));
  theMinusSpecies->addMolecule(aVoxel);
}

// y:width:rows:filaments
// z:length:cols:subunits
void FilamentProcess::connectFilaments(unsigned aStartCoord,
                                          unsigned aRows, unsigned aCols) {
  for(unsigned i(0); i != aCols; ++i)
    {
      for(unsigned j(0); j != aRows; ++j)
        { 
          if(i > 0)
            { 
              //NORTH-SOUTH
              unsigned a(aStartCoord+j*aCols+i);
              unsigned b(aStartCoord+j*aCols+(i-1));
              connectSubunit(a, b, NORTH, SOUTH);
            }
          else if(Periodic)
            {
              //periodic NORTH-SOUTH
              unsigned a(aStartCoord+j*aCols); 
              unsigned b(aStartCoord+j*aCols+aCols-1);
              connectSubunit(a, b, NORTH, SOUTH);
            }
        }
    }
}

// y:width:rows:filaments
// z:length:cols:subunits
void FilamentProcess::connectTrailSubunits(unsigned aStartCoord,
                                              unsigned aRows, unsigned aCols)
{
  for(unsigned i(0); i != aCols; ++i)
    {
      for(unsigned j(1); j != aRows; ++j)
        { 
          //NW-SW
          unsigned a(aStartCoord+j*aCols+i);
          unsigned b(aStartCoord+(j-1)*aCols+i);
          connectSubunit(a, b, NW, SW);
        }
    }
}

void FilamentProcess::elongateFilaments(Species* aVacant,
                                           unsigned aStartCoord,
                                           unsigned aRows,
                                           unsigned aCols,
                                           double aRadius)
{
  for(unsigned i(0); i != aRows; ++i)
    {
      Voxel* startVoxel(&(*theLattice)[aStartCoord+i*aCols]);
      Point A(*startVoxel->point);
      for(unsigned j(1); j != aCols; ++j)
        {
          disp_(A, lengthVector, aRadius*2);
          Voxel* aVoxel(addCompVoxel(i, j, A, aVacant, aStartCoord, aCols));
          if(j == aCols-1)
            {
              thePlusSpecies->addMolecule(aVoxel);
            }
        }
    }
}

//Is inside the parent compartment and confined by the length of the MT:
bool FilamentProcess::isInside(Point& aPoint)
{
  double disp(point2planeDisp(aPoint, lengthVector, dot(lengthVector, Minus)));
  //Use nDifffuseRadius/2 instead of 0 because we don't want additional
  //interface voxels at the edge of the plus or minus end. So only molecules
  //hitting along the normal of the surface of the MT at the ends can bind.
  //This would avoid bias of molecule directly hitting the MT ends from the
  //sides and binding:
  if(disp > nDiffuseRadius/2)
    { 
      disp = point2planeDisp(aPoint, lengthVector, dot(lengthVector, Plus));
      if(disp < -nDiffuseRadius/2)
        {
          return true;
        }
    }
  return false;
}

bool FilamentProcess::isOnAboveSurface(Point& aPoint)
{
  double disp(point2lineDist(aPoint, lengthVector, Minus));
  if(disp >= nRadius)
    {
      return true;
    }
  return false;
}

double FilamentProcess::getDistanceToSurface(Point& aPoint)
{
  return abs(point2lineDist(aPoint, lengthVector, Minus)-nRadius);
}

void FilamentProcess::extendInterfacesOverSurface()
{
  for(unsigned i(intStartIndex); i != theInterfaceSpecies->size(); ++i)
    {
      //Traverse both directions of the first interface:
      if(i == intStartIndex+1 && theInterfaceSpecies->size()-intStartIndex == 2)
        {
          --i;
        }
      unsigned voxelCoord(theInterfaceSpecies->getCoord(i));
      Voxel& anInterface((*theLattice)[voxelCoord]);
      Point intPoint(theSpatiocyteStepper->coord2point(anInterface.coord));
      Point intLinePoint(point2lineIntersect(intPoint, lengthVector, Minus));
      double aSurfaceDisplace(dot(lengthVector, intLinePoint));
      double nearestAdjDist(libecs::INF);
      Voxel* nearestAdj;
      for(unsigned j(0); j != theAdjoiningCoordSize; ++j)
        {
          Voxel& adjoin((*theLattice)[anInterface.adjoiningCoords[j]]);
          if(theSpecies[getID(adjoin)]->getIsCompVacant())
            {
              Point adjPoint(theSpatiocyteStepper->coord2point( adjoin.coord));
              Point adjLinePoint(point2lineIntersect(adjPoint, lengthVector,
                                                     Minus));
              double adjDist(point2lineDist(adjPoint, lengthVector, Minus));
              double adjLineDist(distance(intLinePoint, adjLinePoint));
              double adjDisp(point2planeDisp(adjPoint, lengthVector,
                                             aSurfaceDisplace));
              for(unsigned k(0); k != theAdjoiningCoordSize; ++k)
                {
                  Voxel& subAdjoin((*theLattice)[adjoin.adjoiningCoords[k]]);
                  if(theSpecies[getID(subAdjoin)]->getIsCompVacant())
                    { 
                      Point subPoint(theSpatiocyteStepper->coord2point(
                                                             subAdjoin.coord));
                      Point subLinePoint(point2lineIntersect(subPoint,
                                                         lengthVector, Minus));
                      double subDist(point2lineDist(subPoint, lengthVector,
                                                    Minus));
                      double subLineDist(distance(intLinePoint, subLinePoint));
                      double subDisp(point2planeDisp(subPoint, lengthVector,
                                             aSurfaceDisplace));
                      if(subDist+adjDist < nearestAdjDist && 
                         isInside(adjPoint) && 
                         (adjDisp < 0) == (subDisp < 0) && 
                         subLineDist > adjLineDist)
                        {
                          nearestAdjDist = subDist+adjDist;
                          nearestAdj = &adjoin;
                        }
                    }
                }
            }
        }
      if(nearestAdjDist <= nDiffuseRadius*2 && 
         getID(nearestAdj) != theInterfaceSpecies->getID())
        { 
          addInterfaceVoxel(*nearestAdj);
        }
    }
}


void FilamentProcess::addFilamentIntersectInterfaceVoxel(Voxel& aVoxel,
                                                         Point& aPoint,
                                                         Point& filamentPoint)
{
  //std::cout << "inside add surface, nRadius:" << nRadius << std::endl;
  Point aLinePoint(point2lineIntersect(aPoint, lengthVector, Minus));
  Point aSurfaceNormal(sub(aPoint, aLinePoint));
  double aSurfaceDisplace(dot(aSurfaceNormal, Minus));
  //We have already checked that aVoxel is a compVacant
  double dispA(point2planeDisp(aPoint, aSurfaceNormal, aSurfaceDisplace));
  //std::cout << "dispA:" << dispA << std::endl;
  if(dispA == 0 && getID(aVoxel) != theInterfaceSpecies->getID())
    {
      addInterfaceVoxel(aVoxel);
    }
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      Voxel& adjoin((*theLattice)[aVoxel.adjoiningCoords[i]]);
      Point pointB(theSpatiocyteStepper->coord2point(adjoin.coord));
      double dispB(point2planeDisp(pointB, aSurfaceNormal, aSurfaceDisplace));
      //std::cout << "dispB:" << dispB << std::endl;
      if(dispB == 0 && theSpecies[getID(adjoin)]->getIsCompVacant() &&
         getID(adjoin) != theInterfaceSpecies->getID())
        {
          addInterfaceVoxel(adjoin);
        }
      if((dispA < 0) != (dispB < 0) &&
         getID(aVoxel) != theInterfaceSpecies->getID() &&
         getID(adjoin) != theInterfaceSpecies->getID())
        {
          //If aVoxel is nearer to the plane:
          if(abs(dispA) <= abs(dispB))
            {
              //std::cout << "dispA" << std::endl;
              addInterfaceVoxel(aVoxel);
            }
          //If the adjoin is nearer to the plane:
          else if(theSpecies[getID(adjoin)]->getIsCompVacant())
            {
              //std::cout << "dispB" << std::endl;
              addInterfaceVoxel(adjoin);
            }
        }
    }
}

/*
void FilamentProcess::addSurfaceIntersectInterfaceVoxel(Voxel& aVoxel,
                                                        Point& aPoint)
{
  std::cout << "inside add surface, nRadius:" << nRadius << std::endl;
  Point aLinePoint(point2lineIntersect(aPoint, lengthVector, Minus));
  Point aSurfaceNormal(sub(aPoint, aLinePoint));
  double aSurfaceDisplace(dot(aSurfaceNormal, Minus));
  //We have already checked that aVoxel is a compVacant
  double dispA(point2planeDisp(aPoint, aSurfaceNormal, aSurfaceDisplace));
  std::cout << "dispA:" << dispA << std::endl;
  if(dispA == 0 && getID(aVoxel) != theInterfaceSpecies->getID())
    {
      addInterfaceVoxel(aVoxel);
    }
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      Voxel& adjoin((*theLattice)[aVoxel.adjoiningCoords[i]]);
      Point pointB(theSpatiocyteStepper->coord2point(adjoin.coord));
      double dispB(point2planeDisp(pointB, aSurfaceNormal, aSurfaceDisplace));
      std::cout << "dispB:" << dispB << std::endl;
      if(dispB == 0 && theSpecies[getID(adjoin)]->getIsCompVacant() &&
         getID(adjoin) != theInterfaceSpecies->getID())
        {
          addInterfaceVoxel(adjoin);
        }
      if((dispA < 0) != (dispB < 0) &&
         getID(aVoxel) != theInterfaceSpecies->getID() &&
         getID(adjoin) != theInterfaceSpecies->getID())
        {
          //If aVoxel is nearer to the plane:
          if(abs(dispA) <= abs(dispB))
            {
              std::cout << "dispA" << std::endl;
              addInterfaceVoxel(aVoxel);
            }
          //If the adjoin is nearer to the plane:
          else if(theSpecies[getID(adjoin)]->getIsCompVacant())
            {
              std::cout << "dispB" << std::endl;
              addInterfaceVoxel(adjoin);
            }
        }
    }
}
*/

}
