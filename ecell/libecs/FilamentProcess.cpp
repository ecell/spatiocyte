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
  setSubunitStart();
  if(Length)
    {
      Subunits = (unsigned)rint(Length/(2*DiffuseRadius));
    }
  //Length is from the center point of the first subunit to the center point
  //of the last subunit
  Length = Subunits*2*DiffuseRadius-2*DiffuseRadius;
  Width = Radius*2;
  Height = Radius*2;
  theDimension = 1;
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
  //theSpecies[2]->setIsPopulated();
}

void FilamentProcess::setTrailSize(unsigned start, unsigned end) {
  for(unsigned i(start); i != end; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
      subunit.trailSize = subunit.adjoiningSize;
    }
}

/*
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
*/

void FilamentProcess::initializeVectors() { 
  //subunitStart is the center point of the first vacant species voxel:
  lengthStart = subunitStart;
  
  Point& origin(theComp->centerPoint);
  Point tmp(sub(subunitStart, origin));
  rotate(tmp);
  subunitStart = add(tmp, origin);

  tmp = sub(lengthStart, origin);
  rotate(tmp);
  lengthStart = add(tmp, origin);
  Minus = lengthStart;

  rotate(lengthVector);
  lengthEnd = disp(lengthStart, lengthVector, nLength);
  Plus = lengthEnd;

  rotate(widthVector);
  widthEnd = disp(lengthEnd, widthVector, nWidth);

  rotate(heightVector);
  heightEnd = disp(widthEnd, heightVector, nHeight);

  Point center(lengthStart);
  disp_(center, lengthVector, nLength/2);
  disp_(center, widthVector, nWidth/2);
  theComp->centerPoint = center;

  heightDisplace = dot(heightVector, lengthStart);
  lengthDisplace = dot(lengthVector, lengthStart);
  lengthDisplaceOpp = dot(lengthVector, lengthEnd);
  widthDisplace = dot(widthVector, lengthStart);
  
  //remove the following
  surfaceNormal = heightVector;
  surfaceDisplace = heightDisplace;
}


void FilamentProcess::setSubunitStart()
{
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  //row => z
  //col => x
  //layer => y
  const unsigned minCoord(theSpatiocyteStepper->global2coord(
                                          aComp->minCoord.row,
                                          aComp->minCoord.layer,
                                          aComp->minCoord.col));
  const unsigned maxCoord(theSpatiocyteStepper->global2coord(
                                          aComp->maxCoord.row,
                                          aComp->maxCoord.layer,
                                          aComp->maxCoord.col));
  Point maxPoint(theSpatiocyteStepper->coord2point(maxCoord));
  subunitStart = theSpatiocyteStepper->coord2point(minCoord);
 /* 
  Point maxPoint(aComp->maxPoint);
  subunitStart = aComp->minPoint;
  */
  Point lengths(sub(maxPoint, subunitStart));
  Point& center(theComp->centerPoint);
  center.x = 0;
  center.y = 0;
  center.z = 0;
  if(OriginX == 1)
    {
      subunitStart.x = maxPoint.x;
    }
  else if(OriginX != -1)
    {
      subunitStart.x += lengths.x*0.5*(1+OriginX);
    }
  if(OriginY == 1)
    {
      subunitStart.y = maxPoint.y;
    }
  else if(OriginY != -1)
    {
      subunitStart.y += lengths.y*0.5*(1+OriginY);
    }
  if(OriginZ == 1)
    {
      subunitStart.z = maxPoint.z;
    }
  else if(OriginZ != -1)
    {
      subunitStart.z += lengths.z*0.5*(1+OriginZ);
    }
  if(LineZ)
    {
      lengthVector = Point(0, 0, 1);
      widthVector = Point(1, 0, 0);
      heightVector = Point(0, 1, 0);
      if(!Length && !Subunits)
        {
          Length = lengths.z*VoxelRadius*2;
          center.z = lengths.z*0.5;
        }
      subunitStart.z -= lengths.z*0.5;
    }
  else if(LineY)
    {
      lengthVector = Point(0, 1, 0);
      widthVector = Point(1, 0, 0);
      heightVector = Point(0, 0, 1);
      if(!Length && !Subunits)
        {
          Length = lengths.y*VoxelRadius*2;
          center.y = lengths.y*0.5;
        }
      subunitStart.y -= lengths.y*0.5;
    }
  //LineX 
  else
    {
      lengthVector = Point(1, 0, 0);
      widthVector = Point(0, 0, 1);
      heightVector = Point(0, 1, 0);
      if(!Length && !Subunits)
        {
          Length = lengths.x*VoxelRadius*2;
          center.x = lengths.x*0.5;
        }
      subunitStart.x -= lengths.x*0.5;
    }
  add_(center, subunitStart);
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
  //Use nDifffuseRadius/2 instead of 0 because we don't want additional
  //interface voxels at the edge of the plus or minus end. So only molecules
  //hitting along the normal of the surface of the MT at the ends can bind.
  //This would avoid bias of molecule directly hitting the MT ends from the
  //sides and binding:
  double dispA(point2planeDisp(aPoint, lengthVector, lengthDisplace));
  if(dispA >= 0 && dispA <= nLength)
    {
      return true;
    }
  return false;
}


//Returns positive distance if the point is within the Minus and Plus
//points of the filament, otherwise returns the distance overshot the Minus
//or Plus points in negative.
double FilamentProcess::getMinDistanceFromLineEnd(Point& aPoint)
{
  double dispA(point2planeDisp(aPoint, lengthVector, lengthDisplace));
  return std::min(dispA, nLength-dispA);
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

unsigned FilamentProcess::getAdjoiningInterfaceCnt(Voxel& aVoxel)
{
  unsigned cnt(0);
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      const unsigned coord(aVoxel.adjoiningCoords[i]);
      Voxel& adjoin((*theLattice)[coord]);
      if(getID(adjoin) == theInterfaceSpecies->getID())
        {
          ++cnt;
        }
    }
  return cnt;
}

bool FilamentProcess::isAdjoin(Voxel* aSource, Voxel* aTarget)
{
  for(unsigned i(0); i != theAdjoiningCoordSize; ++i)
    {
      const unsigned coord(aSource->adjoiningCoords[i]);
      if(aTarget == &(*theLattice)[coord])
        {
          return true;
        }
    }
  return false;
}

bool FilamentProcess::getFilamentAdjoin(Voxel* aVoxel,
                                        const bool direction,
                                        const unsigned adjIndex,
                                        const unsigned maxAdjInterface,
                                        const double aSurfaceDisp,
                                        const double aDisp,
                                        Point& adjPoint,
                                        double& adjDist,
                                        double& adjDisp,
                                        Voxel** adjoin)
{
  *adjoin = &(*theLattice)[aVoxel->adjoiningCoords[adjIndex]];
  if(theSpecies[getID(*adjoin)]->getIsCompVacant() &&
     getID(*adjoin) != theInterfaceSpecies->getID())
    { 
      adjPoint = theSpatiocyteStepper->coord2point((*adjoin)->coord);
      adjDist = point2lineDist(adjPoint, lengthVector, Minus);
      adjDisp = point2planeDisp(adjPoint, lengthVector, aSurfaceDisp);
      if((adjDisp < 0) == direction && abs(adjDisp) > std::max(abs(aDisp), 0.2)
         && getAdjoiningInterfaceCnt(**adjoin) <= maxAdjInterface)
        {
          return true;
        }
    }
  return false;
}

void FilamentProcess::extendInterfacesOverSurface()
{
  const double delta(1.6*nDiffuseRadius);
  bool direction(true);
  for(int i(intStartIndex); i != theInterfaceSpecies->size(); ++i)
    {
      unsigned intLatticeCoord(theInterfaceSpecies->getCoord(i));
      Voxel& interface((*theLattice)[intLatticeCoord]);
      if(getAdjoiningInterfaceCnt(interface) > 1)
        {
          continue;
        }
      Point intPoint(theSpatiocyteStepper->coord2point(interface.coord));
      Point intLinePoint(point2lineIntersect(intPoint, lengthVector, Minus));
      double aSurfaceDisp(dot(lengthVector, intLinePoint));
      double firstDist(libecs::INF);
      double secondDist(libecs::INF);
      double thirdDist(libecs::INF);
      Voxel* firstAdj(NULL);
      Voxel* secondAdj(NULL);
      Voxel* thirdAdj(NULL);
      Voxel* secondSub(NULL);
      Voxel* thirdSub(NULL);
      Voxel* thirdSubSub(NULL);
      for(unsigned j(0); j != theAdjoiningCoordSize; ++j)
        {
          Voxel* adjoin;
          double adjDist;
          double adjDisp;
          Point adjPoint;
          if(getFilamentAdjoin(&interface, direction, j, 1, aSurfaceDisp,
                               0.2, adjPoint, adjDist, adjDisp, &adjoin))
            { 
              if(getMinDistanceFromLineEnd(adjPoint) > -nDiffuseRadius)
                {
                  for(unsigned k(0); k != theAdjoiningCoordSize; ++k)
                    {
                      Voxel* sub;
                      double subDist;
                      double subDisp;
                      Point subPoint;
                      if(getFilamentAdjoin(adjoin, direction, k, 0,
                         aSurfaceDisp, adjDisp, subPoint, subDist, subDisp,
                         &sub))
                        {
                          if(getMinDistanceFromLineEnd(subPoint) > 
                             -nDiffuseRadius)
                            {
                              for(unsigned l(0); l != theAdjoiningCoordSize;
                                  ++l)
                                {
                                  Voxel* subSub;
                                  double subSubDist;
                                  double subSubDisp;
                                  Point subSubPoint;
                                  if(getFilamentAdjoin(sub, direction, l, 0,
                                     aSurfaceDisp, subDisp, subSubPoint,
                                     subSubDist, subSubDisp, &subSub))
                                    {
                                      if(!isAdjoin(subSub, adjoin) &&
                                         getMinDistanceFromLineEnd(subSubPoint)
                                         > nDiffuseRadius)
                                        {
                                          if(subSubDist+subDist+adjDist < 
                                             thirdDist && 
                                             adjDist < delta &&
                                             subDist < delta &&
                                             subSubDist < delta)
                                            {
                                              thirdDist = 
                                                subSubDist+subDist+adjDist;
                                              thirdAdj = adjoin;
                                              thirdSub = sub;
                                              thirdSubSub = subSub;
                                            }
                                        }
                                    }
                                }
                              if(adjDist+subDist < secondDist &&
                                 adjDist < delta &&
                                 subDist < delta)
                                {
                                  secondDist = adjDist+subDist;
                                  secondAdj = adjoin;
                                  secondSub = sub;
                                }
                            }
                        }
                    }
                  if(adjDist < firstDist && adjDist < delta)
                    {
                      firstDist = adjDist;
                      firstAdj = adjoin;
                    }
                }
            }
        }
      if(thirdDist != libecs::INF)
        {
          addInterfaceVoxel(*thirdAdj);
          addInterfaceVoxel(*thirdSub);
          //addInterfaceVoxel(*thirdSubSub);
        }
      else if(secondDist != libecs::INF)
        {
          addInterfaceVoxel(*secondAdj);
          addInterfaceVoxel(*secondSub);
        }
      else if(firstDist != libecs::INF)
        {
          addInterfaceVoxel(*firstAdj);
        }
      else if(i == theInterfaceSpecies->size()-1)
        {
          direction = !direction;
          i = intStartIndex-1;
          if(direction == true)
            {
              return;
            }
        }
    }
}


}
