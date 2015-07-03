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

#include <libecs/CompartmentProcess.hpp>
#include <libecs/SpatiocyteVector.hpp>

LIBECS_DM_INIT_STATIC(CompartmentProcess, Process); 


void CompartmentProcess::prepreinitialize()
{
  SpatiocyteProcess::prepreinitialize();
  theInterfaceVariable = createVariable("Interface");
  if(LipidRadius)
    {
      if(LipidRadius < 0)
        {
          LipidRadius = 0;
        }
      else
        {
          theLipidVariable = createVariable("Lipid");
        }
    }
}

void CompartmentProcess::initialize()
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
      //HD Species
      if(aSpecies == NULL)
        {
          continue;
        }
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
          else
            {
              theLipidCompSpecies.push_back(aSpecies);
            }
        }
      else
        {
          theVacantCompSpecies.push_back(aSpecies);
        }
      if(RegularLattice)
        {
          aSpecies->setIsRegularLattice(theDiffuseSize);
        }
      if(Periodic)
        {
          aSpecies->setIsPeriodic();
        }
    }
  if(!theVacantSpecies)
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + 
                      "]: This compartment requires one " +
                      "nonHD variable reference with -1 " +
                      "coefficient as the vacant species, " +
                      "but none is given."); 
    }
  if(LipidRadius)
    {
      theLipidSpecies = theSpatiocyteStepper->addSpecies(
                                                   theLipidVariable);
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
  //Lattice voxel radius:
  VoxelRadius = theSpatiocyteStepper->getVoxelRadius();
  //Normalized off-lattice voxel radius:
  nDiffuseRadius = DiffuseRadius/(VoxelRadius*2);
  //SubunitRadius is the actual radius of a molecule.
  //For a normal molecule the SubunitRadius = DiffuseRadius.
  //For a multiscale molecule, the SubunitRadius can be larger
  //than the DiffuseRadius:
  nSubunitRadius = SubunitRadius/(VoxelRadius*2);
  //Normalized lipid voxel radius:
  nLipidRadius = LipidRadius/(VoxelRadius*2);
  nGridSize = 10*nDiffuseRadius;
}

void CompartmentProcess::initializeFirst()
{
  SpatiocyteProcess::initializeFirst();
  theComp = new Comp;
  theVacantSpecies->setIsCompVacant();
  theVacantSpecies->setIsOffLattice();
  theVacantSpecies->setComp(theComp);
  if(theLipidSpecies)
    {
      theLipidSpecies->setIsCompVacant();
      theLipidSpecies->setIsOffLattice();
      theLipidSpecies->setComp(theComp);
    }
  for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
    {
      theLipidCompSpecies[i]->setIsOffLattice();
      //setVacantSpecies must be declared here since it needs
      //to be overwritten by DiffusionProcess in initializeSecond:
      theLipidCompSpecies[i]->setVacantSpecies(theLipidSpecies);
      theLipidCompSpecies[i]->setComp(theComp);
    }
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      theVacantCompSpecies[i]->setIsOffLattice();
      //setVacantSpecies must be declared here since it needs
      //to be overwritten by DiffusionProcess in initializeSecond:
      theVacantCompSpecies[i]->setVacantSpecies(theVacantSpecies);
      theVacantCompSpecies[i]->setComp(theComp);
      if(theLipidSpecies)
        {
          theVacantCompSpecies[i]->setIsMultiscale();
        }
    }
}

unsigned CompartmentProcess::getLatticeResizeCoord(unsigned aStartCoord)
{
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  if(aComp->diffusiveComp)
    {
      Comp* aDiffusiveComp(aComp->diffusiveComp);
      if(aDiffusiveComp->interfaceID == theSpecies.size())
        {
          aDiffusiveComp->interfaceID = theInterfaceSpecies->getID();
        }
      else
        {
          theInterfaceSpecies = theSpecies[aDiffusiveComp->interfaceID];
        }
    }
  aComp->interfaceID = theInterfaceSpecies->getID();
  *theComp = *aComp;
  theVacantSpecies->resetFixedAdjoins();
  theVacantSpecies->setMoleculeRadius(DiffuseRadius);
  if(theLipidSpecies)
    {
      theLipidSpecies->resetFixedAdjoins();
      theLipidSpecies->setMoleculeRadius(LipidRadius);
    }
  //The compartment center point (origin):
  Origin = aComp->centerPoint;
  setCompartmentDimension();
  theComp->dimension = theDimension;
  theVacantSpecies->setDimension(theDimension);
  setLipidCompSpeciesProperties();
  setVacantCompSpeciesProperties();
  subStartCoord = aStartCoord;
  lipStartCoord = aStartCoord+Filaments*Subunits;
  endCoord = lipStartCoord+LipidRows*LipidCols;
  return endCoord-aStartCoord;
}

void CompartmentProcess::setVacantCompSpeciesProperties()
{
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      theVacantCompSpecies[i]->setDimension(theDimension);
      if(SubunitRadius && 
         theVacantCompSpecies[i]->getMoleculeRadius() == VoxelRadius)
        {
          theVacantCompSpecies[i]->setMoleculeRadius(SubunitRadius);
        }
      theVacantCompSpecies[i]->setDiffuseRadius(DiffuseRadius);
      if(theLipidSpecies)
        {
          theVacantCompSpecies[i]->setMultiscaleVacantSpecies(theLipidSpecies);
        }
    }
}

int CompartmentProcess::getCoefficient(Species* aSpecies)
{
  for(VariableReferenceVector::iterator i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      if(aSpecies->getVariable() == (*i).getVariable()) 
        {
          return (*i).getCoefficient();
        }
    }
  return 0;
}

Species* CompartmentProcess::coefficient2species(int aCoeff)
{
  for(VariableReferenceVector::iterator i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      if((*i).getCoefficient() == aCoeff)
        {
          return theSpatiocyteStepper->variable2species((*i).getVariable());
        }
    }
  return NULL;
}

void CompartmentProcess::setLipidCompSpeciesProperties()
{
  for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
    {
      theLipidCompSpecies[i]->setDimension(theDimension);
      theLipidCompSpecies[i]->setMoleculeRadius(LipidRadius);
    }
}

void CompartmentProcess::updateResizedLattice()
{
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      //TODO: replace subStartCoord with theVacantSpecies->getCoord(0) 
      //TODO: replace lipStartCoord with theLipidSpecies->getCoord(0) 
      theVacantCompSpecies[i]->setVacStartCoord(subStartCoord, Filaments,
                                                Subunits);
      if(theLipidSpecies)
        {
          theVacantCompSpecies[i]->setLipStartCoord(lipStartCoord, LipidRows,
                                                    LipidCols);
        }
      //To support walkRegular in SpatiocyteSpecies
      else if(RegularLattice)
        {
          theVacantCompSpecies[i]->setLipStartCoord(subStartCoord, Filaments,
                                                    Subunits);
        }
    }
  for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
    {
      theLipidCompSpecies[i]->setLipStartCoord(lipStartCoord, LipidRows,
                                                LipidCols);
    }
}


//May not work well when parent compartment is rotated. Modify this according
//to FilamentProcess::setSubunitStart if it doesn'w work well.
void CompartmentProcess::setSubunitStart()
{
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  lengthVector = Point(0, 0, 1);
  widthVector = Point(0, 1, 0);
  heightVector = Point(1, 0, 0);
  if(Autofit && aComp->surfaceSub)
    {
      Point nearest;
      Point farthest;
      getStartVoxelPoint(subunitStart, nearest, farthest);
      double dist(subunitStart.z-nearest.z+nVoxelRadius-nDiffuseRadius);
      subunitStart.z -= int(dist/(nDiffuseRadius*2))*2*nDiffuseRadius;
      dist = subunitStart.y-nearest.y+nVoxelRadius-nDiffuseRadius;
      unsigned cnt(int(dist/(nDiffuseRadius*sqrt(3))));
      subunitStart.y -= cnt*nDiffuseRadius*sqrt(3);
      if(cnt%2 == 1)
        {
          subunitStart.z += nDiffuseRadius;
        }
      Width = (farthest.y-subunitStart.y+nVoxelRadius+nDiffuseRadius)*
        VoxelRadius*2.00001;
      Length = (farthest.z-subunitStart.z+nVoxelRadius+nDiffuseRadius)*
        VoxelRadius*2.00001;
    }
  else
    {
      //row => z
      //col => x
      //layer => y
      //subunitStart = aComp->minPoint;
      const unsigned coord(theSpatiocyteStepper->global2coord(
                                              aComp->minCoord.row,
                                              aComp->minCoord.layer,
                                              aComp->minCoord.col));
      subunitStart = theSpatiocyteStepper->coord2point(coord);
      Point& center(theComp->centerPoint);
      if(PlaneYZ)
        {
          if(!Length && !Subunits)
            {
              Length = aComp->lengthZ*VoxelRadius*2;
              center.z = aComp->lengthZ;
            }
          if(!Width && !Filaments)
            {
              Width = aComp->lengthY*VoxelRadius*2;
              center.y = aComp->lengthY;
            }
          center.x = 2*nDiffuseRadius; 
          subunitStart.x = aComp->maxPoint.x;
        }
      else if(PlaneXZ)
        {
          lengthVector = Point(0, 0, 1);
          widthVector = Point(1, 0, 0);
          heightVector = Point(0, 1, 0);
          if(!Length && !Subunits)
            {
              Length = aComp->lengthZ*VoxelRadius*2;
              center.z = aComp->lengthZ;
            }
          if(!Width && !Filaments)
            {
              Width = aComp->lengthX*VoxelRadius*2;
              center.x = aComp->lengthX;
            }
          center.y = 2*nDiffuseRadius;
          subunitStart.y = aComp->maxPoint.y;
        }
      else
        {
          lengthVector = Point(1, 0, 0);
          widthVector = Point(0, 1, 0);
          heightVector = Point(0, 0, 1);
          if(!Length && !Subunits)
            {
              Length = aComp->lengthX*VoxelRadius*2;
              center.x = aComp->lengthX;
            }
          if(!Width && !Filaments)
            {
              Width = aComp->lengthY*VoxelRadius*2;
              center.y = aComp->lengthY;
            }
          center.z = 2*nDiffuseRadius;
          subunitStart.z = aComp->maxPoint.z;
        }
      center.x = center.x/2+subunitStart.x;
      center.y = center.y/2+subunitStart.y;
      center.z = center.z/2+subunitStart.z;
      const double multX(OriginX*theComp->lengthX/2);
      const double multY(OriginY*theComp->lengthY/2);
      const double multZ(OriginZ*theComp->lengthZ/2);
      center.x += multX;
      center.y += multY;
      center.z += multZ;
      subunitStart.x += multX;
      subunitStart.y += multY;
      subunitStart.z += multZ;
    }
}

// width:rows:filaments
// length:cols:subunits
void CompartmentProcess::setCompartmentDimension()
{
  setSubunitStart();
  if(Length && !Subunits)
    {
      Subunits = (unsigned)(Length/(DiffuseRadius*2));
    }
  if(Width && !Filaments)
    {
      Filaments = (unsigned)((Width-2*DiffuseRadius)/
                                 (DiffuseRadius*sqrt(3)))+1;
    }
  if(Periodic && Filaments%2 != 0)
    {
      ++Filaments;
    }
  //Need to use 2.00001 here to avoid rounding off error when calculating
  //LipidRows below:
  Width = 2.00001*DiffuseRadius+(Filaments-1)*DiffuseRadius*sqrt(3); 
  Height = 2*DiffuseRadius;
  if(!Filaments)
    {
      Filaments = 1;
      Width = 2.00001*DiffuseRadius;
    }
  if(!Subunits)
    {
      Subunits = 1;
      Length = 2.00001*DiffuseRadius;
    }
  if(Filaments == 1)
    {
      theDimension = 1;
      Length = Subunits*DiffuseRadius*2;
    }
  if(Subunits == 1)
    {
      theDimension = 1;
      Width = Filaments*DiffuseRadius*2;
    }
  if(Subunits != 1 && Filaments != 1)
    {
      theDimension = 2;
      //Add DiffuseRadius for the protrusion from hexagonal arrangement:
      Length = Subunits*DiffuseRadius*2+DiffuseRadius;
    }
  if(theLipidSpecies)
    {
      LipidCols = (unsigned)(Length/(LipidRadius*2));
      LipidRows = (unsigned)((Width-2*LipidRadius)/(LipidRadius*sqrt(3)))+1;
    }
  //TODO: make length, width, height consistent with the definitions
  allocateGrid();
}

void CompartmentProcess::allocateGrid()
{
  //in SpatiocyteStepper:
  //Normalized compartment lengths in terms of lattice voxel radius:
  nParentHeight = theComp->lengthX+nVoxelRadius*2;
  nParentWidth = theComp->lengthY+nVoxelRadius*2;
  nParentLength = theComp->lengthZ+nVoxelRadius*2;
  nLength = Length/(VoxelRadius*2);
  nWidth = Width/(VoxelRadius*2);
  nHeight = Height/(VoxelRadius*2);
  theComp->lengthX = nHeight;
  theComp->lengthY = nWidth;
  theComp->lengthZ = nLength;

  parentOrigin = theComp->centerPoint;
  Point s(parentOrigin);
  parentOrigin.x -= nParentHeight/2;
  parentOrigin.y -= nParentWidth/2;
  parentOrigin.z -= nParentLength/2;
  s = parentOrigin;
  gridCols = (unsigned)ceil(nParentLength/nGridSize);
  gridRows = (unsigned)ceil(nParentWidth/nGridSize);
  gridLayers = (unsigned)ceil(nParentHeight/nGridSize);
  theVacGrid.resize(gridCols*gridRows*gridLayers);
  /*
  if(theLipidSpecies)
    {
      theLipGrid.resize(gridCols*gridRows*gridLayers);
    }
    */
  //Actual surface area = Width*Length
  theComp->actualArea = Width*Length;
  theComp->specArea = Width*Length;
}

void CompartmentProcess::initializeCompartment()
{
  if(!isCompartmentalized)
    {
      thePoints.resize(endCoord-subStartCoord);
      vacStartIndex = theVacantSpecies->size();
      intStartIndex = theInterfaceSpecies->size();
      initializeVectors();
      initializeFilaments(subunitStart, Filaments, Subunits, nDiffuseRadius,
                          theVacantSpecies, subStartCoord);
      elongateFilaments(theVacantSpecies, subStartCoord, Filaments, Subunits,
                        nDiffuseRadius);
      connectFilaments(subStartCoord, Filaments, Subunits);
      setDiffuseSize(subStartCoord, lipStartCoord);
      setGrid(theVacantSpecies, theVacGrid, subStartCoord);
      interfaceSubunits();
      if(theLipidSpecies)
        {
          initializeFilaments(lipidStart, LipidRows, LipidCols, nLipidRadius,
                              theLipidSpecies, lipStartCoord);
          elongateFilaments(theLipidSpecies, lipStartCoord, LipidRows,
                            LipidCols, nLipidRadius);
          connectFilaments(lipStartCoord, LipidRows, LipidCols);
          setDiffuseSize(lipStartCoord, endCoord);
          //setGrid(theLipidSpecies, theLipGrid, lipStartCoord);
          setSpeciesIntersectLipids();
        }
      //To support walkRegular in SpatiocyteSpecies
      else if(RegularLattice)
        {
          setSpeciesIntersectVacants();
        }
      isCompartmentalized = true;
    }
  theVacantSpecies->setIsPopulated();
  if(theLipidSpecies)
    {
      theLipidSpecies->setIsPopulated();
    }
  theInterfaceSpecies->setIsPopulated();
  //theSpecies[2]->setIsPopulated();
}

void CompartmentProcess::setGrid(Species* aSpecies,
                                 std::vector<std::vector<unsigned> >& aGrid,
                                 unsigned aStartCoord)
{ 
  if(aSpecies)
    {
      for(unsigned i(vacStartIndex); i != aSpecies->size(); ++i)
        {
          Point& aPoint(*(*theLattice)[aStartCoord+i-vacStartIndex].point);
          const int row((int)((aPoint.y-parentOrigin.y)/nGridSize));
          const int col((int)((aPoint.z-parentOrigin.z)/nGridSize));
          const int layer((int)((aPoint.x-parentOrigin.x)/nGridSize));
          if(row >= 0 && row < gridRows && layer >= 0 && layer < gridLayers &&
             col >= 0 && col < gridCols)
            {
              aGrid[row+
                gridRows*layer+
                gridRows*gridLayers*col].push_back(aStartCoord+i-vacStartIndex);
            }
        }
    }
}

// y:width:rows:filaments
// z:length:cols:subunits
void CompartmentProcess::setSpeciesIntersectVacants()
{
  setAdjoinOffsets(Subunits);
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      theVacantCompSpecies[i]->setAdjoinOffsets(theAdjoinOffsets,
                                                theRowOffsets);
    }
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      theVacantCompSpecies[i]->setIntersectOffsets(theAdjoinOffsets,
                                                   theRowOffsets,
                                                   theVacantSpecies,
                                                   subunitStart, Filaments,
                                                   Subunits, nDiffuseRadius,
                                                   SubunitAngle,
                                                   surfaceNormal, widthVector,
                                                   lengthVector);
    }
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      theVacantCompSpecies[i]->setProductPairOffsets(theAdjoinOffsets,
                                                     theRowOffsets,
                                                     theVacantSpecies,
                                                     subunitStart, Filaments,
                                                     Subunits, nDiffuseRadius,
                                                     SubunitAngle,
                                                     surfaceNormal, widthVector,
                                                     lengthVector);
    }
}
void CompartmentProcess::setSpeciesIntersectLipids()
{
  setAdjoinOffsets(LipidCols);
  if(RegularLattice)
    {
      for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
        {
          theLipidCompSpecies[i]->setAdjoinOffsets(theAdjoinOffsets,
                                                   theRowOffsets);
        }
      for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
        {
          theVacantCompSpecies[i]->setIntersectOffsets(theAdjoinOffsets,
                                                       theRowOffsets,
                                                       theLipidSpecies,
                                                       lipidStart, Filaments,
                                                       Subunits, nLipidRadius,
                                                       SubunitAngle,
                                                       surfaceNormal,
                                                       widthVector,
                                                       lengthVector);
        }
      for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
        {
          theVacantCompSpecies[i]->setProductPairOffsets(theAdjoinOffsets,
                                                         theRowOffsets,
                                                         theLipidSpecies,
                                                         lipidStart, Filaments,
                                                         Subunits, nLipidRadius,
                                                         SubunitAngle,
                                                         surfaceNormal,
                                                         widthVector,
                                                         lengthVector);
        }
    }
  else
    {
      for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
        {
          if(theVacantCompSpecies[i]->getIsMultiscale())
            {
              theVacantCompSpecies[i]->setIntersectLipids(theLipidSpecies,
                                                lipidStart, nGridSize, gridCols,
                                                gridRows, theVacGrid, Filaments,
                                                Subunits);
            }
        }
    }
}

void CompartmentProcess::getStartVoxelPoint(Point& start, Point& nearest,
                                            Point& farthest)
{
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  Species* surface(aComp->surfaceSub->vacantSpecies);
  double dist(0);
  Point origin;
  origin.x = 0;
  origin.y = 0;
  origin.z = 0;
  if(surface->size())
    {
      nearest = theSpatiocyteStepper->coord2point(surface->getCoord(0));
      farthest = nearest; 
      origin.x = nearest.x;
      dist = distance(nearest, origin);
      start = nearest;
    }
  for(unsigned i(1); i < surface->size(); ++i)
    {
      Point aPoint(theSpatiocyteStepper->coord2point(surface->getCoord(i)));
      if(aPoint.x < nearest.x)
        {
          nearest.x = aPoint.x;
          origin.x = aPoint.x;
          dist = distance(aPoint, origin);
          start = nearest;
        }
      else if(aPoint.x == nearest.x)
        {
          origin.x = aPoint.x;
          double aDist(distance(aPoint, origin));
          if(aDist < dist)
            {
              dist = aDist;
              start = aPoint;
            }
        }
      if(aPoint.y < nearest.y)
        {
          nearest.y = aPoint.y;
        }
      if(aPoint.z < nearest.z)
        {
          nearest.z = aPoint.z;
        }
      if(aPoint.x > farthest.x)
        {
          farthest.x = aPoint.x;
        }
      if(aPoint.y > farthest.y)
        {
          farthest.y = aPoint.y;
        }
      if(aPoint.z > farthest.z)
        {
          farthest.z = aPoint.z;
        }
    }
}

void CompartmentProcess::initializeVectors()
{
  lengthStart = subunitStart;
  /*
  //For Lipid start:
  lengthStart.z -= nDiffuseRadius;
  lengthStart.y -= nDiffuseRadius;
  */
  
  Point& origin(theComp->centerPoint);
  Point tmp(sub(subunitStart, origin));
  rotate(tmp);
  subunitStart = add(tmp, origin);

  tmp = sub(lengthStart, origin);
  rotate(tmp);
  lengthStart = add(tmp, origin);

  rotate(lengthVector);
  lengthEnd = disp(lengthStart, lengthVector, nLength);

  rotate(widthVector);
  widthEnd = disp(lengthEnd, widthVector, nWidth);

  rotate(heightVector);
  heightEnd = disp(widthEnd, heightVector, nHeight);

  if(theLipidSpecies)
    {
      lipidStart = lengthStart;
      disp_(lipidStart, lengthVector, nLipidRadius);
      disp_(lipidStart, widthVector, nLipidRadius);
    }

  Point center(lengthStart);
  disp_(center, lengthVector, nLength/2);
  disp_(center, widthVector, nWidth/2);
  theComp->centerPoint = center;

  //Set up surface vectors:
  surfaceNormal = cross(lengthVector, widthVector);
  surfaceNormal = norm(surfaceNormal);
  surfaceDisplace = dot(surfaceNormal, widthEnd);
  lengthDisplace = dot(lengthVector, lengthStart);
  lengthDisplaceOpp = dot(lengthVector, lengthEnd);
  widthDisplace = dot(widthVector, lengthEnd);
  widthDisplaceOpp = dot(widthVector, widthEnd);
}

void CompartmentProcess::rotateAsParent(Point& V)
{
  theSpatiocyteStepper->rotateX(theComp->rotateX, &V, -1);
  theSpatiocyteStepper->rotateY(theComp->rotateY, &V, -1);
  theSpatiocyteStepper->rotateZ(theComp->rotateZ, &V, -1);
}

void CompartmentProcess::rotate(Point& V)
{
  rotateAsParent(V);
  theSpatiocyteStepper->rotateX(RotateX, &V, 1);
  theSpatiocyteStepper->rotateY(RotateY, &V, 1);
  theSpatiocyteStepper->rotateZ(RotateZ, &V, 1);
}

void CompartmentProcess::initializeFilaments(Point& aStartPoint, unsigned aRows,
                                             unsigned aCols, double aRadius,
                                             Species* aVacant,
                                             unsigned aStartCoord)
{
  //The first comp voxel must have the aStartCoord:
  if(aStartCoord != endCoord)
    {
      addCompVoxel(0, 0, aStartPoint, aVacant, aStartCoord, aCols);
    }
  for(unsigned i(1); i < aRows; ++i)
    {
      Point U(aStartPoint);
      disp_(U, widthVector, i*aRadius*sqrt(3)); 
      if(i%2 == 1)
        {
          disp_(U, lengthVector, -aRadius); 
        }
      addCompVoxel(i, 0, U, aVacant, aStartCoord, aCols);
    }
}

Voxel* CompartmentProcess::addCompVoxel(unsigned rowIndex, 
                                        unsigned colIndex,
                                        Point& aPoint,
                                        Species* aVacant,
                                        unsigned aStartCoord,
                                        unsigned aCols)
{
  const unsigned aCoord(aStartCoord+rowIndex*aCols+colIndex);
  Voxel& aVoxel((*theLattice)[aCoord]);
  aVoxel.point = &thePoints[aStartCoord-subStartCoord+rowIndex*aCols+colIndex];
  *aVoxel.point = aPoint;
  if(RegularLattice)
    {
      aVoxel.adjoiningSize = theDiffuseSize;
      aVoxel.diffuseSize = theDiffuseSize;
      aVoxel.adjoiningCoords = new unsigned int[theDiffuseSize];
      for(unsigned i(0); i != theDiffuseSize; ++i)
        {
          aVoxel.adjoiningCoords[i] = theNullCoord;
        }
    }
  else
    {
      aVoxel.adjoiningSize = 0;
    }
  aVacant->addCompVoxel(aCoord);
  return &aVoxel;
}

void CompartmentProcess::elongateFilaments(Species* aVacant,
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
          addCompVoxel(i, j, A, aVacant, aStartCoord, aCols);
        }
    }
}


void CompartmentProcess::connectSubunit(unsigned a, unsigned b, 
                                        unsigned adjoinA, unsigned adjoinB)
{
  Voxel& voxelA((*theLattice)[a]);
  Voxel& voxelB((*theLattice)[b]);
  if(RegularLattice)
    {
      voxelA.adjoiningCoords[adjoinA] = b;
      voxelB.adjoiningCoords[adjoinB] = a;
    }
  else
    {
      addAdjoin(voxelA, b);
      addAdjoin(voxelB, a);
    }
}

/*
 row0   row1    row2
 fil0   fil1    fil2
 [NW] [ NORTH ] [NE] sub0, col0
      [subunit]      sub1, col1
 [SW] [ SOUTH ] [SE] sub2, col2
*/

void CompartmentProcess::setAdjoinOffsets(const unsigned cols)
{
  theAdjoinOffsets.resize(2);
  theAdjoinOffsets[0].resize(theDiffuseSize);
  theAdjoinOffsets[0][NORTH] = -1;
  theAdjoinOffsets[0][SOUTH] = 1;
  theAdjoinOffsets[0][NW] = -cols;
  theAdjoinOffsets[0][SW] = -cols+1;
  theAdjoinOffsets[0][NE] = cols;
  theAdjoinOffsets[0][SE] = cols+1;

  theAdjoinOffsets[1].resize(theDiffuseSize);
  theAdjoinOffsets[1][NORTH] = -1;
  theAdjoinOffsets[1][SOUTH] = 1;
  theAdjoinOffsets[1][NW] = -cols-1;
  theAdjoinOffsets[1][SW] = -cols;
  theAdjoinOffsets[1][NE] = cols-1;
  theAdjoinOffsets[1][SE] = cols;

  theRowOffsets.resize(theDiffuseSize);
  theRowOffsets[NORTH] = 0;
  theRowOffsets[SOUTH] = 0;
  theRowOffsets[NW] = -1;
  theRowOffsets[SW] = -1;
  theRowOffsets[NE] = 1;
  theRowOffsets[SE] = 1;
}

// y:width:rows:filaments
// z:length:cols:subunits
void CompartmentProcess::connectFilaments(unsigned aStartCoord,
                                          unsigned aRows, unsigned aCols)
{
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
              if(j%2 == 1)
                {
                  if(j+1 < aRows)
                    {
                      //periodic NE-SW 
                      b = aStartCoord+(j+1)*aCols+aCols-1; 
                      connectSubunit(a, b, NE, SW);
                    }
                  else if(j == aRows-1)
                    {
                      //periodic NE-SW 
                      b = aStartCoord+aCols-1; 
                      connectSubunit(a, b, NE, SW);
                    }
                  //periodic NW-SE
                  b = aStartCoord+(j-1)*aCols+aCols-1; 
                  connectSubunit(a, b, NW, SE);
                }
            }
          if(j > 0)
            {
              if(j%2 == 1)
                {
                  //SW-NE
                  unsigned a(aStartCoord+j*aCols+i);
                  unsigned b(aStartCoord+(j-1)*aCols+i); 
                  connectSubunit(a, b, SW, NE);
                  if(i > 0)
                    {
                      //NW-SE
                      b = aStartCoord+(j-1)*aCols+(i-1); 
                      connectSubunit(a, b, NW, SE);
                    }
                }
              else
                {
                  //NW-SE
                  unsigned a(aStartCoord+j*aCols+i);
                  unsigned b(aStartCoord+(j-1)*aCols+i); 
                  connectSubunit(a, b, NW, SE);
                  if(i+1 < aCols)
                    {
                      //SW-NE
                      b = aStartCoord+(j-1)*aCols+(i+1); 
                      connectSubunit(a, b, SW, NE);
                    }
                }
            }
        }
      if(Periodic && aRows > 1)
        { 
          //periodic NW-SE
          unsigned a(aStartCoord+i); //row 0
          unsigned b(aStartCoord+(aRows-1)*aCols+i); 
          connectSubunit(a, b, NW, SE);
          if(i+1 < aCols)
            {
              //periodic SW-NE
              b = aStartCoord+(aRows-1)*aCols+(i+1); 
              connectSubunit(a, b, SW, NE);
            }
        }
    }
}

void CompartmentProcess::setDiffuseSize(unsigned start, unsigned end)
{
  for(unsigned i(start); i != end; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
      subunit.diffuseSize = subunit.adjoiningSize;
    }
}

void CompartmentProcess::addAdjoin(Voxel& aVoxel, unsigned coord)
{
  unsigned* temp(new unsigned[aVoxel.adjoiningSize+1]);
  for(unsigned int i(0); i != aVoxel.adjoiningSize; ++i)
    {
      //Avoid duplicated adjoins:
      if(aVoxel.adjoiningCoords[i] == coord)
        {
          delete[] temp;
          return;
        }
      temp[i] = aVoxel.adjoiningCoords[i];
    }
  delete[] aVoxel.adjoiningCoords;
  temp[aVoxel.adjoiningSize++] = coord;
  aVoxel.adjoiningCoords = temp;
}


Voxel* CompartmentProcess::getNearestVoxelToSurface(const unsigned subIndex,
                                                    double& nearestDist,
                                                    const bool isInterface)
{
  Voxel* nearestVoxel;
  Point subPoint(*(*theLattice)[subIndex+subStartCoord].point);
  Point a(subPoint);
  Point b(subPoint);
  disp_(a, lengthVector, 1.5*nVoxelRadius);
  disp_(a, widthVector, 1.5*nVoxelRadius);
  disp_(a, surfaceNormal, 1.5*nVoxelRadius);
  disp_(b, lengthVector, -1.5*nVoxelRadius);
  disp_(b, widthVector, -1.5*nVoxelRadius);
  disp_(b, surfaceNormal, -1.5*nVoxelRadius);
  Point top(Point(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)));
  Point bottom(Point(std::min(a.x, b.x), std::min(a.y, b.y),
                     std::min(a.z, b.z)));
  unsigned blRow(0);
  unsigned blLayer(0);
  unsigned blCol(0);
  theSpatiocyteStepper->point2global(bottom, blRow, blLayer, blCol);
  unsigned trRow(0);
  unsigned trLayer(0);
  unsigned trCol(0);
  theSpatiocyteStepper->point2global(top, trRow, trLayer, trCol);
  for(unsigned i(blRow); i <= trRow; ++i)
    {
      for(unsigned j(blLayer); j <= trLayer; ++j)
        {
          for(unsigned k(blCol); k <= trCol; ++k)
            {
              unsigned aCoord(theSpatiocyteStepper->global2coord(i, j, k)); 
              Voxel& aVoxel((*theLattice)[aCoord]);
              if((isInterface && getID(aVoxel) == 
                  theInterfaceSpecies->getID()) ||
                 (!isInterface && theSpecies[getID(aVoxel)]->getIsCompVacant()))
                { 
                  Point aPoint(theSpatiocyteStepper->coord2point(aCoord));
                  double aDist(std::fabs(getDisplacementToSurface(aPoint)));
                  if(aDist < nearestDist && isInside(aPoint))
                    {
                      nearestDist = aDist;
                      nearestVoxel = &aVoxel;
                    }
                }
            }
        }
    }
  return nearestVoxel;
}

Voxel* CompartmentProcess::getNearestVoxelToSubunit(const unsigned subIndex,
                                                    double& nearestDist,
                                                    const bool isInterface)
{
  Voxel* nearestVoxel;
  Point subPoint(*(*theLattice)[subIndex+subStartCoord].point);
  Point a(subPoint);
  Point b(subPoint);
  disp_(a, lengthVector, 1.5*nVoxelRadius);
  disp_(a, widthVector, 1.5*nVoxelRadius);
  disp_(a, surfaceNormal, 1.5*nVoxelRadius);
  disp_(b, lengthVector, -1.5*nVoxelRadius);
  disp_(b, widthVector, -1.5*nVoxelRadius);
  disp_(b, surfaceNormal, -1.5*nVoxelRadius);
  Point top(Point(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)));
  Point bottom(Point(std::min(a.x, b.x), std::min(a.y, b.y),
                     std::min(a.z, b.z)));
  unsigned blRow(0);
  unsigned blLayer(0);
  unsigned blCol(0);
  theSpatiocyteStepper->point2global(bottom, blRow, blLayer, blCol);
  unsigned trRow(0);
  unsigned trLayer(0);
  unsigned trCol(0);
  theSpatiocyteStepper->point2global(top, trRow, trLayer, trCol);
  for(unsigned i(blRow); i <= trRow; ++i)
    {
      for(unsigned j(blLayer); j <= trLayer; ++j)
        {
          for(unsigned k(blCol); k <= trCol; ++k)
            {
              unsigned aCoord(theSpatiocyteStepper->global2coord(i, j, k)); 
              Voxel& aVoxel((*theLattice)[aCoord]);
              if((isInterface && getID(aVoxel) == 
                  theInterfaceSpecies->getID()) ||
                 (!isInterface && theSpecies[getID(aVoxel)]->getIsCompVacant()))
                { 
                  Point aPoint(theSpatiocyteStepper->coord2point(aCoord));
                  double aDist(distance(subPoint, aPoint));
                  if(aDist < nearestDist && isInside(aPoint))
                    {
                      nearestDist = aDist;
                      nearestVoxel = &aVoxel;
                    }
                }
            }
        }
    }
  return nearestVoxel;
}


//get a voxel that intersects either a subunit or surface, with the
//shortest distance
void CompartmentProcess::addFirstInterface()
{
  const double delta(std::max(nDiffuseRadius, nVoxelRadius));
  double nearestDist(libecs::INF);
  Voxel* nearestVoxel(NULL);
  unsigned subIndex(0);
  while(nearestDist > nDiffuseRadius/2 && 
        subIndex < lipStartCoord-subStartCoord)
    {
      double tmpDist(libecs::INF); 
      Voxel* tmpVoxel(getNearestVoxelToSubunit(subIndex, tmpDist, false));
      if(tmpDist < nearestDist)
        {
          nearestDist = tmpDist;
          nearestVoxel = tmpVoxel;
        }
      tmpDist = libecs::INF; 
      tmpVoxel = getNearestVoxelToSurface(subIndex, tmpDist, false); 
      if(tmpDist < nearestDist)
        {
          nearestDist = tmpDist;
          nearestVoxel = tmpVoxel;
        }
      ++subIndex;
    }
  if(nearestVoxel != NULL && nearestDist <= 1.01*delta)
    {
      std::cout << "Found the first interface voxel, dist:" <<
       nearestDist << " subIndex:" << subIndex << " subMax:" <<
       lipStartCoord-subStartCoord << std::endl;
      addInterfaceVoxel(*nearestVoxel);
    }
  else
    {
      std::cout << "Couldn't find the first interface voxel, dist:" <<
       nearestDist << std::endl;
    }
}

void CompartmentProcess::addInterfaceVoxel(Voxel& aVoxel)
{
  //removeCompVoxel with getIndex(x) doesn't work when you combine multiple
  //volume compartments
  //Species* aCompVacant(theInterfaceSpecies->getVacantSpecies());
  //aCompVacant->removeCompVoxel(aCompVacant->getIndex(&aVoxel));
  theInterfaceSpecies->addMolecule(&aVoxel);
}


/* Summary
 Interface voxels
   adjoins[0 → diffuseSize-1] = original volume voxels
   adjoins[diffuseSize → adjoiningSize-1] = subunit voxels
 Subunit voxels
   adjoins[0 → diffuseSize-1] = original subunit voxels
   adjoins[diffuseSize → adjoiningSize-1] = volume voxels that are not 
                                            interface species
*/
void CompartmentProcess::connectSubunitInterfaceAdjoins()
{
  for(unsigned i(0); i != subunitInterfaces.size(); ++i)
    {
      for(unsigned j(0); j != subunitInterfaces[i].size(); ++j)
        {
          Voxel& subunit((*theLattice)[i+subStartCoord]);
          Voxel& interface(*theInterfaceSpecies->getMolecule(
                                        subunitInterfaces[i][j]));
          addAdjoin(interface, i+subStartCoord);
          for(unsigned k(0); k != interface.diffuseSize; ++k)
            {
              unsigned coord(interface.adjoiningCoords[k]);
              Voxel& adjoin((*theLattice)[coord]);
              if(theSpecies[getID(adjoin)]->getIsCompVacant() && 
                 getID(adjoin) != theInterfaceSpecies->getID() &&
                 isDissociationSide(adjoin.coord))
                {
                  addAdjoin(subunit, coord);
                }
            }
          //Remove adjoins of volume voxels pointing to interface if
          //they are not from the correction binding side of the surface:
          removeAdjoinsFromNonBindingSide(interface);
        }
    }
}

void CompartmentProcess::removeAdjoinsFromNonBindingSide(Voxel& interface)
{
  for(unsigned i(0); i != interface.diffuseSize; ++i)
    {
      unsigned coord(interface.adjoiningCoords[i]);
      Voxel& adjoin((*theLattice)[coord]);
      if(theSpecies[getID(adjoin)]->getIsCompVacant() && 
         getID(adjoin) != theInterfaceSpecies->getID() &&
         !isBindingSide(adjoin.coord))
        {
          for(unsigned j(0); j != adjoin.diffuseSize; ++j)
            {
              if(adjoin.adjoiningCoords[j] == interface.coord)
                {
                  //Point to itself it is pointing to interface:
                  adjoin.adjoiningCoords[j] = adjoin.coord;
                }
            }
        }
    }
}

bool CompartmentProcess::isBindingSide(const unsigned aCoord)
{
  Point aPoint(theSpatiocyteStepper->coord2point(aCoord));
  //Dissociating from subunit to volume voxels:
  if(BindingDirection == 2)
    {
      return true;
    }
  else if(BindingDirection == 0 && isOnAboveSurface(aPoint))
    {
      return true;
    }
  else if(BindingDirection == 1 && !isOnAboveSurface(aPoint))
    {
      return true;
    }
  return false;
}

bool CompartmentProcess::isDissociationSide(const unsigned aCoord)
{
  Point aPoint(theSpatiocyteStepper->coord2point(aCoord));
  //Dissociating from subunit to volume voxels:
  if(DissociationDirection == 2)
    {
      return true;
    }
  else if(DissociationDirection == 0 && isOnAboveSurface(aPoint))
    {
      return true;
    }
  else if(DissociationDirection == 1 && !isOnAboveSurface(aPoint))
    {
      return true;
    }
  return false;
}

void CompartmentProcess::setSubunitInterfaces()
{
  subunitInterfaces.resize(Filaments*Subunits);
  subunitInterfaceDists.resize(Filaments*Subunits);
  for(unsigned i(intStartIndex); i != theInterfaceSpecies->size(); ++i)
    {
      setNearestSubunit(i, 1);
    }
}

void CompartmentProcess::setNearestSubunitForOrphanInterfaces()
{
  std::vector<unsigned> interfaceSubunitPairs(theInterfaceSpecies->size(), 0);
  for(unsigned i(subStartCoord); i != lipStartCoord; ++i)
    {
      std::vector<unsigned>& interfaces(subunitInterfaces[i-subStartCoord]);
      for(unsigned j(0); j != interfaces.size(); ++j)
        {
          ++interfaceSubunitPairs[interfaces[j]];
        }
    }
  for(unsigned i(intStartIndex); i != interfaceSubunitPairs.size(); ++i)
    {
      if(!interfaceSubunitPairs[i])
        {
          setNearestSubunit(i, 5);
        }
    }
}

void CompartmentProcess::addSortedSubunitInterface(const unsigned subIndex,
                                                   const unsigned intIndex,
                                                   const double aDist,
                                       const unsigned maxSubunitInterfaceSize)
{
  std::vector<unsigned>& interfaces(subunitInterfaces[subIndex]);
  std::vector<double>& dists(subunitInterfaceDists[subIndex]);
  if(interfaces.size() < maxSubunitInterfaceSize)
    {
      interfaces.push_back(intIndex);
      dists.push_back(aDist);
      //Sort the list once it is full:
      if(maxSubunitInterfaceSize > 1 && 
         dists.size() == maxSubunitInterfaceSize)
        {
          std::vector<unsigned> tmpInterfaces(dists.size());
          std::vector<double> tmpDists(dists.size());
          for(unsigned i(0); i != dists.size(); ++i)
            {
              unsigned idx(dists.size()-1);
              const double dist(dists[i]);
              for(unsigned j(0); j < i; ++j)
                {
                  if(dist < dists[j])
                    {
                      --idx;
                    }
                }
              for(unsigned j(i+1); j < dists.size(); ++j)
                {
                  if(dist < dists[j])
                    {
                      --idx;
                    }
                  else if(dist == dists[j])
                    {
                      --idx;
                    }
                }
              tmpDists[idx] = dist;
              tmpInterfaces[idx] = interfaces[i];
            }
          subunitInterfaces[subIndex] = tmpInterfaces;
          subunitInterfaceDists[subIndex] = tmpDists;
        }
    }
  //Replace an existing interface if the current one is nearer to the subunit:
  else
    {
      for(int i(dists.size()-1); i >= 0; --i)
        {
          if(aDist < dists[i])
            {
              if(i < dists.size()-1)
                {
                  dists[i+1] = dists[i];
                  interfaces[i+1] = interfaces[i];
                }
              dists[i] = aDist;
              interfaces[i] = intIndex;
            }
          else
            {
              break;
            }
        }
    }
}

void CompartmentProcess::setNearestSubunit(const unsigned intIndex,
                                       const unsigned maxSubunitInterfaceSize)
{ 
  Voxel& aVoxel(*theInterfaceSpecies->getMolecule(intIndex));
  Point aPoint(theSpatiocyteStepper->coord2point(aVoxel.coord));
  const int row((int)((aPoint.y-parentOrigin.y)/nGridSize));
  const int col((int)((aPoint.z-parentOrigin.z)/nGridSize));
  const int layer((int)((aPoint.x-parentOrigin.x)/nGridSize)); 
  double nearestDist(libecs::INF);
  unsigned nearestCoord(0);
  if(row >= 0 && row < gridRows && layer >= 0 && layer < gridLayers &&
     col >= 0 && col < gridCols)
    {
      for(int i(std::max(0, layer-1)); i != std::min(layer+2, gridLayers); ++i)
        {
          for(int j(std::max(0, col-1)); j != std::min(col+2, gridCols); ++j)
            {
              for(int k(std::max(0, row-1)); k != std::min(row+2, gridRows);
                  ++k)
                {
                  const std::vector<unsigned>& coords(theVacGrid[k+gridRows*i+
                                                      gridRows*gridLayers*j]);
                  for(unsigned l(0); l < coords.size(); ++l)
                    {
                      const double dist(distance(
                                   *(*theLattice)[coords[l]].point, aPoint));
                      if(dist < nearestDist)
                        {
                          nearestCoord = coords[l];
                          nearestDist = dist;
                        }
                    }
                }
            }
        }
    }
  const double delta(nDiffuseRadius+nVoxelRadius);
  if(nearestDist != libecs::INF && nearestDist <= delta)
    {
      addSortedSubunitInterface(nearestCoord-subStartCoord, intIndex,
                                nearestDist, maxSubunitInterfaceSize);
    }
  else
    {
      if(Verbose)
        {
          std::cout << "i:" << intIndex << " nearestDist of subunit from " <<
                "orphan interface:" <<  nearestDist << " delta:" << delta <<
                std::endl;
        }
    }
}


void CompartmentProcess::setNearestInterfaceForOrphanSubunits()
{
  const double delta(std::max(nDiffuseRadius, nVoxelRadius)*2);
  for(unsigned i(subStartCoord); i != lipStartCoord; ++i)
    {
      unsigned subIndex(i-subStartCoord);
      if(!subunitInterfaces[subIndex].size())
        {
          double nearestDist(libecs::INF);
          Voxel* interface(getNearestVoxelToSubunit(subIndex, nearestDist,
                                                    true));
          if(nearestDist != libecs::INF && nearestDist <= delta)
            {
              addSortedSubunitInterface(subIndex, 
                                  theInterfaceSpecies->getIndex(interface),
                                  nearestDist, 1);
            }
          else
            {
              if(Verbose)
                {
                  std::cout << "i:" << i << " nearestDist of interface from " <<
                    "orphan subunit:" <<  nearestDist << " delta:" << delta <<
                    std::endl;
                }
            }
        }
    }
}


//To get better accuracy of actual volume size:
void CompartmentProcess::removeInterfaceCompVoxels()
{
  Species* aCompVacant(theInterfaceSpecies->getVacantSpecies());
  for(unsigned i(0); i != aCompVacant->size(); )
    { 
      Voxel* aVoxel(aCompVacant->getMolecule(i));
      if(getID(aVoxel) == theInterfaceSpecies->getID())
        {
          aCompVacant->removeCompVoxel(i);
        }
      else
        {
          ++i;
        }
    }
}


void CompartmentProcess::interfaceSubunits()
{
  addFirstInterface();
  extendInterfacesOverSurface();
  setSubunitInterfaces();


  /*
  std::vector<unsigned> sizes(12,0);
  for(unsigned i(subStartCoord); i != lipStartCoord; ++i)
    {
      sizes[subunitInterfaces[i-subStartCoord].size()]++;
    }
  for(unsigned i(0); i != sizes.size(); ++i)
    {
      std::cout << "before setNearestInterfaceForOrphanSubunits subunit size i:" << i << " " << sizes[i] <<
        std::endl;
    }
    */


  setNearestInterfaceForOrphanSubunits();


  /*
  for(unsigned i(0); i != sizes.size(); ++i)
    {
      sizes[i] = 0;
    }
  for(unsigned i(subStartCoord); i != lipStartCoord; ++i)
    {
      sizes[subunitInterfaces[i-subStartCoord].size()]++;
    }
  for(unsigned i(0); i != sizes.size(); ++i)
    {
      std::cout << "after size i:" << i << " " << sizes[i] << std::endl;
    }
  std::vector<unsigned> interfaceSubunitPairs(theInterfaceSpecies->size(), 0);
  for(unsigned i(subStartCoord); i != lipStartCoord; ++i)
    {
      std::vector<unsigned>& interfaces(subunitInterfaces[i-subStartCoord]);
      for(unsigned j(0); j != interfaces.size(); ++j)
        {
          ++interfaceSubunitPairs[interfaces[j]];
        }
    }
  std::vector<unsigned> pairs(20,0);
  for(unsigned i(intStartIndex); i != interfaceSubunitPairs.size(); ++i)
    {
      pairs[interfaceSubunitPairs[i]]++;
    }
  for(unsigned i(0); i != pairs.size(); ++i)
    {
      std::cout << "before setNearestSubunitForOrphanInterfaces interface size i:" << i << " " << pairs[i] <<
      std::endl;
    }
    */


  setNearestSubunitForOrphanInterfaces();


  /*
  for(unsigned i(0); i!= interfaceSubunitPairs.size(); ++i)
    {
      interfaceSubunitPairs[i] = 0;
    }
  for(unsigned i(subStartCoord); i != lipStartCoord; ++i)
    {
      std::vector<unsigned>& interfaces(subunitInterfaces[i-subStartCoord]);
      for(unsigned j(0); j != interfaces.size(); ++j)
        {
          ++interfaceSubunitPairs[interfaces[j]];
        }
    }
  for(unsigned i(0); i!= pairs.size(); ++i)
    {
      pairs[i] = 0;
    }
  for(unsigned i(intStartIndex); i != interfaceSubunitPairs.size(); ++i)
    {
      pairs[interfaceSubunitPairs[i]]++;
    }
  for(unsigned i(0); i != pairs.size(); ++i)
    {
      std::cout << "after interface size i:" << i << " " << pairs[i] <<
      std::endl;
    }
    */


  connectSubunitInterfaceAdjoins();
  //removeInterfaceCompVoxels();


  /*
  //check subunit adjoins
  std::vector<unsigned> adjoins(30,0);
  for(unsigned i(subStartCoord); i != lipStartCoord; ++i)
    {
      Voxel& subunit((*theLattice)[i]);
      adjoins[subunit.adjoiningSize-subunit.diffuseSize]++;
      if(subunit.adjoiningSize-subunit.diffuseSize<3)
        {
          //theSpecies[2]->addMolecule(&subunit);
        }
    }
  unsigned totalAdjs(0);
  unsigned totalAdjsE(0);
  for(unsigned i(0); i != adjoins.size(); ++i)
    {
      totalAdjsE += adjoins[i]*i;
      totalAdjs += adjoins[i];
      std::cout << "subunit adjoins size i:" << i << " " << adjoins[i] << 
      std::endl;
    }
  std::cout << "total adjs:" << totalAdjs << " extAdjs:" << totalAdjsE <<
    " ave:" << 
  double(totalAdjsE)/theVacantSpecies->size() << std::endl;

  //check interface adjoins
  std::vector<unsigned> adjoinsI(30,0);
  for(unsigned i(intStartIndex); i != theInterfaceSpecies->size(); ++i)
    {
      unsigned voxelCoord(theInterfaceSpecies->getCoord(i));
      Voxel& interface((*theLattice)[voxelCoord]);
      adjoinsI[interface.adjoiningSize-interface.diffuseSize]++;
      if(interface.adjoiningSize-interface.diffuseSize<1)
        {
          //theSpecies[2]->addMolecule(&interface);
        }
    }
  unsigned totalAdjsI(0);
  unsigned totalAdjsEI(0);
  for(unsigned i(0); i != adjoinsI.size(); ++i)
    {
      totalAdjsEI += adjoinsI[i]*i;
      totalAdjsI += adjoinsI[i];
      std::cout << "interface adjoins size i:" << i << " " << adjoinsI[i] <<
      std::endl;
    }
  std::cout << "total adjs:" << totalAdjsI << " extAdjs:" << totalAdjsEI << 
    " ave:" << double(totalAdjsEI)/theVacantSpecies->size() << std::endl;

  //Check interface duplicates:
  for(unsigned i(0); i != theInterfaceSpecies->size(); ++i)
    {
      for(unsigned j(0); j != theInterfaceSpecies->size(); ++j)
        {
          if(i != j && theInterfaceSpecies->getCoord(i) ==
             theInterfaceSpecies->getCoord(j))
            {
              std::cout << "duplicate i:" << i << " j:" << j << std::endl;
            }
        }
    }
    */
}


void CompartmentProcess::extendInterfacesOverSurface()
{
  for(unsigned i(intStartIndex); i != theInterfaceSpecies->size(); ++i)
    {
      unsigned voxelCoord(theInterfaceSpecies->getCoord(i));
      Voxel& anInterface((*theLattice)[voxelCoord]);
      for(unsigned j(0); j != theAdjoiningCoordSize; ++j)
        {
          Voxel& adjoin((*theLattice)[anInterface.adjoiningCoords[j]]);
          if(theSpecies[getID(adjoin)]->getIsCompVacant())
            {
              Point aPoint(theSpatiocyteStepper->coord2point(adjoin.coord));
              addSurfaceIntersectInterfaceVoxel(adjoin, aPoint);
            }
        }
    }
}


void CompartmentProcess::addSurfaceIntersectInterfaceVoxel(Voxel& aVoxel,
                                                           Point& aPoint)
{
  //We have already checked that aVoxel is a compVacant
  const double delta(std::max(nDiffuseRadius, nVoxelRadius)/2);
  double dispA(getDisplacementToSurface(aPoint));
  bool isInsideA(isInside(aPoint));
  if(dispA == 0 && getID(aVoxel) != theInterfaceSpecies->getID() &&
     isInsideA)
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
         getID(adjoin) != theInterfaceSpecies->getID() && isInsideB)
        {
          addInterfaceVoxel(adjoin);
        }
      if((dispA < 0) != (dispB < 0) &&
         getID(aVoxel) != theInterfaceSpecies->getID() &&
         getID(adjoin) != theInterfaceSpecies->getID())
        {
          //If aVoxel is nearer to the plane:
          if(std::fabs(dispA) <= std::fabs(dispB))
            {
              if(isInsideA)
                {
                  addInterfaceVoxel(aVoxel);
                }
              else if(theSpecies[getID(adjoin)]->getIsCompVacant() &&
                      isInsideB && std::fabs(dispB) <= delta)
                {
                  addInterfaceVoxel(adjoin);
                }
            }
          //If the adjoin is nearer to the plane:
          else if(theSpecies[getID(adjoin)]->getIsCompVacant())
            {
              if(isInsideB)
                {
                  addInterfaceVoxel(adjoin);
                }
              else if(isInsideA && std::fabs(dispA) <= delta)
                {
                  addInterfaceVoxel(aVoxel);
                }
            }
        }
    }
}

bool CompartmentProcess::isInside(Point& aPoint)
{
  const double dispA(point2planeDisp(aPoint, lengthVector, lengthDisplace));
  if(dispA >= -nDiffuseRadius && dispA <= nLength-2*nDiffuseRadius)
    {
      const double dispB(point2planeDisp(aPoint, widthVector, widthDisplace));
      if(dispB >= -nDiffuseRadius && dispB <= nWidth-nDiffuseRadius)
        {
          return true;
        }
    }
  return false;
}

double CompartmentProcess::getDisplacementToSurface(Point& aPoint)
{
  return point2planeDisp(aPoint, surfaceNormal, surfaceDisplace);
}


bool CompartmentProcess::isOnAboveSurface(Point& aPoint)
{
  if(getDisplacementToSurface(aPoint) >= 0)
    {
      return true;
    }
  return false;
}


void CompartmentProcess::printParameters()
{
  cout << getPropertyInterface().getClassName() << "[" <<
    getFullID().asString() << "]" << std::endl;
  switch (theComp->dimension)
    { 
      case 1:
          cout << "  Line compartment:" << std::endl;
          cout << "   Actual length:"<< Length << std::endl;
          break;
      case 2:
      default:
          cout << "  Surface compartment:" << std::endl;
          cout << "   Actual length:"<< Length << std::endl;
          cout << "   Actual width:"<< Width << std::endl;
          cout << "   Filaments:" << Filaments << std::endl;
          cout << "   Subunits:" << Subunits << std::endl;
          cout << "   Actual area:"<< theComp->actualArea << std::endl;
          cout << "   Specified area:"<< theComp->specArea << std::endl;
          break;
    }
  cout << "   Interface species size:"<< theInterfaceSpecies->size()
    << std::endl;
  if(theLipidSpecies)
    {
      cout << "  Lipid species:" << getIDString(theLipidSpecies) << 
        " number:" << theLipidSpecies->size() << std::endl;
      for(unsigned i(0); i != theLipidCompSpecies.size(); ++i)
        {
          cout << "    " << getIDString(theLipidCompSpecies[i]) <<
            " number:" << theLipidCompSpecies[i]->size() << std::endl;
        }
    } 
  cout << "  Vacant species:" << getIDString(theVacantSpecies) << 
    " number:" << theVacantSpecies->size() << std::endl;
  for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
    {
      cout << "    " << getIDString(theVacantCompSpecies[i]) <<
        " number:" << theVacantCompSpecies[i]->size() << std::endl;
    }
}
