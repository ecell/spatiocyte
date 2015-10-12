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

#include <algorithm>
#include <gsl/gsl_randist.h>
#include <boost/lexical_cast.hpp>
#include <libecs/Model.hpp>
#include <libecs/System.hpp>
#include <libecs/Stepper.hpp>
#include <libecs/Process.hpp>
#include <libecs/VariableReference.hpp>
#include <libecs/MoleculePopulateProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>
#include <libecs/SpatiocyteProcess.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(MoleculePopulateProcess, Process); 


void MoleculePopulateProcess::initialize()
{
  if(isInitialized)
    {
      return;
    }
  SpatiocyteProcess::initialize();
  isPriorityQueued = true;
  for(VariableReferenceVector::const_iterator
	i(theVariableReferenceVector.begin());
        i != theVariableReferenceVector.end(); ++i)
    {
      Variable* aVariable((*i).getVariable());
      if(aVariable->getName() == "HD")
        {
          THROW_EXCEPTION(ValueError, getPropertyInterface().getClassName() +
            " [" + getFullID().asString() + "]: " +  
            aVariable->getFullID().asString() + " is a HD species and " +
            "therefore cannot be populated");
        }
    }
}

void MoleculePopulateProcess::initializeSecond()
{
  SpatiocyteProcess::initializeSecond();
  for(std::vector<Species*>::const_iterator i(theProcessSpecies.begin());
      i != theProcessSpecies.end(); ++i)
    {
      (*i)->setPopulateProcess(this, GaussianSigma);
    }
}

void MoleculePopulateProcess::initializeFourth()
{
  SpatiocyteProcess::initializeFourth();
  std::vector<Variable*> aVariables;
  std::vector<Process*> aProcesses(theSpatiocyteStepper->getProcessVector());
  for(std::vector<Process*>::const_iterator j(aProcesses.begin());
      j!= aProcesses.end(); ++j)
    {
      MoleculePopulateProcess* aProcess(
                              dynamic_cast<MoleculePopulateProcess*>(*j));
      if(aProcess)
        {
          VariableReferenceVector aVariableReferences(
                                  (*j)->getVariableReferenceVector());
          for(VariableReferenceVector::const_iterator 
              k(aVariableReferences.begin());
              k != aVariableReferences.end(); ++k)
            {
              aVariables.push_back((*k).getVariable());
            }	
        }     
    }
  for(std::vector<Species*>::iterator i(theSpecies.begin());
      i !=theSpecies.end(); ++i)
    {
      if(!(*i)->getIsVacant() && !(*i)->getIsInterface() &&
         !(*i)->getIsPopulated() && (*i)->getVariable() &&
         (*i)->getVariable()->getValue() != 0 &&
         std::find(aVariables.begin(), aVariables.end(),
                   (*i)->getVariable()) == aVariables.end())
        {
          THROW_EXCEPTION(ValueError, getPropertyInterface().getClassName() +
            " [" + (*i)->getVariable()->getFullID().asString() + "]: " +  
            "has a non-zero value but is not populated.");
        }
    }
}

void MoleculePopulateProcess::fire()
{
  for(std::vector<Species*>::const_iterator i(theProcessSpecies.begin());
      i != theProcessSpecies.end(); ++i)
    {
      (*i)->removeMolecules();
      populateUniformSparse(*i);
    }
  theInterval = ResetTime;
  theTime += theInterval; 
  thePriorityQueue->move(theQueueID);
}

void MoleculePopulateProcess::populateGaussian(Species* aSpecies)
{
}

void MoleculePopulateProcess::populateUniformOnMultiscale(Species* aSpecies)
{
  Species* aVacantSpecies(aSpecies->getVacantSpecies());
  aVacantSpecies->updateMolecules();
  const unsigned aSize(aSpecies->getPopulateCoordSize());
  const unsigned aVacantSize(aVacantSpecies->getPopulatableSize());
  cout << "    Populating uniformly on multiscale vacant:" <<
    getIDString(aVacantSpecies) << " available size:" << aVacantSize <<
    std::endl;
  cout << "       " << getIDString(aSpecies) << " current size:" <<
    aSpecies->size() << ", populate size:" << aSize << std::endl;
  if(!aSpecies->getIsPopulated())
    {
      if(UniformLengthX == 1 && UniformLengthY == 1 && UniformLengthZ == 1 &&
         !UniformRadiusXY && !UniformRadiusXZ && !UniformRadiusYZ &&
         !OriginX && !OriginY && !OriginZ && !EdgeX)
        {
          if(aVacantSize < aSize)
            {
              THROW_EXCEPTION(ValueError, String(
                              getPropertyInterface().getClassName()) +
                              "[" + getFullID().asString() + "]: There are " +
                              int2str(aSize) + " " + getIDString(aSpecies) +
                              " molecules that must be uniformly populated," +
                              "\nbut there are only " +
                              int2str(aVacantSize) + 
                              " multiscale vacant voxels of " +
                              getIDString(aSpecies->getVacantSpecies()) +
                              " that can be populated on.");
            }
          for(unsigned i(0); i != aSize; ++i)
            {
              Voxel* aMolecule;
              //After a molecule is added, the diffusive vacant species
              //molecule list is not updated, so we need to check if
              //it really is still vacant:
              aMolecule = aVacantSpecies->getRandomValidMolecule();
              aSpecies->addMolecule(aMolecule);
            }
        }
      else
        {
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
}

void MoleculePopulateProcess::populateUniformOnDiffusiveVacant(Species*
                                                               aSpecies)
{
  Species* aVacantSpecies(aSpecies->getVacantSpecies());
  aVacantSpecies->updateMolecules();
  const unsigned aSize(aSpecies->getPopulateCoordSize());
  const unsigned aVacantSize(aVacantSpecies->getPopulatableSize());
  cout << "    Populating uniformly on diffusive vacant:" <<
    getIDString(aVacantSpecies) << " available size:" << aVacantSize <<
    std::endl;
  cout << "       " << getIDString(aSpecies) << " current size:" <<
    aSpecies->size() << ", populate size:" << aSize << std::endl;
  if(!aSpecies->getIsPopulated())
    {
      if(UniformLengthX == 1 && UniformLengthY == 1 && UniformLengthZ == 1 &&
         !UniformRadiusXY && !UniformRadiusXZ && !UniformRadiusYZ &&
         !OriginX && !OriginY && !OriginZ && !EdgeX)
        {
          if(aVacantSize < aSize)
            {
              THROW_EXCEPTION(ValueError, String(
                              getPropertyInterface().getClassName()) +
                              "[" + getFullID().asString() + "]: There are " +
                              int2str(aSize) + " " + getIDString(aSpecies) +
                              " molecules that must be uniformly populated," +
                              "\nbut there are only " +
                              int2str(aVacantSize) + 
                              " diffuse vacant voxels of " +
                              getIDString(aSpecies->getVacantSpecies()) +
                              " that can be populated on.");
            }
          for(unsigned i(0); i != aSize; ++i)
            {
              Voxel* aMolecule;
              //After a molecule is added, the diffusive vacant species
              //molecule list is not updated, so we need to check if
              //it really is still vacant:
              aMolecule = aVacantSpecies->getRandomValidMolecule();
              aSpecies->addMolecule(aMolecule);
            }
        }
      else
        {
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
}

void MoleculePopulateProcess::populateUniformDense(Species* aSpecies,
                                              unsigned* aList, 
                                              unsigned* aCount)
{
  cout << "      Populating densely:" <<
    getIDString(aSpecies) << " current size:" << aSpecies->size() <<
    ", populate size:" << aSpecies->getPopulateCoordSize() << std::endl;
  Species* aVacantSpecies(aSpecies->getVacantSpecies());
  if(!aSpecies->getIsPopulated())
    {
      if(theLengthBinFractions.size())
        {
          populateBinFractions(aSpecies);
        }
      else if(UniformLengthX == 1 && UniformLengthY == 1 &&
              UniformLengthZ == 1 &&
              !UniformRadiusXY && !UniformRadiusXZ && !UniformRadiusYZ &&
              !OriginX && !OriginY && !OriginZ && !EdgeX)
        {
          unsigned aSize(aSpecies->getPopulateCoordSize());
          for(unsigned j(0); j != aSize; ++j)
            {
              unsigned aCoord;
              do
                {
                  aCoord = aVacantSpecies->getCoord(aList[(*aCount)++]); 
                }
              while(!aSpecies->isPopulatable(&(*theLattice)[aCoord]));
              aSpecies->addMolecule(&(*theLattice)[aCoord]);
            }
        }
      else
        {
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
}

void MoleculePopulateProcess::populateUniformSparse(Species* aSpecies)
{
  cout << "      Populating sparsely:" <<
    getIDString(aSpecies) << " current size:" << aSpecies->size() <<
    ", populate size:" << aSpecies->getPopulateCoordSize() << std::endl;
  Species* aVacantSpecies(aSpecies->getVacantSpecies());
  if(!aSpecies->getIsPopulated())
    {
      if(theLengthBinFractions.size())
        {
          populateBinFractions(aSpecies);
        }
      else if(UniformLengthX == 1 && UniformLengthY == 1 && 
              UniformLengthZ == 1 &&
              UniformLength == 1 && UniformWidth == 1 && UniformHeight == 1 &&
              !UniformRadiusXY && !UniformRadiusXZ && !UniformRadiusYZ &&
              !OriginX && !OriginY && !OriginZ && !EdgeX)
        {
          unsigned aSize(aSpecies->getPopulateCoordSize());
          int availableVoxelSize(aVacantSpecies->size());
          for(unsigned j(0); j != aSize; ++j)
            {
              unsigned aCoord;
              do
                {
                  aCoord = aVacantSpecies->getCoord(gsl_rng_uniform_int(
                                getStepper()->getRng(), availableVoxelSize));
                }
              while(!aSpecies->isPopulatable(&(*theLattice)[aCoord]));
              aSpecies->addMolecule(&(*theLattice)[aCoord]);
            }
        }
      else
        { 
          populateUniformRanged(aSpecies);
        }
      aSpecies->setIsPopulated();
    }
  aSpecies->updateMolecules();
}

void MoleculePopulateProcess::populateBinFractions(Species* aSpecies)
{

  cout << "        Populating bin fractions:" <<
    getIDString(aSpecies) << " current size:" << aSpecies->size() <<
    ", populate size:" << aSpecies->getPopulateCoordSize() << std::endl;
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  Species* aVacantSpecies(aSpecies->getVacantSpecies());
  Point C(aComp->centerPoint);
  Point start(disp(C, aComp->lengthVector, -aComp->nLength/2));
  const double lengthDisplace(dot(aComp->lengthVector, start));
  const unsigned nBins(theLengthBinFractions.size());
  const double binLength(aComp->nLength/nBins);
  std::vector<std::vector<unsigned> > aCoords;
  aCoords.resize(nBins);
  for(unsigned i(0); i != aVacantSpecies->size(); ++i)
    {
      unsigned aCoord(aVacantSpecies->getCoord(i));
      Point aPoint(aVacantSpecies->coord2point(aCoord));
      const double disp(point2planeDisp(aPoint, aComp->lengthVector,
                                        lengthDisplace));
      unsigned bin(std::max(disp/binLength, 0.0));
      bin = std::min(bin, nBins-1);
      Voxel* aVoxel(&(*theLattice)[aCoord]);
      if(aSpecies->isPopulatable(aVoxel))
        {
          aCoords[bin].push_back(aCoord);
        }
    }
  for(unsigned i(0); i != nBins; ++i)
    {
      std::vector<unsigned>& binCoords(aCoords[i]);
      std::random_shuffle(binCoords.begin(), binCoords.end());
      const unsigned size(binCoords.size()*theLengthBinFractions[i]);
      if(size)
        {
          std::cout << "          bin:" << i << " fraction:" <<
            theLengthBinFractions[i] << " populated size:" << size << std::endl;
        }
      for(unsigned j(0); j != size; ++j)
        {
          Voxel* aVoxel(&(*theLattice)[binCoords[j]]);
          aSpecies->addMolecule(aVoxel);
        }
    }
  aSpecies->setInitCoordSize(aSpecies->size());
}

//Doesn't work yet for multiscale
void MoleculePopulateProcess::populateUniformRanged(Species* aSpecies)
{
  cout << "        Populating uniformly ranged:" <<
    getIDString(aSpecies) << " current size:" << aSpecies->size() <<
    ", populate size:" << aSpecies->getPopulateCoordSize();
  Comp* aComp(aSpecies->getComp());
  Species* aVacantSpecies(aSpecies->getVacantSpecies());
  double deltaX(0);
  double deltaY(0);
  double deltaZ(0);
  // Increase the compartment dimensions by delta if it is a surface 
  // compartment:
  if(aComp->dimension == 2)
    {
      deltaX = theSpatiocyteStepper->getNormalizedVoxelRadius()*6/
        aComp->lengthX;
      deltaY = theSpatiocyteStepper->getNormalizedVoxelRadius()*6/
        aComp->lengthY;
      deltaZ = theSpatiocyteStepper->getNormalizedVoxelRadius()*6/
        aComp->lengthZ;
    }
  std::vector<unsigned> aCoords;
  //Get coords at the edge most 
  if(EdgeX)
    {
      unsigned aSize(aSpecies->getPopulateCoordSize());
      std::vector<double> pos;
      for(unsigned i(0); i != aVacantSpecies->size(); ++i)
        {
          unsigned aCoord(aVacantSpecies->getCoord(i));
          if(getID((*theLattice)[aCoord]) == aSpecies->getVacantID())
            {
              Point aPoint(aVacantSpecies->coord2point(aCoord));
              if(aCoords.size() < aSize)
                {
                  aCoords.push_back(aCoord);
                  pos.push_back(aPoint.x);
                }
              else
                {
                  int index(-1);
                  double max(-1);
                  for(unsigned j(0); j != pos.size(); ++j)
                    { 
                      if(EdgeX > 0 && aPoint.x-pos[j] > max)
                        {
                          index = j;
                          max = aPoint.x-pos[j];
                        }
                      else if(EdgeX < 0 && pos[j]-aPoint.x > max)
                        {
                          index = j;
                          max = pos[j]-aPoint.x;
                        }
                    }
                  if(index != -1)
                    {
                      pos[index] = aPoint.x;
                      aCoords[index] = aCoord;
                    }
                }
            }
        }
    }
  else if(!UniformRadiusXY && !UniformRadiusXZ && !UniformRadiusYZ)
    {
      //This is for CompartmentProcess that has defined Comp->widthVector,
      //Comp->heightVector and Comp->lengthVector:
      if(UniformWidth != 1 || UniformLength != 1 || UniformHeight != 1)
        {
          double del(0.0125);
          double maxL(std::min(1.1, OriginX+UniformLength+del));
          double minL(std::max(-1.1, OriginX-UniformLength-del));
          double maxW(std::min(1.1, OriginY+UniformWidth+del));
          double minW(std::max(-1.1, OriginY-UniformWidth-del));
          double maxH(std::min(1.1, OriginZ+UniformHeight+del));
          double minH(std::max(-1.1, OriginZ-UniformHeight-del)); 
          Point C(aComp->centerPoint);
          Point start(disp(C, aComp->widthVector, minW*aComp->nWidth/2));
          start = disp(start, aComp->lengthVector, minL*aComp->nLength/2);
          start = disp(start, aComp->heightVector, minH*aComp->nHeight/2);
          Point end(disp(C, aComp->widthVector, maxW*aComp->nWidth/2));
          end = disp(end, aComp->lengthVector, maxL*aComp->nLength/2);
          end = disp(end, aComp->heightVector, maxH*aComp->nHeight/2);
          double maxX(std::max(start.x, end.x));
          double minX(std::min(start.x, end.x));
          double maxY(std::max(start.y, end.y));
          double minY(std::min(start.y, end.y));
          double maxZ(std::max(start.z, end.z));
          double minZ(std::min(start.z, end.z));
          for(unsigned i(0); i != aVacantSpecies->size(); ++i)
            {
              unsigned aCoord(aVacantSpecies->getCoord(i));
              Point aPoint(aVacantSpecies->coord2point(aCoord));
              if(getID((*theLattice)[aCoord]) == aSpecies->getVacantID() &&
                 aPoint.x <= maxX && aPoint.x >= minX &&
                 aPoint.y <= maxY && aPoint.y >= minY &&
                 aPoint.z <= maxZ && aPoint.z >= minZ)
                {
                  aCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          double maxX(std::min(1.0, OriginX+UniformLengthX));
          double minX(std::max(-1.0, OriginX-UniformLengthX));
          double maxY(std::min(1.0, OriginY+UniformLengthY));
          double minY(std::max(-1.0, OriginY-UniformLengthY));
          double maxZ(std::min(1.0, OriginZ+UniformLengthZ));
          double minZ(std::max(-1.0, OriginZ-UniformLengthZ)); 
          maxX = aComp->centerPoint.x + maxX*aComp->lengthX/2*(1+deltaX);
          minX = aComp->centerPoint.x + minX*aComp->lengthX/2*(1+deltaX);
          maxY = aComp->centerPoint.y + maxY*aComp->lengthY/2*(1+deltaY);
          minY = aComp->centerPoint.y + minY*aComp->lengthY/2*(1+deltaY);
          maxZ = aComp->centerPoint.z + maxZ*aComp->lengthZ/2*(1+deltaZ);
          minZ = aComp->centerPoint.z + minZ*aComp->lengthZ/2*(1+deltaZ);
          for(unsigned i(0); i != aVacantSpecies->size(); ++i)
            {
              unsigned aCoord(aVacantSpecies->getCoord(i));
              Point aPoint(aVacantSpecies->coord2point(aCoord));
              if(getID((*theLattice)[aCoord]) == aSpecies->getVacantID() &&
                 aPoint.x < maxX && aPoint.x > minX &&
                 aPoint.y < maxY && aPoint.y > minY &&
                 aPoint.z < maxZ && aPoint.z > minZ)
                {
                  aCoords.push_back(aCoord);
                }
            }
        }
    }
  else
    {
      double nRadius(0);
      if(UniformRadiusXY)
        {
          nRadius = UniformRadiusXY/(theSpatiocyteStepper->getVoxelRadius()*2);
        }
      else if(UniformRadiusXZ)
        {
          nRadius = UniformRadiusXZ/(theSpatiocyteStepper->getVoxelRadius()*2);
        }
      else
        {
          nRadius = UniformRadiusYZ/(theSpatiocyteStepper->getVoxelRadius()*2);
        }
      double nStart(0);
      if(UniformRadiusWidth > 0)
        {
          nStart = nRadius-UniformRadiusWidth/
            (theSpatiocyteStepper->getVoxelRadius()*2);
          nStart = std::max(0.0, nStart);
        }
      Point C(aComp->centerPoint);
      disp_(C, aComp->widthVector, OriginY*aComp->nWidth/2);
      disp_(C, aComp->lengthVector, OriginX*aComp->nLength/2);
      disp_(C, aComp->heightVector, OriginZ*aComp->nHeight/2);
      for(unsigned i(0); i != aVacantSpecies->size(); ++i)
        {
          const unsigned aCoord(aVacantSpecies->getCoord(i));
          Point aPoint(aVacantSpecies->coord2point(aCoord));
          if(UniformRadiusXY)
            {
              aPoint.z = C.z;
            }
          else if(UniformRadiusXZ)
            {
              aPoint.y = C.y;
            }
          else
            {
              aPoint.x = C.x;
            }
          const double aDist(distance(aPoint, C));
          if(getID((*theLattice)[aCoord]) == aSpecies->getVacantID() &&
             aDist >= nStart && aDist <= nRadius) 
            {
              aCoords.push_back(aCoord);
            }
        }
    }
  unsigned aSize(aSpecies->getPopulateCoordSize());
  cout << ", vacant size:" << aCoords.size() << std::endl;
  if(aCoords.size() < aSize)
    {
      THROW_EXCEPTION(ValueError, String(
                      getPropertyInterface().getClassName()) +
                      "[" + getFullID().asString() + "]: There are " +
                      int2str(aSize) + " " + getIDString(aSpecies) +
                      " molecules that must be uniformly populated in a " +
                      "given range,\n but there are only " +
                      int2str(aCoords.size()) + " vacant voxels of " +
                      getIDString(aSpecies->getVacantSpecies()) +
                      " that can be populated.");
    }
  std::random_shuffle(aCoords.begin(), aCoords.end());
  for(unsigned i(0); i != aCoords.size(); ++i)
    {
      Voxel* aVoxel(&(*theLattice)[aCoords[i]]);
      if(aSpecies->size() < aSize && aSpecies->isPopulatable(aVoxel))
        {
          aSpecies->addMolecule(aVoxel);
        }
    }
}

}

