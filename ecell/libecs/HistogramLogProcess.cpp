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

#include <HistogramLogProcess.hpp>


void HistogramLogProcess::initialize()
{
  if(isInitialized)
    {
      return;
    }
  IteratingLogProcess::initialize();
  theProcessSpecies.resize(0);
  for(VariableReferenceVector::iterator
      i(theSortedVariableReferences.begin());
      i != theSortedVariableReferences.end(); ++i)
    {
      Variable* aVariable((*i).getVariable()); 
      if(aVariable->getName() != "HD")
        {
          Species* aSpecies(theSpatiocyteStepper->addSpecies(aVariable));
          if((*i).getCoefficient() >= 0)
            {
              theProcessSpecies.push_back(aSpecies);
            }
          else
            {
              if((*i).getCoefficient() == -1)
                {
                  theCenterSpecies = aSpecies;
                }
              else
                {
                  theMarkerSpecies = aSpecies;
                }
            }
        }
    }
}

void HistogramLogProcess::initializeLastOnce()
{
  theTotalIterations = Iterations;
  SuperOrigin = theSpatiocyteStepper->system2Comp(getSuperSystem()
                                                  )->centerPoint;
  VoxelDiameter = theSpatiocyteStepper->getVoxelRadius()*2;
  nLength = Length/VoxelDiameter;
  nWidth = Width/VoxelDiameter;
  nHeight = RadialHeight/VoxelDiameter;
  nRadius = Radius/VoxelDiameter;
  binInterval = nLength/(Bins+1);
  if(RadialHeight)
    {
      nLength = nHeight;
      nWidth = nHeight;
      binInterval = 2*M_PI/(Bins);
    }
  initializeVectors();
  if(Density)
    {
      setVacantSizes();
    }
  if(theCenterSpecies)
    {
      populateCenterSpecies();
    }
  if(theMarkerSpecies)
    {
      populateMarkerSpecies();
    }
}

void HistogramLogProcess::saveFileData(std::ofstream& aFile,
                                       const unsigned iterations)
{
  double aTime(LogStart);
  for(unsigned i(0); i != theLogValues.size(); ++i)
    {
      saveATimePoint(aFile, aTime, iterations, i);
      aTime += LogInterval;
    }
  aFile.close();
}

void HistogramLogProcess::saveATimePoint(std::ofstream& aFile,
                                          const double aTime,
                                          const unsigned anIteration,
                                          const unsigned aCnt)
{
  for(unsigned i(0); i != Bins; ++i)
    {
      aFile << aTime << "," << i;
      for(unsigned j(0); j != theLogValues[aCnt][i].size(); ++j)
        {
          if(Density)
            {
              aFile << "," << theLogValues[aCnt][i][j]/anIteration/
                theVacantSizes[j][i];
              /*
              aFile << "," << theLogValues[aCnt][i][j]/anIteration/
                ((aTime-LogStart)/LogInterval)/theVacantSizes[j][i];
                */
            }
          else
            {
              aFile << "," << theLogValues[aCnt][i][j]/anIteration;
            }
        }
      aFile << std::endl << std::flush;
    }
}

void HistogramLogProcess::saveBackup()
{
  if(SaveCounts > 0 && 
     Iterations%(int)rint(theTotalIterations/SaveCounts) == 0)
    {
      std::string aFileName(FileName.c_str());
      aFileName = aFileName + ".back";
      cout << "Saving backup data in: " << aFileName << std::endl;
      std::ofstream aFile;
      aFile.open(aFileName.c_str(), std::ios::trunc);
      saveFileHeader(aFile);
      double aTime(LogStart);
      int completedIterations(theTotalIterations-Iterations);
      for(unsigned i(0); i != theLogValues.size(); ++i)
        {
          saveATimePoint(aFile, aTime, completedIterations, i);
          aTime += LogInterval;
        }
      aFile.close();
    }
}

void HistogramLogProcess::saveFileHeader(std::ofstream& aFile)
{
  if(theTotalIterations != 1)
    {
      aFile << std::setprecision(15);
    }
  aFile << "Time," << std::scientific;
  if(RadialHeight)
    {
      aFile << binInterval*180/M_PI;
    }
  else
    {
      aFile << binInterval*VoxelDiameter;
    }
  for(unsigned i(0); i != theProcessSpecies.size(); ++i)
    {
      aFile << "," << getIDString(theProcessSpecies[i]);
    }
  aFile << std::endl << std::flush;
}

void HistogramLogProcess::logValues()
{
  if(Iterations == theTotalIterations)
    {
      initLogValues();
    }
  if(Collision)
    {
      logCollision();
    }
  else
    {
      logDensity();
    }
}

void HistogramLogProcess::initLogValues()
{
  theLogValues.resize(timePointCnt+1);
  theLogValues[timePointCnt].resize(Bins);
  for(unsigned i(0); i != Bins; ++i)
    {
      theLogValues[timePointCnt][i].resize(theProcessSpecies.size(), 0);
    }
}

void HistogramLogProcess::logCollision()
{
  for(unsigned i(0); i != theProcessSpecies.size(); ++i)
    {
      Species* aSpecies(theProcessSpecies[i]);
      for(unsigned j(0); j != aSpecies->size(); ++j)
        {
          unsigned bin;
          if(isInside(bin, aSpecies->getPoint(j)))
            {
              theLogValues[timePointCnt][bin][i] += 
                aSpecies->getCollisionCnt(j);
            }
        }
    }
}

void HistogramLogProcess::setVacantSizes()
{
  theVacantSizes.resize(theProcessSpecies.size());
  for(unsigned i(0); i != theProcessSpecies.size(); ++i)
    {
      theVacantSizes[i].resize(Bins);
      Species* aSpecies(theProcessSpecies[i]->getVacantSpecies());
      for(unsigned j(0); j != aSpecies->size(); ++j)
        {
          unsigned bin;
          if(isInside(bin, aSpecies->getPoint(j)))
            {
              theVacantSizes[i][bin] += 1;
            }
        }
      /*
      for(unsigned j(0); j != Bins; ++j)
        {
          std::cout << "species:" << getIDString(theProcessSpecies[i]) <<
            " bin:" << j << " size:" << theVacantSizes[i][j] << std::endl;
        }
        */
    }
}

void HistogramLogProcess::populateMarkerSpecies()
{
  if(theMarkerSpecies)
    {
      theMarkerSpecies->clearMolecules();
      for(unsigned i(0); i != theLattice->size(); ++i)
        {
          Point aPoint(theMarkerSpecies->coord2point(i));
          unsigned bin;
          if(isInside(bin, aPoint))
            {
              theMarkerSpecies->softAddMolecule(&(*theLattice)[i]);
            }
        }
      theMarkerSpecies->setIsPopulated();
    }
}

void HistogramLogProcess::populateCenterSpecies()
{
  if(theCenterSpecies)
    {
      theCenterSpecies->clearMolecules();
      for(unsigned i(0); i != theLattice->size(); ++i)
        {
          Point aPoint(theCenterSpecies->coord2point(i));
          if(distance(CompOrigin, aPoint) < nHeight/2)
            {
              theCenterSpecies->softAddMolecule(&(*theLattice)[i]);
            }
        }
      theCenterSpecies->setIsPopulated();
    }
}

void HistogramLogProcess::logDensity()
{
  for(unsigned i(0); i != theProcessSpecies.size(); ++i)
    {
      Species* aSpecies(theProcessSpecies[i]);
      for(unsigned j(0); j != aSpecies->size(); ++j)
        {
          unsigned bin;
          if(isInside(bin, aSpecies->getPoint(j)))
            {
              theLogValues[timePointCnt][bin][i] += 1;
            }
        }
    }
}

void HistogramLogProcess::initializeVectors()
{
  //Minus end
  MinusLength.x = -nLength/2;
  MinusLength.y = 0;
  MinusLength.z = 0;
  MinusHeight.x = 0;
  MinusHeight.y = -nHeight/2;
  MinusHeight.z = 0;
  MinusWidth.x = 0;
  MinusWidth.y = 0;
  MinusWidth.z = -nWidth/2;
  Comp* aComp(theSpatiocyteStepper->system2Comp(getSuperSystem()));
  CompOrigin.x = OriginX*aComp->lengthX/2;
  CompOrigin.y = OriginY*aComp->lengthY/2;
  CompOrigin.z = OriginZ*aComp->lengthZ/2;
  //Rotated Minus end
  theSpatiocyteStepper->rotateX(aComp->rotateX, &MinusLength, -1);
  theSpatiocyteStepper->rotateY(aComp->rotateY, &MinusLength, -1);
  theSpatiocyteStepper->rotateZ(aComp->rotateZ, &MinusLength, -1);
  theSpatiocyteStepper->rotateX(RotateX, &MinusLength, 1);
  theSpatiocyteStepper->rotateY(RotateY, &MinusLength, 1);
  theSpatiocyteStepper->rotateZ(RotateZ, &MinusLength, 1);
  theSpatiocyteStepper->rotateX(aComp->rotateX, &MinusHeight, -1);
  theSpatiocyteStepper->rotateY(aComp->rotateY, &MinusHeight, -1);
  theSpatiocyteStepper->rotateZ(aComp->rotateZ, &MinusHeight, -1);
  theSpatiocyteStepper->rotateX(RotateX, &MinusHeight, 1);
  theSpatiocyteStepper->rotateY(RotateY, &MinusHeight, 1);
  theSpatiocyteStepper->rotateZ(RotateZ, &MinusHeight, 1);
  theSpatiocyteStepper->rotateX(aComp->rotateX, &MinusWidth, -1);
  theSpatiocyteStepper->rotateY(aComp->rotateY, &MinusWidth, -1);
  theSpatiocyteStepper->rotateZ(aComp->rotateZ, &MinusWidth, -1);
  theSpatiocyteStepper->rotateX(RotateX, &MinusWidth, 1);
  theSpatiocyteStepper->rotateY(RotateY, &MinusWidth, 1);
  theSpatiocyteStepper->rotateZ(RotateZ, &MinusWidth, 1);
  theSpatiocyteStepper->rotateX(aComp->rotateX, &CompOrigin, -1);
  theSpatiocyteStepper->rotateY(aComp->rotateY, &CompOrigin, -1);
  theSpatiocyteStepper->rotateZ(aComp->rotateZ, &CompOrigin, -1);

  add_(CompOrigin, SuperOrigin);
  add_(MinusLength, CompOrigin);
  add_(MinusHeight, CompOrigin);
  add_(MinusWidth, CompOrigin);
  lengthVector = sub(CompOrigin, MinusLength);
  heightVector = sub(CompOrigin, MinusHeight);
  widthVector = sub(CompOrigin, MinusWidth);
  norm_(lengthVector);
  norm_(heightVector);
  norm_(widthVector);
}

bool HistogramLogProcess::isInside(unsigned& bin, Point N)
{
  if(RadialHeight)
    {
      return isInsideRadial(bin, N);
    }
  else
    {
      return isInsideLength(bin, N);
    }
}

bool HistogramLogProcess::isInsideRadial(unsigned& bin, Point N)
{
  Point D(heightVector);
  Point E(MinusHeight);
  double t((D.x*N.x+D.y*N.y+D.z*N.z-D.x*E.x-D.y*E.y-D.z*E.z)/
           (D.x*D.x +D.y*D.y+D.z*D.z));
  if(t < 0 || t > nHeight)
    {
      return false;
    }
  double dist(distance(CompOrigin, N));
  if(InnerRadius > 0 && dist < InnerRadius/VoxelDiameter)
    {
      return false;
    }
  if(OuterRadius != libecs::INF && dist > OuterRadius/VoxelDiameter)
    {
      return false;
    }
  D = widthVector;
  E = CompOrigin;
  t = (D.x*N.x+D.y*N.y+D.z*N.z-D.x*E.x-D.y*E.y-D.z*E.z)/
    (D.x*D.x +D.y*D.y+D.z*D.z);
  Point radialVector(sub(N, CompOrigin));
  norm_(radialVector);
  double angle(acos(dot(radialVector, lengthVector)));
  if(t < 0)
    {
      angle = 2*M_PI-angle;
    }
  bin = (unsigned)floor(angle/binInterval);
  if(bin >= StartBin)
    {
      bin = bin-StartBin;
    }
  else
    {
      bin = Bins+bin-StartBin;
    }
  return true;
}

bool HistogramLogProcess::isInsideLength(unsigned& bin, Point N)
{
  Point D(lengthVector);
  Point E(MinusLength);
  double t((D.x*N.x+D.y*N.y+D.z*N.z-D.x*E.x-D.y*E.y-D.z*E.z)/
           (D.x*D.x +D.y*D.y+D.z*D.z));
  if(t < 0)
    {
      return false;
    }
  bin = (unsigned)floor(t/binInterval);
  if(bin >= Bins)
    {
      return false;
    }
  double dist(sqrt(pow(-N.x+E.x+D.x*t, 2)+pow(-N.y+E.y+D.y*t, 2)+
                   pow(-N.z+E.z+D.z*t, 2)));
  if(dist > nRadius)
    {
      return false;
    }
  return true;
}

LIBECS_DM_INIT_STATIC(HistogramLogProcess, Process); 
