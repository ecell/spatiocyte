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

#include <libecs/VisualizationLogProcess.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(VisualizationLogProcess, Process); 

void VisualizationLogProcess::initializeLog()
{
  unsigned int aLatticeType(theSpatiocyteStepper->getLatticeType());
  theLogFile.write((char*)(&aLatticeType), sizeof(aLatticeType));
  theLogFile.write((char*)(&theMeanCount), sizeof(theMeanCount));
  unsigned int aStartMol(0);
  theLogFile.write((char*)(&aStartMol), sizeof(aStartMol));
  unsigned int aRowSize(theSpatiocyteStepper->getRowSize());
  theLogFile.write((char*)(&aRowSize), sizeof(aRowSize));
  unsigned int aLayerSize(theSpatiocyteStepper->getLayerSize());
  theLogFile.write((char*)(&aLayerSize), sizeof(aLayerSize));
  unsigned int aColSize(theSpatiocyteStepper->getColSize());
  theLogFile.write((char*)(&aColSize), sizeof(aColSize));
  Point aCenterPoint(theSpatiocyteStepper->getCenterPoint());
  double aRealRowSize(aCenterPoint.z*2);
  theLogFile.write((char*)(&aRealRowSize), sizeof(aRealRowSize));
  double aRealLayerSize(aCenterPoint.y*2);
  theLogFile.write((char*)(&aRealLayerSize), sizeof(aRealLayerSize));
  double aRealColSize(aCenterPoint.x*2);
  theLogFile.write((char*)(&aRealColSize), sizeof(aRealColSize));
  unsigned int aLatticeSpSize(theLatticeSpecies.size());
  theLogFile.write((char*)(&aLatticeSpSize), sizeof(aLatticeSpSize));
  unsigned int aPolymerSize(thePolymerSpecies.size());
  theLogFile.write((char*)(&aPolymerSize), sizeof(aPolymerSize));
  unsigned int aReservedSize(0);
  theLogFile.write((char*)(&aReservedSize), sizeof(aReservedSize));
  unsigned int anOffLatticeSpSize(theOffLatticeSpecies.size());
  theLogFile.write((char*)(&anOffLatticeSpSize), sizeof(anOffLatticeSpSize));
  //theLogMarker is a constant throughout the simulation:
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
  double aVoxelRadius(theSpatiocyteStepper->getVoxelRadius());
  theLogFile.write((char*)(&aVoxelRadius), sizeof(aVoxelRadius));
  for(unsigned int i(0); i != theLatticeSpecies.size(); ++i)
    {
      unsigned int aStringSize(
       theLatticeSpecies[i]->getVariable()->getFullID().asString().size());
      theLogFile.write((char*)(&aStringSize), sizeof(aStringSize));
      theLogFile.write(
       theLatticeSpecies[i]->getVariable()->getFullID().asString().c_str(),
       aStringSize);
      double aRadius(theLatticeSpecies[i]->getMolRadius());
      theLogFile.write((char*)(&aRadius), sizeof(aRadius));
    }
  for(unsigned int i(0); i!=thePolymerSpecies.size(); ++i)
    {
      unsigned int aPolymerIndex(thePolymerIndex[i]);
      theLogFile.write((char*) (&aPolymerIndex), sizeof(aPolymerIndex));
      double aRadius(thePolymerSpecies[i]->getMolRadius());
      theLogFile.write((char*)(&aRadius), sizeof(aRadius));
    }
  for(unsigned int i(0); i != theOffLatticeSpecies.size(); ++i)
    {
      unsigned int aStringSize(
       theOffLatticeSpecies[i]->getVariable()->getFullID().asString().size());
      theLogFile.write((char*)(&aStringSize), sizeof(aStringSize));
      theLogFile.write(
       theOffLatticeSpecies[i]->getVariable()->getFullID().asString().c_str(),
       aStringSize);
      double aRadius(theOffLatticeSpecies[i]->getMolRadius());
      theLogFile.write((char*)(&aRadius), sizeof(aRadius));
    }
}

void VisualizationLogProcess::logMols(int anIndex)
{
  Species* aSpecies(theLatticeSpecies[anIndex]);
  //No need to log lipid or non diffusing vacant molecules since we have
  //already logged them once during initialization:
  if(aSpecies->getIsCompVacant())
    {
      return;
    }
  //The remaining vacant molecules must be diffusive, so we need to update
  //them before logging their position:
  //Also update the molecules of the tag species:
  aSpecies->updateMols();
  theLogFile.write((char*)(&anIndex), sizeof(anIndex));
  //The species molecule size:
  int aSize(aSpecies->size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(unsigned i(0); i != theIDs->size(); ++i)
    {
      for(unsigned j(0); j != aSpecies->size(i); ++j)
        {
          unsigned aCoord((*theInfo)[i][aSpecies->getMol(i, j)].coord);
          theLogFile.write((char*)(&aCoord), sizeof(aCoord));
        }
    }
}  

void VisualizationLogProcess::logOffLattice(int anIndex)
{
  Species* aSpecies(theOffLatticeSpecies[anIndex]);
  if(aSpecies->getIsVacant() && !aSpecies->getIsDiffusiveVacant() &&
     !aSpecies->getIsReactiveVacant())
    {
      return;
    }
  aSpecies->updateMols();
  theLogFile.write((char*)(&anIndex), sizeof(anIndex));
  //The species molecule size:
  int aSize(aSpecies->size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(unsigned i(0); i != theIDs->size(); ++i)
    {
      for(unsigned j(0); j != aSpecies->size(i); ++j)
        {
          Point aPoint(aSpecies->getPoint(i, j));
          theLogFile.write((char*)(&aPoint), sizeof(aPoint));
        }
    }
}  

void VisualizationLogProcess::logPolymers(int anIndex)
{
  Species* aSpecies(thePolymerSpecies[anIndex]);
  theLogFile.write((char*)(&anIndex), sizeof(anIndex));
  //The species molecule size:
  int aSize(aSpecies->size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(unsigned i(0); i != theIDs->size(); ++i)
    {
      for(unsigned j(0); j != aSpecies->size(i); ++j)
        {
          Point aPoint(aSpecies->getPoint(i, j));
          theLogFile.write((char*)(&aPoint), sizeof(aPoint));
        }
    }
}  

void VisualizationLogProcess::logSourceMols(int anIndex)
{
  /*
  Species* aSpecies(thePolymerSpecies[anIndex]);
  int aSourceIndex(theLatticeSpecies.size()+anIndex);
  theLogFile.write((char*)(&aSourceIndex), sizeof(aSourceIndex));
  const std::vector<unsigned int> aMols(aSpecies->getSourceMols());
  int aSize(aMols.size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(unsigned int i(0); i != aMols.size(); ++i)
    {
      unsigned int aMol(aMols[i]);
      theLogFile.write((char*)(&aMol), sizeof(aMol));
    }
    */
}  

void VisualizationLogProcess::logTargetMols(int anIndex)
{
  /*
  Species* aSpecies(thePolymerSpecies[anIndex]);
  int aTargetIndex(theLatticeSpecies.size()+thePolymerSpecies.size()+anIndex);
  theLogFile.write((char*)(&aTargetIndex), sizeof(aTargetIndex));
  const std::vector<unsigned int> aMols(aSpecies->getTargetMols());
  int aSize(aMols.size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(unsigned int i(0); i != aMols.size(); ++i)
    {
      unsigned int aMol(aMols[i]);
      theLogFile.write((char*)(&aMol), sizeof(aMol));
    }
    */
}  

void VisualizationLogProcess::logSharedMols(int anIndex)
{
  /*
  Species* aSpecies(thePolymerSpecies[anIndex]);
  int aSharedIndex(theLatticeSpecies.size()+thePolymerSpecies.size()*2+anIndex);
  theLogFile.write((char*)(&aSharedIndex), sizeof(aSharedIndex));
  const std::vector<unsigned int> aMols(aSpecies->getSharedMols());
  int aSize(aMols.size());
  theLogFile.write((char*)(&aSize), sizeof(aSize)); 
  for(unsigned int i(0); i != aMols.size(); ++i)
    {
      unsigned int aMol(aMols[i]);
      theLogFile.write((char*)(&aMol), sizeof(aMol));
    }
    */
}  

void VisualizationLogProcess::logSpecies()
{
  double aCurrentTime(theSpatiocyteStepper->getCurrentTime());
  theLogFile.write((char*)(&aCurrentTime), sizeof(aCurrentTime));
  for(unsigned int i(0); i != theLatticeSpecies.size(); ++i)
    {
      logMols(i);
    }
  for(unsigned int i(0); i != thePolymerSpecies.size(); ++i)
    {
      logSourceMols(i);
    }
  for(unsigned int i(0); i != thePolymerSpecies.size(); ++i)
    {
      logTargetMols(i);
    }
  for(unsigned int i(0); i != thePolymerSpecies.size(); ++i)
    {
      logSharedMols(i);
    }
  /*
  for(unsigned int i(0); i!=theReservedSpecies.size(); ++i)
    {
      logReservedMols(i);
    }
    */
  //theLogMarker is a constant throughout the simulation:
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
  for(unsigned int i(0); i != thePolymerSpecies.size(); ++i)
    {
      logPolymers(i);
    }
  for(unsigned int i(0); i != theOffLatticeSpecies.size(); ++i)
    {
      logOffLattice(i);
    }
  //theLogMarker is a constant throughout the simulation:
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
}

void VisualizationLogProcess::logCompVacant()
{
  /*
  double aCurrentTime(theSpatiocyteStepper->getCurrentTime());
  theLogFile.write((char*)(&aCurrentTime), sizeof(aCurrentTime));
  for(unsigned int i(0); i != theLatticeSpecies.size(); ++i)
    {
      if(theLatticeSpecies[i]->getIsCompVacant())
        {
          Species* aVacantSpecies(theLatticeSpecies[i]);
          //The species index in the process:
          theLogFile.write((char*)(&i), sizeof(i));
          //The species molecule size:
          unsigned int aSize(aVacantSpecies->size());
          theLogFile.write((char*)(&aSize), sizeof(aSize)); 
          for(unsigned int j(0); j != aSize; ++j)
            {
              unsigned int aMol(aVacantSpecies->getMol(j));
              theLogFile.write((char*)(&aMol), sizeof(aMol));
            }  
        }
    }
    */
  //theLogMarker is a constant throughout the simulation:
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
  for(unsigned int k(0); k != theOffLatticeSpecies.size(); ++k)
    {
      if(theOffLatticeSpecies[k]->getIsCompVacant())
        {
          Species* aSpecies(theOffLatticeSpecies[k]);
          theLogFile.write((char*)(&k), sizeof(k));
          //The species molecule size:
          int aSize(aSpecies->size());
          theLogFile.write((char*)(&aSize), sizeof(aSize)); 
          for(unsigned i(0); i != theIDs->size(); ++i)
            {
              for(unsigned j(0); j != aSpecies->size(i); ++j)
                {
                  Point aPoint(aSpecies->getPoint(i, j));
                  theLogFile.write((char*)(&aPoint), sizeof(aPoint));
                }
            }
        }
    }
  theLogFile.write((char*)(&theLogMarker), sizeof(theLogMarker));
}

}
