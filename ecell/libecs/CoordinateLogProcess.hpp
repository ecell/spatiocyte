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


#ifndef __CoordinateLogProcess_hpp
#define __CoordinateLogProcess_hpp

#include <fstream>
#include <libecs/IteratingLogProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>

namespace libecs
{

LIBECS_DM_CLASS(CoordinateLogProcess, IteratingLogProcess)
{ 
public:
  LIBECS_DM_OBJECT(CoordinateLogProcess, Process)
    {
      INHERIT_PROPERTIES(IteratingLogProcess);
      PROPERTYSLOT_SET_GET(Integer, MoleculeID);
      PROPERTYSLOT_SET_GET(Integer, SpeciesID);
    }
  SIMPLE_SET_GET_METHOD(Integer, MoleculeID);
  SIMPLE_SET_GET_METHOD(Integer, SpeciesID);
  CoordinateLogProcess():
    MoleculeID(0),
    SpeciesID(0),
    theMoleculeSize(0)
  {
    FileName = "CoordinateLog.csv";
  }
  virtual ~CoordinateLogProcess() {}
  virtual void initializeLastOnce()
    {
      for(unsigned i(0); i != theProcessSpecies.size(); ++i)
        {
          theMoleculeSize += theProcessSpecies[i]->size();
        }
      theLogFile.open(FileName.c_str(), std::ios::trunc);
      theStartCoord = theSpatiocyteStepper->getStartCoord();
      saveFileHeader(theLogFile);
    }
  virtual void fire()
    {
      if(theTime <= LogEnd)
        {
          logSpecies();
          theTime += theInterval;
        }
      else
        {
          theTime = libecs::INF;
          theLogFile.flush();
          theLogFile.close();
        }
      thePriorityQueue->moveTop();
    }
  void logSpecies()
    {
      for(unsigned i(0); i != theProcessSpecies.size(); ++i)
        {
          theLogFile << getStepper()->getCurrentTime();
          logMolecules(i);
          theLogFile << std::endl;
        }
    }
protected:
  void logMolecules(int anIndex)
    {
      Species* aSpecies(theProcessSpecies[anIndex]);
      //For GFP species, whose molecules are not always up to date:
      aSpecies->updateMolecules();
      for(unsigned i(0); i != aSpecies->size(); ++i)
        {
          if(SpeciesID)
            {
              //For GFP tagged molecule, we get its actual species:
              theLogFile << "," << aSpecies->getID(aSpecies->getMolecule(i));
            }
          if(MoleculeID)
            {
              theLogFile << "," << i;
            }
          Point aPoint(aSpecies->getPoint(i));
          theLogFile << "," << aPoint.x << "," << aPoint.y << "," <<
            aPoint.z;
        }
    }
private:
  unsigned MoleculeID;
  unsigned SpeciesID;
  unsigned theStartCoord;
  double theMoleculeSize;
};

}

#endif /* __CoordinateLogProcess_hpp */
