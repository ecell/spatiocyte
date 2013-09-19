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


#ifndef __MechanicsProcess_hpp
#define __MechanicsProcess_hpp

#include <fstream>
#include <iostream>
#include <string.h>
#include <libecs/SpatiocyteProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>

namespace libecs
{

extern "C" void extractnumeric(char *ipf, int &idot);
extern "C" void openfile(char *ipf,int &idiv,int &idot);
extern "C" void timedump(double &delt,double &logInt);
extern "C" double dscp[][12];
extern "C" int idfrm,idphi,idvis,idvfx,idgam,idpsi,idsfr,idbfr,iddrg,idtrc,idhyc;
extern "C" void initsvec();
extern "C" void openunit12();
extern "C" void initarea(double &area);
extern "C" void clchm(double &area);
extern "C" void cmstpsiz(double *tstep, double &j);
extern "C" double tnext,time_,tstop,tstp;
extern "C" void dfdriver(const char *eul);
extern "C" void chksurmol();
extern "C" void wrfile(int &k,int &isve, char *ipf,int &idot, double &delt, double &logInt);

//extern "C" double hvec[2048][3][12];

LIBECS_DM_CLASS(MechanicsProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(MechanicsProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(String, FileName);
    }
  MechanicsProcess():
    FileName("fitz.000") {}
  virtual ~MechanicsProcess() {}
  SIMPLE_SET_GET_METHOD(String, FileName); 
  virtual void prepreinitialize()
    {
      SpatiocyteProcess::prepreinitialize();
      theVacantVariable = createVariable("Vacant");
    }
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      SpatiocyteProcess::initialize();
      theVacantSpecies = theSpatiocyteStepper->addSpecies(theVacantVariable);
      for(VariableReferenceVector::iterator
          i(theVariableReferenceVector.begin());
          i != theVariableReferenceVector.end(); ++i)
        {
          Species* aSpecies(theSpatiocyteStepper->variable2species(
                                   (*i).getVariable())); 
          if(!(*i).getCoefficient())
            {
              theVacantCompSpecies.push_back(aSpecies);
            }
        }
      isPriorityQueued = true;
    }	
  virtual void initializeFirst()
    {
      SpatiocyteProcess::initializeFirst();
      theComp = new Comp;
      theVacantSpecies->setIsCompVacant();
      theVacantSpecies->setIsOffLattice();
      theVacantSpecies->setComp(theComp);
      for(unsigned i(0); i != theVacantCompSpecies.size(); ++i)
        {
          theVacantCompSpecies[i]->setIsOffLattice();
          //setVacantSpecies must be declared here since it needs
          //to be overwritten by DiffusionProcess in initializeSecond:
          theVacantCompSpecies[i]->setVacantSpecies(theVacantSpecies);
          theVacantCompSpecies[i]->setComp(theComp);
        }
    }

  virtual void initializeThird();
  virtual void initializeFifth();
  virtual void fire();
  virtual bool isOnAboveSurface(Point&);
  void assignQuad();

protected:
  String FileName;
   char *idot1;
   char *ipf;
   int idot;
   int idebug;
   int isve;
   int id1;
   int id2;
   double id3;
   double surfaceDisplace;
   double logInt;
   double delt;
   double area;
   double *cmdt;
   char *eul;
   Comp* theComp;  
   Point bottomLeft;
   Point topRight;
   Point surfaceNormal; 
   Point AB;
   Point AC;
   std::vector<unsigned> surfaceCoords;
   Variable* theVacantVariable;
   Species* theVacantSpecies;
   std::vector<Species*> theVacantCompSpecies;
   //std::vector<unsigned> quadIndex;
};



}


#endif /* __MechanicsProcess_hpp */
