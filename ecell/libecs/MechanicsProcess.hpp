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
extern "C" int idfrm,idphi,idvis,idvfx,idgam,idpsi,idsfr,idbfr,iddrg,idtrc,idhyc,nq,isoq[][4],
iqos[],lqos[],kqos[],ns;
extern "C" double dscp[][12],hvec[][3][12];
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
  MechanicsProcess() {}
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
  virtual bool isOnBelowSideSurface(Point&,Point&,Point&);
  void assignQuad();
  void assignNeigh();
  void fitMechanotoSpatio();
  void getBLTR(int);
  void getSurfaceCoords();
  void calculateSurfaceNormal(Point&,Point&,Point&);
  void populateSurface();

private:
  String FileName;
  char *idot1;
  char *ipf;
  char *eul;
  int idot;
  int idebug;
  int isve;
  int id1;
  int id2;
  double id3;
  double fixsurfaceDisplace;
  double logInt;
  double delt;
  double area;
  double *cmdt;
  double realHvec[345][3][12];
  double minhvecX;
  double minhvecY;
  double minhvecZ;
  double maxhvecX;
  double maxhvecY;
  double maxhvecZ;
  double voxelRadius; 
  double normVoxelRadius; 
  double lengthX;
  double lengthY;
  double lengthZ;
  std::vector<unsigned> surfaceCoords;
  std::vector<unsigned>row;
  std::vector<unsigned>col;
  std::vector<unsigned>lay;
  std::vector<unsigned>corn;
  std::vector<Species*> theVacantCompSpecies;
  std::vector<Point> newNode;
  std::vector<std::vector<int> > quadIndex;
  std::vector<std::vector<int> > neigh;
  Comp* theComp;  
  Point bottomLeft;
  Point topRight;
  Point fixsurfaceNormal; 
  Point AB;
  Point AC;
  Point vectorCut;
  Point planeCut;
  Variable* theVacantVariable;
  Species* theVacantSpecies;

};



}


#endif /* __MechanicsProcess_hpp */
