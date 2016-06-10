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


#ifndef __ReactionProcess_hpp
#define __ReactionProcess_hpp

#include <fstream>
#include <libecs/SpatiocyteProcess.hpp>

namespace libecs
{

LIBECS_DM_CLASS(ReactionProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(ReactionProcess, Process)
    {
      INHERIT_PROPERTIES(SpatiocyteProcess);
      PROPERTYSLOT_SET_GET(Real, k);
      PROPERTYSLOT_SET_GET(Real, p);
      PROPERTYSLOT_SET_GET(Real, LogStart);
      PROPERTYSLOT_SET_GET(Integer, LogEvent);
      PROPERTYSLOT_SET_GET(Integer, SearchVacant);
      PROPERTYSLOT_GET_NO_LOAD_SAVE(Integer, Order);
      PROPERTYSLOT_SET_GET(String, FileName);
    }
  ReactionProcess():
    coefficientA(0),
    coefficientB(0),
    coefficientC(0),
    coefficientD(0),
    coefficientE(0),
    coefficientF(0),
    coefficientG(0),
    coefficientH(0),
    LogStart(0),
    LogEvent(0),
    SearchVacant(-1),
    theOrder(0),
    theEventCnt(0),
    k(-1),
    p(-1),
    A(NULL),
    B(NULL),
    C(NULL),
    D(NULL),
    E(NULL),
    F(NULL),
    G(NULL),
    H(NULL),
    variableA(NULL),
    variableB(NULL),
    variableC(NULL),
    variableD(NULL), 
    variableE(NULL),
    variableF(NULL),
    variableG(NULL),
    variableH(NULL),
    moleculeA(NULL),
    moleculeB(NULL),
    moleculeC(NULL),
    moleculeD(NULL),
    moleculeE(NULL),
    moleculeF(NULL),
    moleculeG(NULL),
    moleculeH(NULL),
    moleculeP(NULL),
    moleculeS(NULL),
    FileName("LogEvent.csv") {}
  virtual ~ReactionProcess() {}
  SIMPLE_SET_GET_METHOD(Real, k);
  SIMPLE_SET_GET_METHOD(Real, p);
  SIMPLE_SET_GET_METHOD(Real, LogStart);
  SIMPLE_SET_GET_METHOD(Integer, LogEvent);
  SIMPLE_SET_GET_METHOD(Integer, SearchVacant);
  SIMPLE_SET_GET_METHOD(String, FileName);
  virtual void fire()
    {
      requeue(); //theTop in thePriorityQueue is still this process since
      //we have not interrupted other processes to update their queue. 
      //So it is valid to call requeue, which only requeues theTop process, 
      //assuming it to be this process.
      for(unsigned i(0); i != theInterruptedProcesses.size(); ++i)
        { 
          theSpatiocyteStepper->addInterruptedProcess(
                                                theInterruptedProcesses[i]);
        }
    }
  virtual void interruptProcessesPre()
    {
      for(unsigned i(0); i != theInterruptedProcessesPre.size(); ++i)
        {
          theInterruptedProcessesPre[i]->interruptedPre(this);
        }
    }
  virtual void interruptProcessesPost()
    {
      logEvent();
      for(unsigned i(0); i != theInterruptedProcessesPost.size(); ++i)
        {
          theInterruptedProcessesPost[i]->interruptedPost(this);
        }
    }
  GET_METHOD(Integer, Order)
    {
      return theOrder;
    }
  void logEvent();
  virtual void initialize()
    {
      if(isInitialized)
        {
          return;
        }
      SpatiocyteProcess::initialize();
      //SearchVacant of each reaction is set according to the master
      //SearchVacant property of the SpatiocyteStepper if it is not
      //specified at the reaction level:
      if(SearchVacant == -1)
        {
          SearchVacant = theSpatiocyteStepper->getSearchVacant();
        }
      declareUnidirectional();
      calculateOrder();
    }
  virtual void initializeSecond()
    {
      SpatiocyteProcess::initializeSecond();
      theInterruptedProcesses.resize(0);
      theInterruptedProcessesPre.resize(0);
      theInterruptedProcessesPost.resize(0);
    }
  virtual void initializeLastOnce()
    {
      /*
      std::cout << getIDString() << std::endl;
      std::cout << "\tinterruptedProcesses:" << std::endl;
      for(unsigned i(0); i != theInterruptedProcesses.size(); ++i)
        {
          std::cout << "\t\t" << theInterruptedProcesses[i]->getIDString() << std::endl;
        }
      std::cout << "\tinterruptedProcessesPre:" << std::endl;
      for(unsigned i(0); i != theInterruptedProcessesPre.size(); ++i)
        {
          std::cout << "\t\t" << theInterruptedProcessesPre[i]->getIDString() << std::endl;
        }
      std::cout << "\tinterruptedProcessesPost:" << std::endl;
      for(unsigned i(0); i != theInterruptedProcessesPost.size(); ++i)
        {
          std::cout << "\t\t" << theInterruptedProcessesPost[i]->getIDString() << std::endl;
        }
        */
    }
  //Only ReactionProcesses can interrupt other processes because only 
  //they can change the number of molecules. 
  //This method is called to set the list of processes which will be
  //interrupted by this process. To determine if this process
  //needs to interrupt another process X, this process will call the
  //isDependentOn method of X with this process as the argument:
  virtual void setInterruption(std::vector<Process*> const &aProcessList)
    {
      for(std::vector<Process*>::const_iterator i(aProcessList.begin());
          i != aProcessList.end(); ++i)
        {
          ReactionProcess*
            aReactionProcess(dynamic_cast<ReactionProcess*>(*i));
          SpatiocyteProcess*
            aSpatiocyteProcess(dynamic_cast<SpatiocyteProcess*>(*i));
          const Process* me(dynamic_cast<Process*>(this));
          if(this != aReactionProcess && aSpatiocyteProcess->isDependentOn(me))
            {
              theInterruptedProcesses.push_back(aSpatiocyteProcess);
            }
          //A ReactionProcess can interruptPre and interruptPost itself
          //so we don't exclude this ReactionProcess here:
          if(aSpatiocyteProcess->isDependentOnPre(me))
            {
              theInterruptedProcessesPre.push_back(aSpatiocyteProcess);
            }
          if(aSpatiocyteProcess->isDependentOnPost(me))
            {
              theInterruptedProcessesPost.push_back(aSpatiocyteProcess);
            }
        }
    }
  virtual Species* getA() const
    {
      return A;
    }
  virtual Species* getB() const
    {
      return B;
    }
  virtual Species* getC() const
    {
      return C;
    }
  virtual Species* getD() const
    {
      return D;
    }
  virtual Species* getE() const
    {
      return E;
    }
  virtual Voxel* getMoleculeA()
    {
      return moleculeA;
    }
  virtual Voxel* getMoleculeB()
    {
      return moleculeB;
    }
  virtual Voxel* getMoleculeC()
    {
      return moleculeC;
    }
  virtual Voxel* getMoleculeD()
    {
      return moleculeD;
    }
  virtual Voxel* getMoleculeE()
    {
      return moleculeE;
    }
  virtual Voxel* getMoleculeP()
    {
      return moleculeP;
    }
  virtual Voxel* getMoleculeS()
    {
      return moleculeS;
    }
  virtual void addSubstrateInterrupt(Species* aSpecies, Voxel* aMolecule) {}
  virtual void removeSubstrateInterrupt(Species* aSpecies, Voxel* aMolecule) {}
protected:
  virtual void calculateOrder();
protected:
  int coefficientA;
  int coefficientB;
  int coefficientC;
  int coefficientD;
  int coefficientE;
  int coefficientF;
  int coefficientG;
  int coefficientH;
  int SearchVacant;
  int theOrder;
  unsigned theEventCnt;
  unsigned LogEvent;
  double k;
  double p;
  double LogStart;
  //Species are for non HD species:
  Species* A;
  Species* B;
  Species* C;
  Species* D;
  Species* E;
  Species* F;
  Species* G;
  Species* H;
  //Variables are for HD species:
  Variable* variableA;
  Variable* variableB;
  Variable* variableC;
  Variable* variableD;
  Variable* variableE;
  Variable* variableF;
  Variable* variableG;
  Variable* variableH;
  Voxel* moleculeA;
  Voxel* moleculeB;
  Voxel* moleculeC;
  Voxel* moleculeD;
  Voxel* moleculeE;
  Voxel* moleculeF;
  Voxel* moleculeG;
  Voxel* moleculeH;
  Voxel* moleculeP;
  Voxel* moleculeS;
  String FileName;
  std::ofstream theLogFile;
  std::vector<unsigned> theAdjoinSubstratesA; //-10
  std::vector<unsigned> theAdjoinProductsA;   //10
  std::vector<unsigned> theAdjoinSubstratesB; //-20
  std::vector<unsigned> theAdjoinProductsB;   //20
  std::vector<SpatiocyteProcess*> theInterruptedProcesses;
  std::vector<SpatiocyteProcess*> theInterruptedProcessesPre;
  std::vector<SpatiocyteProcess*> theInterruptedProcessesPost;
};

}

#endif /* __ReactionProcess_hpp */
