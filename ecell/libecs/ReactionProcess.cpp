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

#include <libecs/ReactionProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(ReactionProcess, Process);

void ReactionProcess::calculateOrder()
{ 
  theOrder = 0;
  for(VariableReferenceVector::iterator 
      i(theSortedVariableReferences.begin());
      i != theSortedVariableReferences.end(); ++i)
    {
      const int aCoefficient((*i).getCoefficient());
      Variable* aVariable((*i).getVariable());
      if(aCoefficient < 0)
        {
          if(aCoefficient <= -10)
            {
              if(aVariable->getName() == "HD")
                {
                  THROW_EXCEPTION(ValueError,
                     String(getPropertyInterface().getClassName()) +
                     "[" + getFullID().asString() + "]: For <= -10" +
                      " coefficients, only nonHD species allowed.");
                }
              else
                {
                  Species* species(theSpatiocyteStepper->getSpecies(aVariable));
                  if(aCoefficient == -10)
                    {
                      theAdjoinSubstratesPreA.push_back(species->getID());
                    }
                  else if(aCoefficient == -20)
                    {
                      theAdjoinSubstratesPreB.push_back(species->getID());
                    }
                  if(aCoefficient == -11)
                    {
                      theAdjoinSubstratesPostA.push_back(species->getID());
                    }
                  else if(aCoefficient == -21)
                    {
                      theAdjoinSubstratesPostB.push_back(species->getID());
                    }
                }
            }
          else
            {
              theOrder -= aCoefficient; 
              //The first reactant, A:
              if(A == NULL && variableA == NULL)
                {
                  coefficientA = aCoefficient;
                  if(aVariable->getName() == "HD")
                    {
                      variableA = aVariable;
                    }
                  else
                    {
                      A = theSpatiocyteStepper->getSpecies(aVariable);
                    }
                }
              //The second reactant, B:
              else if(B == NULL && variableB == NULL)
                {
                  coefficientB = aCoefficient;
                  if(aVariable->getName() == "HD")
                    {
                      variableB = aVariable;
                    }
                  else
                    {
                      B = theSpatiocyteStepper->getSpecies(aVariable);
                    }
                }
              //The third reactant, G:
              else
                {
                  coefficientG = aCoefficient;
                  if(aVariable->getName() == "HD")
                    {
                      variableG = aVariable;
                    }
                  else
                    {
                      G = theSpatiocyteStepper->getSpecies(aVariable);
                    }
                }
            }
        }
      else if(aCoefficient > 0)
        {
          if(aCoefficient >= 10)
            {
              if(aVariable->getName() == "HD")
                {
                  THROW_EXCEPTION(ValueError,
                     String(getPropertyInterface().getClassName()) +
                     "[" + getFullID().asString() + "]: For >= 10" +
                      " coefficients, only nonHD species allowed.");
                }
              else
                {
                  Species* species(theSpatiocyteStepper->getSpecies(aVariable));
                  if(aCoefficient == 10)
                    {
                      theAdjoinProductsPreA.push_back(species->getID());
                    }
                  else if(aCoefficient == 20)
                    {
                      theAdjoinProductsPreB.push_back(species->getID());
                    }
                  if(aCoefficient == 11)
                    {
                      theAdjoinProductsPostA.push_back(species->getID());
                    }
                  else if(aCoefficient == 21)
                    {
                      theAdjoinProductsPostB.push_back(species->getID());
                    }
                }
            }
          else
            {
              //The first product, C:
              if(C == NULL && variableC == NULL)
                {
                  coefficientC = aCoefficient;
                  if(aVariable->getName() == "HD")
                    {
                      variableC = aVariable;
                    }
                  else
                    {
                      C = theSpatiocyteStepper->getSpecies(aVariable);
                    }
                }
              //The second product, D:
              else if(D == NULL && variableD == NULL)
                {
                  coefficientD = aCoefficient;
                  if(aVariable->getName() == "HD")
                    {
                      variableD = aVariable;
                    }
                  else
                    {
                      D = theSpatiocyteStepper->getSpecies(aVariable);
                    }
                }
              //The third product, F:
              else if(F == NULL && variableF == NULL)
                {
                  coefficientF = aCoefficient;
                  if(aVariable->getName() == "HD")
                    {
                      variableF = aVariable;
                    }
                  else
                    {
                      F = theSpatiocyteStepper->getSpecies(aVariable);
                    }
                }
            }
        }
      //aCoefficient == 0:
      else
        {
          //The first non-changed species, E:
          if(E == NULL && variableE == NULL)
            {
              coefficientE = aCoefficient;
              if(aVariable->getName() == "HD")
                {
                  variableE = aVariable;
                }
              else
                {
                  E = theSpatiocyteStepper->getSpecies(aVariable);
                }
            }
          //The second non-changed species, H:
          else if(H == NULL && variableH == NULL)
            {
              coefficientH = aCoefficient;
              if(aVariable->getName() == "HD")
                {
                  variableH = aVariable;
                }
              else
                {
                  H = theSpatiocyteStepper->getSpecies(aVariable);
                }
            }
        }
    }
  std::cout << "theOrder:" << theOrder << std::endl;
} 


void ReactionProcess::logEvent()
{
  if(LogEvent && theSpatiocyteStepper->getCurrentTime() >= LogStart)
    {
      if(!theEventCnt)
        {
          if(FileName == "LogEvent.csv")
            {
              FileName = String(getFullID().getID()) + String(".csv");
              std::cout << "LogEvent for reaction " << getIDString() << ":" <<
                FileName << std::endl;
            }
          theLogFile.open(FileName.c_str(), std::ios::trunc);
          theLogFile << "CurrentTime-LogStart,Events,";
          String strA;
          String strB;
          if(variableA)
            {
              strA = getIDString(variableA);
            }
          else if(A)
            {
              strA = getIDString(A);
            }
          if(variableB)
            {
              strB = getIDString(variableB);
            }
          else if(B)
            {
              strB = getIDString(B);
            }
          if(!strB.empty())
            {
              theLogFile << strA << "," << strB << 
                ",Rate=Events/((CurrentTime-LogStart)*" << strA << "*" <<
                strB << std::endl;
            }
          else
            {
              theLogFile << strA << ",Rate=Events/((CurrentTime-LogStart)*" <<
                strA << std::endl;
            }
        }
      unsigned numberA(0);
      unsigned numberB(0);
      if(variableA)
        {
          numberA = variableA->getValue()+1;
        }
      else if(A)
        {
          numberA = A->getVariable()->getValue()+1;
          if(A->getIsCompVacant())
            {
              for(unsigned i(0); i != theSpecies.size(); ++i)
                {
                  if(theSpecies[i] != A &&
                     theSpecies[i]->getVacantSpecies() == A)
                    {
                      numberA -= theSpecies[i]->size();
                    }
                }
            }
        }
      if(variableB)
        {
          numberB = variableB->getValue()+1;
        }
      else if(B)
        {
          numberB = B->getVariable()->getValue()+1;
          if(B->getIsCompVacant())
            {
              for(unsigned i(0); i != theSpecies.size(); ++i)
                {
                  if(theSpecies[i] != B && 
                     theSpecies[i]->getVacantSpecies() == B)
                    {
                      numberB -= theSpecies[i]->size();
                    }
                }
            }
        }
      theLogFile << theSpatiocyteStepper->getCurrentTime()-LogStart<< "," <<
        ++theEventCnt << "," << numberA;
      if(variableB || B)
        {
          theLogFile << "," << numberB << "," << theEventCnt/
            (theSpatiocyteStepper->getCurrentTime()-LogStart)/(numberA*numberB);
        }
      else
        {
          theLogFile << "," << theEventCnt/
            (theSpatiocyteStepper->getCurrentTime()-LogStart)/(numberA);
        }
      theLogFile << std::endl;
    }
}

}
