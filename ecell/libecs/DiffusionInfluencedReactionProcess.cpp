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

#include <libecs/DiffusionInfluencedReactionProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(DiffusionInfluencedReactionProcess, Process);

void DiffusionInfluencedReactionProcess::checkSubstrates()
{
  //HD_A or HD_B:
	if(variableA)
    {
      THROW_EXCEPTION(ValueError, String(
        getPropertyInterface().getClassName()) + " [" + getFullID().asString() +
		    "]: A DiffusionInfluencedReactionProcess cannot have a HD " +
        "substrate species: " + getIDString(variableA));
    }
  if(variableB)
    {
      THROW_EXCEPTION(ValueError, String(
        getPropertyInterface().getClassName()) + " [" + getFullID().asString() +
		    "]: A DiffusionInfluencedReactionProcess cannot have a HD " +
        "substrate species: " + getIDString(variableB));
    }
}

void DiffusionInfluencedReactionProcess::initializeSecond()
{
  ReactionProcess::initializeSecond(); 
  A->setCollision(Collision);
  B->setCollision(Collision);
}

void DiffusionInfluencedReactionProcess::initializeThird()
{
  ReactionProcess::initializeThird();
  A->setDiffusionInfluencedReactantPair(B); 
  B->setDiffusionInfluencedReactantPair(A); 
  r_v = theSpatiocyteStepper->getVoxelRadius();
  D_A = A->getDiffusionCoefficient();
  D_B = B->getDiffusionCoefficient();
  calculateReactionProbability();
  if(A->getIsDiffusing())
    {
      A->setDiffusionInfluencedReaction(this, B->getID(), p); 
    }
  if(B->getIsDiffusing())
    {
      B->setDiffusionInfluencedReaction(this, A->getID(), p); 
    }
}

//Do the reaction A + B -> C + D. So that A <- C and B <- D.
//We need to consider that the source molecule can be either A or B.
//If A and C belong to the same Comp, A <- C.
//Otherwise, find a vacant adjoin voxel of A, X which is the same Comp
//as C and X <- C.
//Similarly, if B and D belong to the same Comp, B <- D.
//Otherwise, find a vacant adjoin voxel of C, Y which is the same Comp
//as D and Y <- D.
bool DiffusionInfluencedReactionProcess::react(unsigned molA, unsigned molB,
                                               unsigned indexA, unsigned indexB)
{
  /*
  moleculeA = molA;
  moleculeB = molB;
  //nonHD_A + nonHD_B -> nonHD_C + HD_D:
  //nonHD_A + nonHD_B -> HD_C + nonHD_D:
  if((variableC && D) || (C && variableD))
    {
      Variable* HD_p(variableC);
      Species* nonHD_p(D);
      if(variableD)
        {
          HD_p = variableD;
          nonHD_p = C;
        }
      if(A->isReplaceable(moleculeA, nonHD_p))
        {
          moleculeP = moleculeA;
          //Hard remove the B molecule, since nonHD_p is in a different Comp:
          (*theIDs)[moleculeB] = B->getVacantID();
        }
      else if(B->isReplaceable(moleculeB, nonHD_p))
        {
          moleculeP = moleculeB;
          //Hard remove the A molecule, since nonHD_p is in a different Comp:
          (*theIDs)[moleculeA] = A->getVacantID();
        }
      else
        { 
          moleculeP = nonHD_p->getRandomAdjoin(moleculeA, SearchVacant);
          //Only proceed if we can find an adjoin vacant voxel
          //of A which can be occupied by C:
          if(moleculeP == theNullMol)
            {
              moleculeP = nonHD_p->getRandomAdjoin(moleculeB, SearchVacant);
              if(moleculeP == theNullMol)
                {
                  return false;
                }
            }
          //Hard remove the A molecule, since nonHD_p is in a different Comp:
          (*theIDs)[moleculeA] = A->getVacantID();
          //Hard remove the B molecule, since nonHD_p is in a different Comp:
          (*theIDs)[moleculeB] = B->getVacantID();
        }
      HD_p->addValue(1);
      nonHD_p->addMol(moleculeP, A->getTag(indexA));
      return true;
    }
  //nonHD_A + nonHD_B -> HD_C:
  else if(variableC && !D && !variableD)
    {

      //Hard remove the A molecule, since nonHD_p is in a different Comp:
      (*theIDs)[moleculeA] = A->getVacantID();
      //Hard remove the B molecule, since nonHD_p is in a different Comp:
      (*theIDs)[moleculeB] = B->getVacantID();
      variableC->addValue(1);
      return true;
    }

  if(A->isReplaceable(moleculeA, C))
    {
      moleculeC = moleculeA;
      if(D)
        {
          if(B->isReplaceable(moleculeB, D))
            {
              moleculeD = moleculeB;
            }
          else
            {
              moleculeD = D->getRandomAdjoin(moleculeC, moleculeC,
                                             SearchVacant);
              if(moleculeD == theNullMol)
                {
                  return false;
                }
              (*theIDs)[moleculeB] = B->getVacantID();
            }
          D->addMol(moleculeD, B->getTag(indexB));
        }
      else
        {
          //Hard remove the B molecule since it is not used:
          (*theIDs)[moleculeB] = B->getVacantID();
        }
    }
  else if(B->isReplaceable(moleculeB, C))
    {
      moleculeC = moleculeB;
      if(D)
        {
          if(A->isReplaceable(moleculeA, D))
            {
              moleculeD = moleculeA;
            }
          else
            {
              moleculeD = D->getRandomAdjoin(moleculeC, moleculeC,
                                             SearchVacant);
              if(moleculeD == theNullMol)
                {
                  return false;
                }
              (*theIDs)[moleculeA] = A->getVacantID();
            }
          D->addMol(moleculeD, B->getTag(indexB));
        }
      else
        {
          //Hard remove the A molecule since it is not used:
          (*theIDs)[moleculeA] = A->getVacantID();
        }
    }
  else
    {
      moleculeC = C->getRandomAdjoin(moleculeA, SearchVacant);
      if(moleculeC == theNullMol)
        {
          moleculeC = C->getRandomAdjoin(moleculeB, SearchVacant);
          if(moleculeC == theNullMol)
            {
              //Only proceed if we can find an adjoin vacant voxel
              //of A or B which can be occupied by C:
              return false;
            }
        }
      if(D)
        {
          moleculeD = D->getRandomAdjoin(moleculeC, moleculeC,
                                                 SearchVacant);
          if(moleculeD == theNullMol)
            {
              return false;
            }
          D->addMol(moleculeD, B->getTag(indexB));
        }
      //Hard remove the A molecule since it is not used:
      (*theIDs)[moleculeA] = A->getVacantID();
      //Hard remove the B molecule since it is not used:
      (*theIDs)[moleculeB] = B->getVacantID();
    }
  C->addMol(moleculeC, A->getTag(indexA));
  addMolE();
  addMolF();
  */
  return true;
}

//positive-coefficient F
void DiffusionInfluencedReactionProcess::addMolF()
{
  /*
  if(!F)
    {
      return;
    }
  moleculeF = F->getRandomAdjoin(moleculeC, SearchVacant);
  if(moleculeF == theNullMol)
    {
      moleculeF = F->getRandomAdjoin(moleculeD, SearchVacant);
      if(moleculeF == theNullMol)
        {
          return;
        }
    }
  F->addMol(moleculeF);
  */
}

//zero-coefficient E
//we create a molecule E at random location in the compartment to avoid
//rebinding effect, useful to maintain the concentration of a substrate species
//even after the reaction:
void DiffusionInfluencedReactionProcess::addMolE()
{
  /*
  if(!E)
    {
      return;
    }
  moleculeE = E->getRandomCompMol(1);
  if(moleculeE == theNullMol)
    {
      std::cout << getFullID().asString() << " unable to add molecule E" <<
        std::endl;
      return;
    }
  E->addMol(moleculeE);
  */
}

void DiffusionInfluencedReactionProcess::finalizeReaction()
{
  //The number of molecules may have changed for both reactant and product
  //species. We need to update SpatiocyteNextReactionProcesses which are
  //dependent on these species:
  for(std::vector<SpatiocyteProcess*>::const_iterator 
      i(theInterruptedProcesses.begin());
      i!=theInterruptedProcesses.end(); ++i)
    {
      (*i)->substrateValueChanged(theSpatiocyteStepper->getCurrentTime());
    }
}

void DiffusionInfluencedReactionProcess::calculateReactionProbability()
{
  //Refer to the paper for the description of the variables used in this
  //method.
  if(A->getDimension() == 3 && B->getDimension() == 3)
    {
      if(A != B)
        {
          if(p == -1)
            {
              p = k/(6*sqrt(2)*(D_A+D_B)*r_v);
            }
          else
            {
              k = p*(6*sqrt(2)*(D_A+D_B)*r_v);
            }
        }
      else
        {
          if(p == -1)
            {
              p = k/(6*sqrt(2)*D_A*r_v);
            }
          else
            {
              k = p*(6*sqrt(2)*D_A*r_v);
            }
        }
    }
  else if(A->getDimension() != 3 && B->getDimension() != 3)
    {
      //Inter-surface Comp reaction.
      //For surface edge absorbing reactions:
      if(A->getComp() != B->getComp())
        {
          k = p;
        }
      else if(A != B)
        {
          if(p == -1)
            {
              p = pow(2*sqrt(2)+4*sqrt(3)+3*sqrt(6)+sqrt(22), 2)*k/
                (72*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))*(D_A+D_B));
            }
          else
            {
              k = p*(72*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))*(D_A+D_B))/
                pow(2*sqrt(2)+4*sqrt(3)+3*sqrt(6)+sqrt(22), 2);
            }
        }
      else
        {
          if(p == -1)
            {
              p = pow(2*sqrt(2)+4*sqrt(3)+3*sqrt(6)+sqrt(22), 2)*k/
                (72*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))*(D_A));
            }
          else
            {
              k = p*(72*(6*sqrt(2)+4*sqrt(3)+3*sqrt(6))*(D_A))/
                pow(2*sqrt(2)+4*sqrt(3)+3*sqrt(6)+sqrt(22), 2);
            }
        }
    }
  else if(A->getDimension() == 3 && B->getIsLipid())
    {
      if(p == -1)
        {
          p = 24*k*r_v/((6+3*sqrt(3)+2*sqrt(6))*D_A);
        }
      else
        {
          k = p*((6+3*sqrt(3)+2*sqrt(6))*D_A)/(24*r_v);
        }
    }
  else if(A->getIsLipid() && B->getDimension() == 3)
    {
      if(p == -1)
        {
          p = 24*k*r_v/((6+3*sqrt(3)+2*sqrt(6))*D_B);
        }
      else
        {
          k = p*((6+3*sqrt(3)+2*sqrt(6))*D_B)/(24*r_v);
        }
    }
  else if(A->getDimension() == 3 && B->getDimension() != 3)
    {
      if(p == -1)
        {
          p = sqrt(2)*k/(3*D_A*r_v);
        }
      else
        {
          k = p*(3*D_A*r_v)/sqrt(2);
        }
    }
  else if(A->getDimension() != 3 && B->getDimension() == 3)
    {
      if(p == -1)
        {
          p = sqrt(2)*k/(3*D_B*r_v);
        }
      else
        {
          k = p*(3*D_B*r_v)/sqrt(2);
        }
    }
  else
    {
      THROW_EXCEPTION(ValueError, 
                      String(getPropertyInterface().getClassName()) + 
                      " [" + getFullID().asString() + 
                      "]: Error in type of second order reaction.");
    }
}

void DiffusionInfluencedReactionProcess::printParameters()
{
  String aProcess(String(getPropertyInterface().getClassName()) + 
                                      "[" + getFullID().asString() + "]");
  std::cout << aProcess << std::endl;
  std::cout << "  " << getIDString(A) << " + " <<  getIDString(B) << " -> ";
  if(C)
    {
      std::cout << getIDString(C);
    }
  else
    {
      std::cout << getIDString(variableC);
    }
  if(D)
    {
      std::cout << " + " << getIDString(D);
    }
  else if(variableD)
    {
      std::cout << " + " << getIDString(variableD);
    }
  std::cout << ": k=" << k << ", p=" << p << 
    ", p_A=" << A->getReactionProbability(B->getID()) <<
    ", p_B=" << B->getReactionProbability(A->getID()) << std::endl; 
}

}
