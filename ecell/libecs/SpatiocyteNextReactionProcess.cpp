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

#include <libecs/SpatiocyteNextReactionProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>
#include <libecs/ReactionProcess.hpp>

namespace libecs
{

LIBECS_DM_INIT_STATIC(SpatiocyteNextReactionProcess, Process);

void SpatiocyteNextReactionProcess::fire()
{
  if(A)
    {
      A->updateMols();
    }
  if(B)
    {
      B->updateMols();
    }
  if(theOrder == 0)
    {
      if(C)
        { 
          moleculeC = C->getRandomCompMol(SearchVacant);
          if(moleculeC == theNullMol)
            {
              requeue();
              return;
            }
          C->addMol(moleculeC);
        }
      else if(variableC)
        {
          variableC->addValue(coefficientC);
        }
    }
  else if(theOrder == 1)
    { 
      //nonHD_A -> nonHD_C + nonHD_D:
      if(A && C && D)
        {
          if(BindingSite == -1)
            {
              if(!reactACD(A, C, D))
                {
                  return;
                }
            }
          else if(!reactACDbind(A, C, D))
            {
              return;
            }
        }
      //nonHD_A -> nonHD_C:
      else if(A && C && !D && !variableD)
        {
          if(BindingSite == -1)
            {
              if(!reactAC(A, C))
                {
                  return;
                }
            }
          else if(!reactACbind(A, C))
            {
              return;
            }
        }
      //nonHD_A -> HD_C + HD_D:
      else if(A && variableC && variableD)
        {
          moleculeA = A->getRandomMol();
          A->removeMol(moleculeA);
          variableC->addValue(coefficientC);
          variableD->addValue(coefficientD);
        }
      //nonHD_A -> HD_C:
      else if(A && variableC && !D && !variableD)
        {
          moleculeA = A->getRandomMol();
          A->removeMol(moleculeA);
          variableC->addValue(coefficientC);
        }
      //nonHD_A -> nonHD_C + HD_D:
      //nonHD_A -> HD_C + nonHD_D:
      else if(A && ((variableC && D) || (C && variableD)))
        {
          Variable* HD_p(variableC);
          int coefficient(coefficientC);
          Species* nonHD_p(D);
          if(variableD)
            {
              HD_p = variableD;
              coefficient = coefficientD;
              nonHD_p = C;
            }
          if(reactAC(A, nonHD_p))
             {
               HD_p->addValue(coefficient);
             }
          else
            {
              return;
            }
        }
      //HD_A -> nonHD_C:
      else if(variableA && C && !D && !variableD)
        {
          moleculeC = reactvAC(variableA, C);
          if(moleculeC == theNullMol)
            {
              requeue();
              return;
            }
          else
            {
              variableA->addValue(coefficientA);
              C->addMol(moleculeC);
            }
        }
      //HD_A -> nonHD_C + nonHD_D:
      else if(variableA && C && D)
        {
          moleculeC = theNullMol;
          moleculeD = theNullMol;
          Comp* compA(theSpatiocyteStepper->system2Comp(
                         variableA->getSuperSystem()));
          //Occupy C in a voxel of compartment C that adjoins compartment A
          //if A is a surface compartment:
          if(compA != C->getComp() && compA->dimension != 3)
            {
              moleculeC = C->getRandomAdjoinCompMol(compA, SearchVacant);
              if(moleculeC)
                {
                  moleculeD = D->getRandomAdjoin(moleculeC, moleculeC,
                                                         SearchVacant);
                }
            }
          else if(compA != D->getComp() && compA->dimension != 3)
            {
              moleculeD = D->getRandomAdjoinCompMol(compA, SearchVacant);
              if(moleculeD)
                {
                  moleculeC = C->getRandomAdjoin(moleculeD, moleculeD,
                                                         SearchVacant);
                }
            }
          else
            {
              moleculeC = C->getRandomCompMol(SearchVacant);
              if(moleculeC)
                {
                  moleculeD = D->getRandomAdjoin(moleculeC, moleculeC,
                                                         SearchVacant);
                }
            }
          if(moleculeC == theNullMol || moleculeD == theNullMol)
            {
              requeue();
              return;
            }
          variableA->addValue(coefficientA);
          C->addMol(moleculeC);
          D->addMol(moleculeD);
        }
      //HD_A -> HD_C + HD_D:
      else if(variableA && variableC && variableD)
        {
          variableA->addValue(coefficientA);
          variableC->addValue(coefficientC);
          variableD->addValue(coefficientD);
        }
      //HD_A -> HD_C:
      else if(variableA && variableC && !D && !variableD)
        {
          variableA->addValue(coefficientA);
          variableC->addValue(coefficientC);
        }
      //HD_A -> nonHD_C + HD_D:
      //HD_A -> HD_C + nonHD_D:
      else if(variableA && ((variableC && D) || (C && variableD)))
        {
          Variable* HD_p(variableC);
          int coefficient(coefficientC);
          Species* nonHD_p(D);
          if(variableD)
            {
              HD_p = variableD;
              coefficient = coefficientD;
              nonHD_p = C;
            }
          moleculeP = reactvAC(variableA, nonHD_p);
          if(moleculeP == theNullMol)
            {
              requeue();
              return;
            }
          variableA->addValue(coefficientA);
          nonHD_p->addMol(moleculeP);
          HD_p->addValue(coefficient);
        }
    }
  //number of substrate species = 2:
  //coefficients could be more.
  else
    {
      //HD + HD -> product(s)
      if(variableA && variableB)
        {
          //HD + HD -> HD: 
          if(variableC && !variableD && !D)
            {
              variableA->addValue(coefficientA);
              variableB->addValue(coefficientB);
              variableC->addValue(coefficientC);
            }
          //HD + HD -> nonHD: 
          else if(C && !variableD && !D)
            { 
              moleculeC = reactvAvBC(C);
              if(moleculeC == theNullMol)
                {
                  requeue();
                  return;
                }
              variableA->addValue(coefficientA);
              variableB->addValue(coefficientB);
              C->addMol(moleculeC);
            }
          //HD + HD -> HD + HD: 
          else if(variableC && variableD)
            {
              variableA->addValue(coefficientA);
              variableB->addValue(coefficientB);
              variableC->addValue(coefficientC);
              variableD->addValue(coefficientD);
            }
          //HD + HD -> HD + nonHD: 
          //HD + HD -> nonHD + HD: 
          else if((variableC && D) || (C && variableD))
            {
              Variable* HD_p(variableC);
              int coefficient(coefficientC);
              Species* nonHD_p(D);
              if(variableD)
                {
                  HD_p = variableD;
                  coefficient = coefficientD;
                  nonHD_p = C;
                }
              moleculeP = reactvAvBC(nonHD_p);
              if(moleculeP == theNullMol)
                {
                  requeue();
                  return;
                }
              variableA->addValue(coefficientA);
              variableB->addValue(coefficientB);
              nonHD_p->addMol(moleculeP);
              HD_p->addValue(coefficient);
            }
          //HD + HD -> nonHD + nonHD: 
          else if(C && D)
            {
              moleculeC = reactvAvBC(C);
              moleculeD = theNullMol;
              if(moleculeC == theNullMol)
                {
                  moleculeD = reactvAvBC(D);
                  if(moleculeD)
                    {
                      moleculeC = C->getRandomAdjoin(moleculeD,
                                                     moleculeD, SearchVacant);
                    }
                }
              else
                { 
                  moleculeD = D->getRandomAdjoin(moleculeC, moleculeC,
                                                         SearchVacant);
                }
              if(moleculeC == theNullMol || moleculeD == theNullMol)
                {
                  requeue();
                  return;
                }
              variableA->addValue(coefficientA);
              variableB->addValue(coefficientB);
              C->addMol(moleculeC);
              D->addMol(moleculeD);
            }
        }
      //HD + nonHD -> product(s)
      //nonHD + HD -> product(s)
      else if(variableA || variableB)
        {
          Species* nonHD(A);
          Variable* HD(variableB);
          int coefficient(coefficientB);
          if(B)
            {
              nonHD = B;
              HD = variableA;
              coefficient = coefficientA;
            }
          //nonHD + HD -> nonHD + nonHD: 
          //HD + nonHD -> nonHD + nonHD: 
          if(C && D)
            { 
              if(!reactACD(nonHD, C, D))
                {
                  return;
                }
              HD->addValue(coefficient);
            }
          //nonHD + HD -> nonHD:
          //HD + nonHD -> nonHD:
          else if(C && !D && !variableD)
            {
              if(!reactAC(nonHD, C))
                {
                  return;
                }
              HD->addValue(coefficient);
            }
          //nonHD + HD -> HD:
          //HD + nonHD -> HD:
          else if(variableC && !D && !variableD)
            {
              moleculeS = nonHD->getRandomMol();
              nonHD->removeMol(moleculeS);
              HD->addValue(coefficient);
              variableC->addValue(coefficientC);
            }
          //HD + nonHD -> HD + nonHD:
          //HD + nonHD -> nonHD + HD:
          //nonHD + HD -> HD + nonHD:
          //nonHD + HD -> nonHD + HD:
          else if((variableC && D) || (C && variableD))
            {
              Variable* HD_p(variableC);
              int coefficient_p(coefficientC);
              Species* nonHD_p(D);
              if(variableD)
                {
                  HD_p = variableD;
                  coefficient_p = coefficientD;
                  nonHD_p = C;
                }
              if(!reactAC(nonHD, nonHD_p))
                {
                  return;
                }
              HD->addValue(coefficient);
              HD_p->addValue(coefficient_p);
            }
        }
      //nonHD + nonHD -> product(s)
      else
        {
          //nonHD + nonHD -> nonHD + nonHD
          if(C && D)
            {
              reactABCD();
            }
          //nonHD + nonHD -> nonHD
          else
            {
              //multiNonHD + nonHD -> nonHD
              if(A->getIsMultiscale())
                {
                  if(!reactMultiABC())
                    {
                      return;
                    }
                }
            }
        }
    }
  ReactionProcess::fire();
}

//MultiNonHD.nonHD -> nonHD
//MultiA.B -> C
bool SpatiocyteNextReactionProcess::reactMultiABC()
{
  std::cout << "react Multi ABC:" << nextIndexA << std::endl;
  moleculeA = A->getMol(nextIndexA);
  moleculeC = C->getRandomAdjoin(moleculeA, SearchVacant);
  if(moleculeC == theNullMol)
    {
      requeue();
      return false;
    }
  A->removeMolIndex(nextIndexA);
  std::cout << "adding molecule" << std::endl;
  C->addMol(moleculeC);
  return true;
}

//nonHD.nonHD -> nonHD + nonHD
//Both A and B are immobile nonHD
void SpatiocyteNextReactionProcess::reactABCD()
{
  unsigned rand(gsl_rng_uniform_int(getStepper()->getRng(), moleculesA.size()));
  moleculeA = moleculesA[rand];
  moleculeB = A->getRandomAdjoin(moleculeA, B, SearchVacant);
  if(A != C)
    {
      unsigned indexA(A->getMolIndex(moleculeA));
      Tag tagA(A->getTag(indexA));
      A->removeMolIndex(indexA);
      C->addMol(moleculeA, tagA);
    }
  if(B != D)
    { 
      unsigned indexB(B->getMolIndex(moleculeB));
      Tag tagB(B->getTag(indexB));
      B->removeMolIndex(indexB);
      D->addMol(moleculeB, tagB);
    }
}

//nonHD -> nonHD + nonHD
bool SpatiocyteNextReactionProcess::reactACD(Species* a, Species* c, Species* d)
{
  unsigned indexA(a->getRandomIndex());
  moleculeA = a->getMol(indexA);
  moleculeC = theNullMol;
  moleculeD = theNullMol;
  if(a->getVacantID() == c->getVacantID() || a->getID() == c->getVacantID())
    {
      moleculeC = moleculeA;
      moleculeD = d->getRandomAdjoin(moleculeC, moleculeC, SearchVacant);
      if(moleculeD == theNullMol)
        {
          requeue();
          return false;
        }
    }
  else if(a->getVacantID() == d->getVacantID() ||
          a->getID() == d->getVacantID())
    {
      moleculeD = moleculeA;
      moleculeC = c->getRandomAdjoin(moleculeD, moleculeD,
                                             SearchVacant);
      if(moleculeC == theNullMol)
        {
          requeue();
          return false;
        }
    }
  else
    {
      moleculeC = c->getRandomAdjoin(moleculeA, SearchVacant);
      if(moleculeC == theNullMol)
        {
          //Only proceed if we can find an adjoin vacant voxel
          //of nonND which can be occupied by C:
          requeue();
          return false;
        }
      moleculeD = d->getRandomAdjoin(moleculeC, moleculeC,
                                             SearchVacant);
      if(moleculeD == theNullMol)
        {
          requeue();
          return false;
        }
    }
  Tag tagA(a->getTag(indexA));
  a->removeMolIndex(indexA);
  c->addMol(moleculeC, tagA);
  d->addMol(moleculeD);
  return true;
}

//nonHD (+ E) -> nonHD
bool SpatiocyteNextReactionProcess::reactAC(Species* a, Species* c)
{
  unsigned indexA(a->getRandomIndex());
  moleculeA = a->getMol(indexA);
  moleculeC = theNullMol;
  if(a->getVacantID() == c->getVacantID() || a->getID() == c->getVacantID())
    {
      moleculeC = moleculeA;
    }
  else
    {
      moleculeC = c->getRandomAdjoin(moleculeA, SearchVacant);
      if(moleculeC == theNullMol)
        {
          //Only proceed if we can find an adjoin vacant voxel
          //of nonND which can be occupied by C:
          requeue();
          return false;
        }
    }
  Tag tagA(a->getTag(indexA));
  a->removeMolIndex(indexA);
  c->addMol(moleculeC, tagA);
  removeMolE();
  return true;
}

//zero-coefficient E
//we remove a molecule E at random location in the compartment to allow
//rebinding effect of the dissociated nonHD molecule while maintain the 
//concentration of a substrate species even after the reaction:
void SpatiocyteNextReactionProcess::removeMolE()
{
  if(!E)
    {
      return;
    }
  moleculeE = E->getRandomMol();
  if(moleculeE == theNullMol)
    {
      std::cout << getFullID().asString() << " unable to remove molecule E" <<
        std::endl;
    }
  E->removeMol(moleculeE);
}

//nonHD (+Vacant[BindingSite]) -> nonHD[BindingSite]
bool SpatiocyteNextReactionProcess::reactACbind(Species* a, Species* c)
{
  unsigned indexA(a->getRandomIndex());
  moleculeA = a->getMol(indexA);
  moleculeC = c->getBindingSiteAdjoin(moleculeA, BindingSite);
  if(moleculeC == theNullMol)
    {
      //Only proceed if we can find an adjoin vacant voxel
      //of nonND which can be occupied by C:
      requeue();
      return false;
    }
  Tag tagA(a->getTag(indexA));
  a->removeMolIndex(indexA);
  c->addMol(moleculeC, tagA);
  return true;
}

//nonHD(a) -> nonHD(c[BindingSite]) if vacant[BindingSite]
//else if BindingSite is not vacant:
//nonHD(a) -> nonHD(d)
bool SpatiocyteNextReactionProcess::reactACDbind(Species* a, Species* c,
                                                 Species* d)
{
  unsigned indexA(a->getRandomIndex());
  moleculeA = a->getMol(indexA);
  moleculeC = c->getBindingSiteAdjoin(moleculeA, BindingSite);
  if(moleculeC == theNullMol)
    {
      moleculeD = theNullMol;
      moleculeD = d->getRandomAdjoin(moleculeA, SearchVacant);
      if(moleculeD == theNullMol)
        {
          requeue();
          return false;
        }
      Tag tagA(a->getTag(indexA));
      a->removeMolIndex(indexA);
      d->addMol(moleculeD, tagA);
      return true;
    }
  Tag tagA(a->getTag(indexA));
  a->removeMolIndex(indexA);
  c->addMol(moleculeC, tagA);
  return true;
}

//HD -> nonHD
unsigned SpatiocyteNextReactionProcess::reactvAC(Variable* vA, Species* c)
{
  moleculeC = theNullMol;
  Comp* compA(theSpatiocyteStepper->system2Comp(vA->getSuperSystem()));
  //Occupy C in a voxel of compartment C that adjoins compartment A
  //if A is a surface compartment:
  if(compA != c->getComp() && compA->dimension != 3)
    {
      moleculeC = c->getRandomAdjoinCompMol(compA, SearchVacant);
    }
  else
    {
      moleculeC = c->getRandomCompMol(SearchVacant);
    }
  return moleculeC;
}

Comp* SpatiocyteNextReactionProcess::getComp2D(Species* c)
{
  Comp* compA(theSpatiocyteStepper->system2Comp(variableA->getSuperSystem()));
  Comp* compB(theSpatiocyteStepper->system2Comp(variableB->getSuperSystem()));
  Comp* comp2D(NULL);
  if(compA->dimension == 2)
    {
      comp2D = compA;
    }
  else if(compB->dimension == 2)
    {
      comp2D = compB;
    }
  //Occupy C in a voxel of compartment C that adjoins compartment A
  //if A is a surface compartment:
  if(comp2D != c->getComp() && comp2D != NULL)
    {
      return comp2D;
    }
  return NULL;
}

unsigned SpatiocyteNextReactionProcess::reactvAvBC(Species* c)
{
  moleculeC = theNullMol;
  Comp* aComp2D(getComp2D(c));
  if(aComp2D)
    {
      moleculeC = C->getRandomAdjoinCompMol(aComp2D, SearchVacant);
    }
  else
    {
      moleculeC = C->getRandomCompMol(SearchVacant);
    }
  return moleculeC;
}


Real SpatiocyteNextReactionProcess::getPropensity_ZerothOrder() 
{
  return p;
}

Real SpatiocyteNextReactionProcess::getPropensity_FirstOrder() 
{
  if(A)
    {
      A->updateMolSize();
    }
  Real aValue(theVariableReferenceVector[0].getVariable()->getValue());
  if(aValue > 0.0)
    {
      //std::cout << "p1:" << p << " v:" << aValue << std::endl;
      return p*aValue;
    }
  else
    {
      //std::cout << "p1:0" << std::endl;
      return 0.0;
    }
}

//Need to solve homodimerization reaction of two substrate species (Size-1):
Real SpatiocyteNextReactionProcess::getPropensity_SecondOrder_TwoSubstrates() 
{
  if(A)
    {
      A->updateMolSize();
    }
  if(B)
    {
      B->updateMolSize();
    }
  double sizeA(theVariableReferenceVector[0].getVariable()->getValue());
  double sizeB(theVariableReferenceVector[1].getVariable()->getValue());
  //Required for HD species when |coefficient| is != 1
  if(variableA)
    {
      sizeA = pow(sizeA, sqrt(coefficientA*coefficientA));
    }
  if(variableB)
    {
      sizeB = pow(sizeB, sqrt(coefficientB*coefficientB));
    }
  if(sizeA > 0.0 && sizeB > 0.0)
    {
      return p*sizeA*sizeB;
    }
  else
    {
      return 0.0;
    }
}

double SpatiocyteNextReactionProcess::getIntervalUnbindMultiAB()
{
  if(!A->size() || !p)
    {
      return libecs::INF;
    }
  nextIndexA = 0;
  double fraction(A->getMultiscaleBoundFraction(nextIndexA,
                                            B->getVacantSpecies()->getID())); 
  double rand(gsl_rng_uniform_pos(getStepper()->getRng()));
  double denom((p*fraction)*(-log(rand)));
  double nextInterval(libecs::INF);
  if(denom)
    {
      nextInterval = 1.0/denom;
    }
  for(unsigned i(1); i != A->size(); ++i)
    {
      fraction = A->getMultiscaleBoundFraction(i,
                                               B->getVacantSpecies()->getID());
      rand = gsl_rng_uniform_pos(getStepper()->getRng());
      denom = (p*fraction)*(-log(rand));
      double interval(libecs::INF);
      if(denom)
        {
          interval = 1.0/denom;
        }
      if(interval < nextInterval)
        {
          nextIndexA = i;
          nextInterval = interval;
        }
    }
  return nextInterval;
}

double SpatiocyteNextReactionProcess::getIntervalUnbindAB()
{
  A->updateMolSize();
  B->updateMolSize();
  if(A->getIsMultiscale())
    {
      return getIntervalUnbindMultiAB();
    }
  const double sizeB(updateSizesAB());
  const double sizeA(moleculesA.size());
  const double rand(gsl_rng_uniform_pos(getStepper()->getRng()));
  const double denom((p*sizeA*sizeB)*(-log(rand)));
  if(denom)
    {
      return 1.0/denom;
    }
  return libecs::INF;
}

unsigned SpatiocyteNextReactionProcess::updateSizesAB()
{
  unsigned sizeB(0);
  moleculesA.resize(0);
  for(unsigned i(0); i != A->size(); ++i)
    {
      moleculeA = A->getMol(i);
      unsigned cnt(A->getAdjoinMolCnt(moleculeA, B));
      if(cnt)
        {
          moleculesA.push_back(moleculeA);
          sizeB += cnt;
        }
    }
  return sizeB;
}

Real SpatiocyteNextReactionProcess::getPropensity_SecondOrder_OneSubstrate() 
{
  if(A)
    {
      A->updateMolSize();
    }
  Real aValue(theVariableReferenceVector[0].getVariable()->getValue());
  //There must be two or more molecules:
  if(aValue > 1.0)
    {
      return p*aValue*(aValue-1.0);
    }
  else
    {
      return 0.0;
    }
}

void SpatiocyteNextReactionProcess::initializeSecond()
{
  ReactionProcess::initializeSecond();
  //if second order, with A and B substrates, both of them
  //must be immobile, or A must be multiscale:
  if (A && B)
    {
      if(A->getDiffusionCoefficient() && !A->getIsMultiscale())
        { 
          THROW_EXCEPTION(ValueError, String(
                            getPropertyInterface().getClassName()) +
                            "[" + getFullID().asString() + 
                            "]: A SpatiocyteNextReactionProcess can have two " +
                            "nonHD substrates (second order) only when both " +
                            "of the species are immobile, or substrate is " +
                            "multiscale and substrate B diffuses within it. " +
                            "However, " + getIDString(A) + " has nonzero " +
                            "diffusion coefficient and not a multiscale " +
                            "species. Use DiffusionInfluencedReaction " +
                            "instead.");
        }
      if(B->getDiffusionCoefficient() && !A->getIsMultiscale())
        {
          THROW_EXCEPTION(ValueError, String(
                            getPropertyInterface().getClassName()) +
                            "[" + getFullID().asString() + 
                            "]: A SpatiocyteNextReactionProcess can have two " +
                            "nonHD substrates (second order) only when both " +
                            "of the species are immobile, or substrate is " +
                            "multiscale and substrate B diffuses within it. " +
                            "However, " + getIDString(B) + " has nonzero " +
                            "diffusion coefficient and not a multiscale " +
                            "species. Use DiffusionInfluencedReaction " +
                            "instead.");
        }
    }
}

//Cannot put the setIsReactiveVacant in initializeSecond because some
//species will only be initialized as vacant in the initializeSecond method
//of other processes (eg. MicrotubuleProcess).
void SpatiocyteNextReactionProcess::initializeThird()
{
  ReactionProcess::initializeThird();
  if(A)
    {
      if(A->getIsVacant())
        {
          A->setIsReactiveVacant();
        }
    }
  else if(variableA)
    {
      variableA->setValue(initSizeA);
    }
  if(B)
    {
      if(B->getIsVacant())
        {
          B->setIsReactiveVacant();
        }
    }
  else if(variableB)
    {
      variableB->setValue(initSizeB);
    }
  if(variableC)
    {
      variableC->setValue(initSizeC);
    }
  if(variableD)
    {
      variableD->setValue(initSizeD);
    }
}

void SpatiocyteNextReactionProcess::initializeFourth()
{
  ReactionProcess::initializeFourth();
  if(p != -1)
    {
      return;
    }
  Comp* compA(NULL);
  Comp* compB(NULL);
  Comp* compC(NULL);
  Comp* compD(NULL);
  if(A)
    {
      compA = A->getComp();
    }
  else if(variableA)
    {
      compA = theSpatiocyteStepper->system2Comp(
                         variableA->getSuperSystem());
    }
  if(B)
    {
      compB = B->getComp();
    }
  else if(variableB)
    {
      compB = theSpatiocyteStepper->system2Comp(
                         variableB->getSuperSystem());
    }
  if(C)
    {
      compC = C->getComp();
    }
  else if(variableC)
    {
      compC = theSpatiocyteStepper->system2Comp(
                         variableC->getSuperSystem());
    }
  if(D)
    {
      compD = D->getComp();
    }
  else if(variableD)
    {
      compD = theSpatiocyteStepper->system2Comp(
                         variableD->getSuperSystem());
    }
  double aVolume(0);
  double anArea(0);
  if(theOrder == 0)
    {
      double aSpace(0);
      if(SpaceC > 0)
        {
          aSpace = SpaceC;
          pFormula << "[aSpace:SpaceC:" << aSpace << "]";
        }
      else if(compC->dimension == 2)
        {
          aSpace = compC->actualArea;
          pFormula << "[aSpace:compC.Area:" << aSpace << "]";
        }
      else
        {
          aSpace = compC->actualVolume;
          pFormula << "[aSpace:compC.Volume:" << aSpace << "]";
        }
      p = k*aSpace;
      pFormula << "[k*aSpace:" << k << "*" << aSpace << "]";
    }
  //Used also by A + B -> product(s) since A and B are in bound form, which
  //is a first order dissociation reaction:
  else if(theOrder == 1 || (A && B)) 
    {
      //Convert the unit m/s of k to 1/s for p if the reaction is a surface
      //adsorption reaction:
      if(compA->dimension == 3 && compC->dimension == 2)
        { 
          if(SpaceA > 0)
            {
              aVolume = SpaceA;
              pFormula << "[aVolume:SpaceA:" << aVolume << "]";
            }
          else
            {
              aVolume = compA->actualVolume;
              pFormula << "[aVolume:compA.Volume:" << aVolume << "]";
            }
          if(SpaceC > 0)
            {
              anArea = SpaceC;
              pFormula << "[anArea:SpaceC:" << anArea << "]";
            }
          else
            {
              anArea = compC->actualArea;
              pFormula << "[anArea:compC.Area:" << anArea << "]";
            }
          p = k*anArea/aVolume;
          pFormula << "[k*anArea/aVolume:" << k << "*" << anArea << "/"
            << aVolume << "]";
          return;
        }
      p = k;
      pFormula << "[k:" << k << "]";
    }
  else
    {
      //If there are two products that don't belong to the same compartment,
      //the reactants must also belong to different compartments:
      if((compD && compD != compC) && (compA == compB))
        {
          NEVER_GET_HERE;
        }
      //If volume + surface <= k(volume)(surface) or
      //   volume + surface <= k(surface)(volume) or
      //   surface + volume <= k(volume)(surface) or
      //   surface + volume <= k(surface)(volume)
      if((compD && (
        (compC->dimension == 3 && compD->dimension == 2 &&
         compA->dimension == 3 && compB->dimension == 2) ||
        (compC->dimension == 3 && compD->dimension == 2 &&
         compA->dimension == 2 && compB->dimension == 3) ||
        (compC->dimension == 2 && compD->dimension == 3 &&
         compA->dimension == 3 && compB->dimension == 2) ||
        (compC->dimension == 2 && compD->dimension == 3 &&
         compA->dimension == 2 && compB->dimension == 3))) ||
      //If volume (+volume) <= k(volume)(volume) or
      //   surface (+surface) <= k(volume)(surface) or
      //   surface (+surface) <= k(surface)(volume)
         ((compC->dimension == 3 && compA->dimension == 3
          && compB->dimension == 3) ||
         (compC->dimension == 2 && compA->dimension == 3 
          && compB->dimension == 2) ||
         (compC->dimension == 2 && compA->dimension == 2 
          && compB->dimension == 3)))
        {
          if(compA->dimension == 3)
            {
              if(SpaceA > 0)
                {
                  aVolume = SpaceA;
                  pFormula << "[aVolume:SpaceA:" << aVolume << "]";
                }
              else
                {
                  aVolume = compA->actualVolume;
                  pFormula << "[aVolume:compA.Volume:" << aVolume << "]";
                }
            }
          else
            {
              if(SpaceB > 0)
                {
                  aVolume = SpaceB;
                  pFormula << "[aVolume:SpaceB:" << aVolume << "]";
                }
              else
                {
                  aVolume = compB->actualVolume;
                  pFormula << "[aVolume:compB.Volume:" << aVolume << "]";
                }
            }
          //unit of k is in m^3/s
          p = k/aVolume;
          pFormula << "[k/aVolume:" << k << "/" << aVolume << "]";
        }
      //If surface (+surface) <= k(surface)(surface) or
      //   volume (+volume) <= k(volume)(surface) or
      //   volume (+volume) <= k(surface)(volume)
      else if((compC->dimension == 2 && compA->dimension == 2 
               && compB->dimension == 2) ||
              (compC->dimension == 3 && compA->dimension == 3 
               && compB->dimension == 2) ||
              (compC->dimension == 3 && compA->dimension == 2 
               && compB->dimension == 3))
        {
          if(compA->dimension == 2)
            {
              if(SpaceA > 0)
                {
                  anArea = SpaceA;
                  pFormula << "[anArea:SpaceA:" << anArea << "]";
                }
              else
                {
                  anArea = compA->actualArea;
                  pFormula << "[anArea:compA.Area:" << anArea << "]";
                }
            }
          else
            {
              if(SpaceB > 0)
                {
                  anArea = SpaceB;
                  pFormula << "[anArea:SpaceB:" << anArea << "]";
                }
              else
                {
                  anArea = compB->actualArea;
                  pFormula << "[anArea:compB.Area:" << anArea << "]";
                }
            }
          //unit of k is in m^2/s
          p = k/anArea;
          pFormula << "[k/anArea:" << k << "/" << anArea << "]";
        }
      else
        {
          NEVER_GET_HERE;
        }
    }
}

void SpatiocyteNextReactionProcess::printParameters()
{
  String aProcess(String(getPropertyInterface().getClassName()) + 
                                      "[" + getFullID().asString() + "]");
  std::cout << aProcess << std::endl;
  if(A)
    {
      std::cout << "  " << getIDString(A);
    }
  else if(variableA)
    {
      std::cout << "  " << getIDString(variableA);
    }
  if(B)
    {
      std::cout << " + " << getIDString(B);
    }
  else if(variableB)
    {
      std::cout << " + " << getIDString(variableB);
    }
  if(!A && !variableA)
    {
      if(C)
        {
          std::cout << "0 -> " << getIDString(C);
        }
      else if(variableC)
        {
          std::cout << "0 -> " << getIDString(variableC);
        }
    }
  else
    {
      if(C)
        {
          std::cout << " -> " << getIDString(C);
        }
      else if(variableC)
        {
          std::cout << " -> " << getIDString(variableC);
        }
    }
  if(D)
    {
      std::cout << " + " << getIDString(D);
    }
  else if(variableD)
    {
      std::cout << " + " << getIDString(variableD);
    }
  std::cout << " k:" << k << " p = " << pFormula.str() << " = " << p
    << " nextTime:" << getInterval() << " propensity:" << getPropensity_R()
    << std::endl;
}

double SpatiocyteNextReactionProcess::getInterval()
{
  if(A && B)
    {
      double interval(getIntervalUnbindAB());
      //std::cout << "interval:" << interval << std::endl;
      return interval;
    }
  double step(getPropensity_R()*(-log(gsl_rng_uniform_pos(getStepper()->getRng()))));
  //std::cout << getFullID().asString() << " " << theTime <<  " next:" << theTime+step << " interval:" << step << std::endl; 
  return step;
}

//Find out if this process is interrupted by the aReactionProcess
//by checking if any of the modified variables of aReactionProcess is a
//substrate of this process:
bool SpatiocyteNextReactionProcess::isInterrupted(ReactionProcess*
                                                  aReactionProcess)
{
  //First get the unique Variables of the aReactionProcess:
  std::vector<Variable*> aVariables;
  VariableReferenceVector
    aVariableReferences(aReactionProcess->getVariableReferenceVector()); 
  for(VariableReferenceVector::iterator i(aVariableReferences.begin());
      i != aVariableReferences.end(); ++i)
    {
      Variable* aVariable((*i).getVariable());
      if(std::find(aVariables.begin(), aVariables.end(), aVariable) ==
         aVariables.end())
        {
          aVariables.push_back(aVariable);
        }
    }
  //Find out if the values of the unique variables will be changed
  //by the ReactionProcess aReactionProcess, i.e., netCoefficient != 0:
  std::vector<int> aNetCoefficients;
  aNetCoefficients.resize(aVariables.size());
  for(std::vector<int>::iterator i(aNetCoefficients.begin());
      i!=aNetCoefficients.end(); ++i)
    {
      (*i) = 0;
    }
  for(VariableReferenceVector::iterator i(aVariableReferences.begin());
      i != aVariableReferences.end(); ++i)
    {
      for(std::vector<Variable*>::const_iterator j(aVariables.begin());
          j != aVariables.end(); ++j)
        {
          if((*i).getVariable() == (*j))
            {
              aNetCoefficients[j-aVariables.begin()] += (*i).getCoefficient();
            }
        }
    }
  //Check if any variable with netCoefficient != 0 is a substrate
  //of this process:
  for(VariableReferenceVector::iterator i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      if((*i).isAccessor())
        {
          Species* aSpecies(theSpatiocyteStepper->variable2species(
                                                 (*i).getVariable()));
          if(aSpecies && aSpecies->getIsMultiscale())
            {
              aSpecies->addInterruptedProcess(this);
            }
          for(std::vector<Variable*>::const_iterator j(aVariables.begin());
              j != aVariables.end(); ++j)
            {
              if((*i).getVariable() == (*j) && 
                 aNetCoefficients[j-aVariables.begin()])
                {
                  return true;
                }
            }
        }
    }
  return false;
}

void SpatiocyteNextReactionProcess::calculateOrder()
{
  ReactionProcess::calculateOrder();
  // set theGetPropensityMethodPtr
  if(theOrder == 0) // no substrate
    {
      theGetPropensityMethodPtr = RealMethodProxy::create<
        &SpatiocyteNextReactionProcess::getPropensity_ZerothOrder>();
    }
  else if(theOrder == 1) // one substrate, first order.
    {
      theGetPropensityMethodPtr = RealMethodProxy::create<
        &SpatiocyteNextReactionProcess::getPropensity_FirstOrder>();
    }
  else
    { 
      //Two unique substrate species:
      //A + B -> products:
      if(getZeroVariableReferenceOffset() == 2)
        {  
          theGetPropensityMethodPtr = RealMethodProxy::
            create<&SpatiocyteNextReactionProcess::
            getPropensity_SecondOrder_TwoSubstrates>();
        }
      //One substrate species, second order
      //A + A -> products:
      else
        {
          theGetPropensityMethodPtr = RealMethodProxy::
            create<&SpatiocyteNextReactionProcess::
            getPropensity_SecondOrder_OneSubstrate>();
        }
    }
}

}
