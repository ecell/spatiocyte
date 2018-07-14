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
  if(react())
    {
      reactAdjoins();
      interruptProcessesPost();
      ReactionProcess::fire();
      return;
    }
  requeue();
}

/*
void SpatiocyteNextReactionProcess::reactAdjoins()
{
  if(theAdjoinSubstratesA.size())
    {
      reactAdjoins(moleculeA, NULL, theAdjoinSubstratesA,
                   theAdjoinProductsA);
    }
}

void SpatiocyteNextReactionProcess::reactAdjoins(Voxel* source,
                                      Voxel* excluded,
                                      std::vector<unsigned>& adjoinSubstrates,
                                      std::vector<unsigned>& adjoinProducts)
{
  for(unsigned i(0); i != source->diffuseSize; ++i)
    { 
      Voxel* mol(&(*theLattice)[source->adjoiningCoords[i]]);
      if(mol != excluded)
        {
          Species* substrate(theSpecies[getID(mol)]);
          std::vector<unsigned>::iterator iterator(std::find(
            adjoinSubstrates.begin(), adjoinSubstrates.end(),
            substrate->getID())); 
          if(iterator != adjoinSubstrates.end())
            {
              Species* product(theSpecies[adjoinProducts[iterator-
                               adjoinSubstrates.begin()]]);

              unsigned index(substrate->getIndex(mol));
              product->addMolecule(mol, substrate->getTag(index));
              substrate->softRemoveMolecule(index);
            }
        }
    }
}
*/


void SpatiocyteNextReactionProcess::updateSubstrateSize()
{
  if(A && A->getIsReactiveVacant())
    {
      A->updateMoleculeSize();
    }
  if(B && B->getIsReactiveVacant())
    {
      B->updateMoleculeSize();
    }
}

void SpatiocyteNextReactionProcess::updateSubstrateMolecules()
{
  if(A && A->getIsReactiveVacant())
    {
      A->updateMolecules();
    }
  if(B && B->getIsReactiveVacant())
    {
      B->updateMolecules();
    }
}

bool SpatiocyteNextReactionProcess::react()
{
  updateSubstrateMolecules();
  if(theOrder == 0)
    {
      if(C)
        { 
          if(C->getIsMultiscale())
            {
              moleculeC = newMultiC();
            }
          else
            {
              moleculeC = C->getRandomPopulatableVoxel(SearchVacant);
            }
          if(moleculeC == NULL)
            {
              return false;
            }
          C->addMolecule(moleculeC);
        }
      else if(variableC)
        {
          variableC->addValue(coefficientC);
        }
    }
  else if(theOrder == 1)
    { 
      //nonHD_A -> nonHD_C + nonHD_D + nonHD_E + nonHD_F + nonHD_H:
      if(BindingSite != -1 && A && C && D && E && F && H)
        {
          return reactACDorFbind(A, C, D, E, F, H);
        }
      //nonHD_A -> nonHD_C + nonHD_D + nonHD_F:
      else if(A && C && D && F)
        {
          return reactACDF(A, C, D, F);
        }
      //nonHD_A -> nonHD_C:
      else if(A && C && !D && !variableD)
        {
          if(BindingSite == -1)
            {
              if(theDeoligomerIndex)
                {
                  return reactDeoligomerize(A, C);
                }
              else if(isMultiAC)
                {
                  return reactMultiAC();
                }
              else
                {
                  return reactAC(A, C);
                }
            }
          else
            {
              return reactACbind(A, C);
            }
        }
      //nonHD_A -> HD_C + HD_D:
      else if(A && variableC && variableD)
        {
          moleculeA = A->getRandomMolecule();
          interruptProcessesPre();
          A->removeMolecule(moleculeA);
          variableC->addValue(coefficientC);
          variableD->addValue(coefficientD);
        }
      //nonHD_A -> HD_C:
      else if(A && variableC && !D && !variableD)
        {
          moleculeA = A->getRandomMolecule();
          interruptProcessesPre();
          A->removeMolecule(moleculeA);
          variableC->addValue(coefficientC);
        }
      //nonHD_A -> nonHD_C + HD_D + HD_F:
      //nonHD_A -> HD_C + nonHD_D + HD_F:
      //nonHD_A -> HD_C + HD_D + nonHD_F:
      else if(A && ((C && variableD && variableF) ||
                    (variableC && D && variableF) ||
                    (variableC && variableD && F)))
        {
          Variable* HD_p1(variableD);
          Variable* HD_p2(variableF);
          int coefficient1(coefficientD);
          int coefficient2(coefficientF);
          Species* nonHD_p(C);
          if(D)
            {
              HD_p1 = variableC;
              coefficient1 = coefficientC;
              nonHD_p = D;
            }
          else if(F)
            {
              HD_p2 = variableC;
              coefficient2 = coefficientC;
              nonHD_p = F;
            }
          if(reactAC(A, nonHD_p))
             {
               HD_p1->addValue(coefficient1);
               HD_p2->addValue(coefficient2);
             }
          else
            {
              return false;
            }
        }
      //nonHD_A -> HD_C + nonHD_D + nonHD_F:
      //nonHD_A -> nonHD_C + HD_D + nonHD_F:
      //nonHD_A -> nonHD_C + nonHD_D + HD_F:
      else if(A && ((variableC && D && F) ||
                    (C && variableD && F) ||
                    (C && D && variableF)))
        {
          Variable* HD_p(variableC);
          int coefficient(coefficientC);
          Species* nonHD_p1(D);
          Species* nonHD_p2(F);
          if(variableD)
            {
              HD_p = variableD;
              coefficient = coefficientD;
              nonHD_p1 = C;
            }
          else if(variableF)
            {
              HD_p = variableF;
              coefficient = coefficientF;
              nonHD_p2 = C;
            }
          if(reactACD(A, nonHD_p1, nonHD_p2))
             {
               HD_p->addValue(coefficient);
             }
          else
            {
              return false;
            }
        }
      //nonHD_A -> nonHD_C + nonHD_D:
      else if(A && C && D)
        {
          if(BindingSite == -1)
            {
              return reactACD(A, C, D);
            }
          else
            {
              return reactACDbind(A, C, D);
            }
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
              return false;
            }
        }
      //HD_A -> nonHD_C:
      else if(variableA && C && !D && !variableD)
        {
          moleculeC = reactvAC(variableA, C);
          if(moleculeC == NULL)
            {
              return false;
            }
          else
            {
              variableA->addValue(coefficientA);
              C->addMolecule(moleculeC);
            }
        }
      //HD_A -> nonHD_C + nonHD_D:
      else if(variableA && C && D)
        {
          moleculeC = NULL;
          moleculeD = NULL;
          Comp* compA(theSpatiocyteStepper->system2Comp(
                         variableA->getSuperSystem()));
          //Occupy C in a voxel of compartment C that adjoins compartment A
          //if A is a surface compartment:
          if(compA != C->getComp() && compA->dimension != 3)
            {
              moleculeC = C->getRandomAdjoiningCompVoxel(compA, SearchVacant);
              if(moleculeC)
                {
                  moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC,
                                                         SearchVacant);
                }
            }
          else if(compA != D->getComp() && compA->dimension != 3)
            {
              moleculeD = D->getRandomAdjoiningCompVoxel(compA, SearchVacant);
              if(moleculeD)
                {
                  moleculeC = C->getRandomAdjoiningVoxel(moleculeD, moleculeD,
                                                         SearchVacant);
                }
            }
          else
            {
              moleculeC = C->getRandomPopulatableVoxel(SearchVacant);
              if(moleculeC)
                {
                  moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC,
                                                         SearchVacant);
                }
            }
          if(moleculeC == NULL || moleculeD == NULL)
            {
              return false;
            }
          variableA->addValue(coefficientA);
          C->addMolecule(moleculeC);
          D->addMolecule(moleculeD);
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
          if(moleculeP == NULL)
            {
              return false;
            }
          variableA->addValue(coefficientA);
          nonHD_p->addMolecule(moleculeP);
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
              if(moleculeC == NULL)
                {
                  return false;
                }
              variableA->addValue(coefficientA);
              variableB->addValue(coefficientB);
              C->addMolecule(moleculeC);
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
              if(moleculeP == NULL)
                {
                  return false;
                }
              variableA->addValue(coefficientA);
              variableB->addValue(coefficientB);
              nonHD_p->addMolecule(moleculeP);
              HD_p->addValue(coefficient);
            }
          //HD + HD -> nonHD + nonHD: 
          else if(C && D)
            {
              moleculeC = reactvAvBC(C);
              moleculeD = NULL;
              if(moleculeC == NULL)
                {
                  moleculeD = reactvAvBC(D);
                  if(moleculeD)
                    {
                      moleculeC = C->getRandomAdjoiningVoxel(moleculeD,
                                                     moleculeD, SearchVacant);
                    }
                }
              else
                { 
                  moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC,
                                                         SearchVacant);
                }
              if(moleculeC == NULL || moleculeD == NULL)
                {
                  return false;
                }
              variableA->addValue(coefficientA);
              variableB->addValue(coefficientB);
              C->addMolecule(moleculeC);
              D->addMolecule(moleculeD);
            }
        }
      //Homodimerization 2HD -> product(s)
      else if(variableA && !variableB && !A && !B)
        {
          //HD + HD -> HD: 
          if(variableC && !variableD && !D)
            {
              variableA->addValue(coefficientA);
              variableC->addValue(coefficientC);
            }
          //HD + HD -> nonHD: 
          else if(C && !variableD && !D)
            { 
              moleculeC = reactvAvBC(C);
              if(moleculeC == NULL)
                {
                  return false;
                }
              variableA->addValue(coefficientA);
              C->addMolecule(moleculeC);
            }
          //HD + HD -> HD + HD: 
          else if(variableC && variableD)
            {
              variableA->addValue(coefficientA);
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
              if(moleculeP == NULL)
                {
                  return false;
                }
              variableA->addValue(coefficientA);
              nonHD_p->addMolecule(moleculeP);
              HD_p->addValue(coefficient);
            }
          //HD + HD -> nonHD + nonHD: 
          else if(C && D)
            {
              moleculeC = reactvAvBC(C);
              moleculeD = NULL;
              if(moleculeC == NULL)
                {
                  moleculeD = reactvAvBC(D);
                  if(moleculeD)
                    {
                      moleculeC = C->getRandomAdjoiningVoxel(moleculeD,
                                                     moleculeD, SearchVacant);
                    }
                }
              else
                { 
                  moleculeD = D->getRandomAdjoiningVoxel(moleculeC, moleculeC,
                                                         SearchVacant);
                }
              if(moleculeC == NULL || moleculeD == NULL)
                {
                  return false;
                }
              variableA->addValue(coefficientA);
              C->addMolecule(moleculeC);
              D->addMolecule(moleculeD);
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
                  return false;
                }
              HD->addValue(coefficient);
            }
          //nonHD + HD -> nonHD:
          //HD + nonHD -> nonHD:
          else if(C && !D && !variableD)
            {
              if(!reactAC(nonHD, C))
                {
                  return false;
                }
              HD->addValue(coefficient);
            }
          //nonHD + HD -> HD:
          //HD + nonHD -> HD:
          else if(variableC && !D && !variableD)
            {
              moleculeS = nonHD->getRandomMolecule();
              nonHD->removeMolecule(moleculeS);
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
                  return false;
                }
              HD->addValue(coefficient);
              HD_p->addValue(coefficient_p);
            }
        }
      //nonHD + nonHD -> product(s)
      else
        {
          //nonHD + nonHD -> nonHD + nonHD
          if(C && D && !F)
            {
              //always true reaction:
              reactABCD();
              if(variableG)
                {
                  variableG->addValue(coefficientG);
                }
              if(variableF)
                {
                  variableF->addValue(coefficientF);
                }
            }
          //nonHD + nonHD -> nonHD
          else
            {
              //nonHD + nonHD -> nonHD
              //always true reaction:
              reactABC();
            }
        }
    }
  return true;
}


Voxel* SpatiocyteNextReactionProcess::newMultiC()
{
  std::vector<unsigned> inMultiCnts;
  moleculeC = C->getRandomPopulatableMulti(SearchVacant, inMultiCnts);
  if(moleculeC == NULL)
    {
      return NULL;
    }
  if(E)
    {
      unsigned cnt(inMultiCnts[E->getID()]);
      double prob(cnt*theRates[1]);
      if(H)
        {
          cnt += inMultiCnts[H->getID()];
          prob += inMultiCnts[H->getID()]*theRates[2];
        }
      prob += (C->getMultiCoordSize()-cnt)*theRates[0];
      prob /= C->getMultiCoordSize();
      if(theRng->Fixed() > prob)
        {
          return NULL;
        }
    }
  return moleculeC;
}

//multiNonHD -> nonHD
bool SpatiocyteNextReactionProcess::reactMultiAC()
{
  const unsigned indexA(theNextIndex);
  moleculeA = A->getMolecule(indexA);
  moleculeC = NULL;
  if(A->getVacantID() == C->getVacantID() || A->getID() == C->getVacantID())
    {
      moleculeC = moleculeA;
    }
  else
    {
      moleculeC = C->getRandomAdjoiningVoxel(moleculeA, SearchVacant);
      if(moleculeC == NULL)
        {
          return false;
        }
    }
  interruptProcessesPre();
  if(A->getIsOnMultiscale() && C->getIsOnMultiscale())
    {
      C->addMoleculeInMulti(moleculeC, A->getTag(indexA)[0].multiIdx);
      A->softRemoveMolecule(indexA);
    }
  else
    {
      std::vector<Tag> tagA(A->getTag(indexA));
      A->removeMolecule(indexA);
      C->addMolecule(moleculeC, tagA);
    }
  return true;
}


//nonHD.nonHD -> nonHD
//Both A and B are immobile nonHD
void SpatiocyteNextReactionProcess::reactABC()
{
  const unsigned rand(theRng->Integer(theCoordsA.size()));
  moleculeA = &(*theLattice)[theCoordsA[rand]];
  moleculeB = A->getRandomDiffuseVoxel(moleculeA, B, 1);
  interruptProcessesPre();
  if(A != C)
    {
      unsigned indexA(A->getIndex(moleculeA));
      std::vector<Tag> tagA(A->getTag(indexA));
      A->removeMolecule(indexA);
      C->addMolecule(moleculeA, tagA);
    }
  unsigned indexB(B->getIndex(moleculeB));
  B->removeMolecule(indexB);
}

//nonHD.nonHD -> nonHD + nonHD
//Both A and B are immobile nonHD
void SpatiocyteNextReactionProcess::reactABCD()
{
  const unsigned rand(theRng->Integer(theCoordsA.size()));
  moleculeA = &(*theLattice)[theCoordsA[rand]];
  moleculeB = A->getRandomDiffuseVoxel(moleculeA, B, 1);
  interruptProcessesPre();
  if(A != C)
    {
      unsigned indexA(A->getIndex(moleculeA));
      std::vector<Tag> tagA(A->getTag(indexA));
      A->removeMolecule(indexA);
      C->addMolecule(moleculeA, tagA);
    }
  if(B != D)
    { 
      unsigned indexB(B->getIndex(moleculeB));
      std::vector<Tag> tagB(B->getTag(indexB));
      B->removeMolecule(indexB);
      D->addMolecule(moleculeB, tagB);
    }
}

//nonHD -> nonHD + nonHD + nonHD
bool SpatiocyteNextReactionProcess::reactACDF(Species* a, Species* c,
                                              Species* d, Species* f)
{
  unsigned indexA(a->getRandomValidIndex());
  if(indexA == a->size())
    {
      return false;
    }
  moleculeA = a->getMolecule(indexA);
  //This is needed when a is species B, used by interruptedPre:
  moleculeB = moleculeA;
  moleculeC = NULL;
  moleculeD = NULL;
  moleculeF = NULL;
  if(a->getVacantID() == c->getVacantID() || a->getID() == c->getVacantID())
    {
      moleculeC = moleculeA;
      moleculeD = d->getRandomAdjoiningVoxel(moleculeC, SearchVacant);
      if(moleculeD == NULL)
        {
          return false;
        }
      moleculeF = f->getRandomAdjoiningVoxel(moleculeC, moleculeD,
                                             SearchVacant);
      if(moleculeF == NULL)
        {
          return false;
        }
    }
  else if(a->getVacantID() == d->getVacantID() ||
          a->getID() == d->getVacantID())
    {
      moleculeD = moleculeA;
      moleculeC = c->getRandomAdjoiningVoxel(moleculeD, SearchVacant);
      if(moleculeC == NULL)
        {
          return false;
        }
      moleculeF = f->getRandomAdjoiningVoxel(moleculeC, moleculeD,
                                             SearchVacant);
      if(moleculeF == NULL)
        {
          return false;
        }
    }
  else if(a->getVacantID() == f->getVacantID() ||
          a->getID() == f->getVacantID())
    {
      moleculeF = moleculeA;
      moleculeC = c->getRandomAdjoiningVoxel(moleculeF, SearchVacant);
      if(moleculeC == NULL)
        {
          return false;
        }
      moleculeD = d->getRandomAdjoiningVoxel(moleculeC, moleculeF,
                                             SearchVacant);
      if(moleculeD == NULL)
        {
          return false;
        }
    }
  else
    {
      moleculeC = c->getRandomAdjoiningVoxel(moleculeA, SearchVacant);
      if(moleculeC == NULL)
        {
          //Only proceed if we can find an adjoining vacant voxel
          //of nonND which can be occupied by C:
          return false;
        }
      moleculeD = d->getRandomAdjoiningVoxel(moleculeC, SearchVacant);
      if(moleculeD == NULL)
        {
          return false;
        }
      moleculeF = f->getRandomAdjoiningVoxel(moleculeC, moleculeD,
                                             SearchVacant);
      if(moleculeF == NULL)
        {
          return false;
        }
    }
  interruptProcessesPre();
  std::vector<Tag> tagA(a->getTag(indexA));
  a->removeMolecule(indexA);
  c->addMolecule(moleculeC, tagA);
  d->addMolecule(moleculeD);
  f->addMolecule(moleculeF);
  return true;
}

//nonHD -> nonHD + nonHD
bool SpatiocyteNextReactionProcess::reactACD(Species* a, Species* c, Species* d)
{
  unsigned indexA(a->getRandomValidIndex());
  if(indexA == a->size())
    {
      return false;
    }
  moleculeA = a->getMolecule(indexA);
  //This is needed when a is species B, used by interruptedPre:
  moleculeB = moleculeA;
  moleculeC = NULL;
  moleculeD = NULL;
  if(a->getVacantID() == c->getVacantID() || a->getID() == c->getVacantID())
    {
      moleculeC = moleculeA;
      moleculeD = d->getRandomAdjoiningVoxel(A, C, moleculeA, SearchVacant);
      if(moleculeD == NULL)
        {
          return false;
        }
    }
  else if(a->getVacantID() == d->getVacantID() ||
          a->getID() == d->getVacantID())
    {
      moleculeD = moleculeA;
      moleculeC = c->getRandomAdjoiningVoxel(A, D, moleculeA, SearchVacant);
      if(moleculeC == NULL)
        {
          return false;
        }
    }
  else
    {
      moleculeC = c->getRandomAdjoiningVoxel(moleculeA, SearchVacant);
      if(moleculeC == NULL)
        {
          //Only proceed if we can find an adjoining vacant voxel
          //of nonND which can be occupied by C:
          return false;
        }
      moleculeD = d->getRandomAdjoiningVoxel(moleculeC, SearchVacant);
      if(moleculeD == NULL)
        {
          return false;
        }
    }
  interruptProcessesPre();
  std::vector<Tag> tagA(a->getTag(indexA));
  std::vector<Tag> tagC;
  std::vector<Tag> tagD;
  const unsigned oligomerSizeC(c->getOligomerSize());
  const unsigned oligomerSizeD(d->getOligomerSize());
  for(unsigned i(0); i < oligomerSizeC && i < tagA.size(); ++i)
    {
      tagC.push_back(tagA[i]);
    }
  for(unsigned i(oligomerSizeC); i < oligomerSizeC+oligomerSizeD &&
      i < tagA.size(); ++i)
    {
      tagD.push_back(tagA[i]);
    }
  a->removeMolecule(indexA);
  c->addMolecule(moleculeC, tagC);
  d->addMolecule(moleculeD, tagD);
  return true;
}

//nonHD (+ E) -> nonHD
bool SpatiocyteNextReactionProcess::reactDeoligomerize(Species* a, Species* c)
{
  if(!A->getBoundCnt(theDeoligomerIndex))
    {
      return false;
    }
  unsigned indexA(a->getRandomOligomerIndex(theDeoligomerIndex));
  moleculeA = a->getMolecule(indexA);
  moleculeC = NULL;
  if(a->getVacantID() == c->getVacantID() || a->getID() == c->getVacantID())
    {
      moleculeC = moleculeA;
    }
  else
    {
      moleculeC = c->getRandomAdjoiningVoxel(moleculeA, SearchVacant);
      if(moleculeC == NULL)
        {
          //Only proceed if we can find an adjoining vacant voxel
          //of nonND which can be occupied by C:
          return false;
        }
    }
  interruptProcessesPre();
  if(a->getIsOnMultiscale() && c->getIsOnMultiscale())
    {
      c->addMoleculeInMulti(moleculeC, a->getTag(indexA)[0].multiIdx);
      a->softRemoveMolecule(indexA);
    }
  else
    {
      std::vector<Tag> tagA(a->getTag(indexA));
      a->removeMolecule(indexA);
      c->addMolecule(moleculeC, tagA);
    }
  return true;
}

//nonHD (+ E) -> nonHD
bool SpatiocyteNextReactionProcess::reactAC(Species* a, Species* c)
{
  unsigned indexA(a->getRandomValidIndex());
  if(indexA == a->size())
    {
      return false;
    }
  moleculeA = a->getMolecule(indexA);
  //This is needed when a is species B, used by interruptedPre:
  moleculeB = moleculeA;
  if(ImplicitUnbind && 
     E->getRandomAdjoiningVoxel(moleculeA, E, SearchVacant) == NULL)
    {
      return false;
    }
  moleculeC = NULL;
  if(a->getVacantID() == c->getVacantID() || a->getID() == c->getVacantID())
    {
      moleculeC = moleculeA;
    }
  else
    {
      moleculeC = c->getRandomAdjoiningVoxel(moleculeA, SearchVacant);
      if(moleculeC == NULL)
        {
          //Only proceed if we can find an adjoining vacant voxel
          //of nonND which can be occupied by C:
          return false;
        }
    }
  interruptProcessesPre();
  if(a->getIsOnMultiscale() && c->getIsOnMultiscale())
    {
      c->addMoleculeInMulti(moleculeC, a->getTag(indexA)[0].multiIdx);
      a->softRemoveMolecule(indexA);
    }
  else
    {
      std::vector<Tag> tagA(a->getTag(indexA));
      a->removeMolecule(indexA);
      c->addMolecule(moleculeC, tagA);
    }
  removeMoleculeE();
  return true;
}

//zero-coefficient E
//we remove a molecule E at random location in the compartment to allow
//rebinding effect of the dissociated nonHD molecule while maintaining (fixing)
//the concentration of a substrate species even after the reaction:
void SpatiocyteNextReactionProcess::removeMoleculeE()
{
  if(!E || ImplicitUnbind)
    {
      return;
    }
  moleculeE = E->getRandomMolecule();
  if(moleculeE == NULL)
    {
      cout << getFullID().asString() << " unable to remove molecule E" <<
        std::endl;
    }
  E->removeMolecule(moleculeE);
}

//A (+Vacant[BindingSite]) -> nonHD[BindingSite]
bool SpatiocyteNextReactionProcess::reactACbind(Species* a, Species* c)
{
  unsigned indexA(a->getRandomValidIndex());
  if(indexA == a->size())
    {
      return false;
    }
  moleculeA = a->getMolecule(indexA);
  moleculeC = c->getBindingSiteAdjoiningVoxel(moleculeA, BindingSite);
  if(moleculeC == NULL)
    {
      //Only proceed if we can find an adjoining vacant voxel
      //of nonND which can be occupied by C:
      return false;
    }
  interruptProcessesPre();
  std::vector<Tag> tagA(a->getTag(indexA));
  a->removeMolecule(indexA);
  c->addMolecule(moleculeC, tagA);
  return true;
}

//molA <- SpeciesC
//
bool SpatiocyteNextReactionProcess::reactACDorFbind(Species* a, Species* c,
                                                    Species* d, Species* e,
                                                    Species* f, Species* h)
{
  unsigned indexA(a->getRandomValidIndex());
  if(indexA == a->size())
    {
      return false;
    }
  //moleculeA will be replaced with species C:
  //we must define moleculeA since it will be used in reactAdjoins:
  moleculeA = a->getMolecule(indexA);
  moleculeC = moleculeA;
  //Get valid molecule D or F to create a product at the binding site:
  moleculeD = d->getBindingSiteAdjoiningVoxel(moleculeC, BindingSite, e);
  if(moleculeD == NULL)
    {
      moleculeF = f->getBindingSiteAdjoiningVoxel(moleculeC, BindingSite, h);
      if(moleculeF == NULL)
        {
          return false;
        }
      //Need to specify moleculeB since it will be used in reactAdjoins:
      if(theAdjoinSubstratesB.size())
        {
          moleculeB = moleculeF;
        }
      interruptProcessesPre();
      std::vector<Tag> tagA(a->getTag(indexA));
      a->removeMolecule(indexA);
      c->addMolecule(moleculeC);
      h->softRemoveMolecule(moleculeF);
      f->addMolecule(moleculeF, tagA);
      return true;
    }
  //Need to specify moleculeB since it will be used in reactAdjoins:
  if(theAdjoinSubstratesB.size())
    {
      moleculeB = moleculeD;
    }
  interruptProcessesPre();
  std::vector<Tag> tagA(a->getTag(indexA));
  a->removeMolecule(indexA);
  c->addMolecule(moleculeC);
  e->softRemoveMolecule(moleculeD);
  d->addMolecule(moleculeD, tagA);
  return true;
}


//A (+Vacant[BindingSite] || +D[BindingSite]) -> 
//[molA <- D_trailA] + (Vacant[BindingSite] || D[BindingSite] <- C)
bool SpatiocyteNextReactionProcess::reactACDbind(Species* a, Species* c,
                                                 Species* d)
{
  unsigned indexA(a->getRandomValidIndex());
  if(indexA == a->size())
    {
      return false;
    }
  moleculeA = a->getMolecule(indexA);
  //Look for Vacant[BindingSite]:
  moleculeC = c->getBindingSiteAdjoiningVoxel(moleculeA, BindingSite);
  if(moleculeC == NULL)
    {
      //Look for D[BindingSite]:
      moleculeC = c->getBindingSiteAdjoiningVoxel(moleculeA, BindingSite, d);
      if(moleculeC == NULL)
        {
          return false;
        }
      interruptProcessesPre();
      std::vector<Tag> tagA(a->getTag(indexA));
      //molA <- D
      if(a->isTrailSpecies(d))
        {
          a->softRemoveMolecule(indexA);
          d->addMolecule(moleculeA);
        }
      else
        {
          a->removeMolecule(indexA);
        }
      //D[BindingSite] <- C
      d->softRemoveMolecule(moleculeC);
      c->addMolecule(moleculeC, tagA);
      return true;
    }
  interruptProcessesPre();
  std::vector<Tag> tagA(a->getTag(indexA));
  a->removeMolecule(indexA);
  if(a->isTrailSpecies(d))
    {
      //molA <- D
      d->addMolecule(moleculeA);
    }
  //Vacant[BindingSite] <- C
  c->addMolecule(moleculeC, tagA);
  return true;
}

//HD -> nonHD
Voxel* SpatiocyteNextReactionProcess::reactvAC(Variable* vA, Species* c)
{
  moleculeC = NULL;
  Comp* compA(theSpatiocyteStepper->system2Comp(vA->getSuperSystem()));
  //Occupy C in a voxel of compartment C that adjoins compartment A
  //if A is a surface compartment:
  if(compA != c->getComp() && compA->dimension != 3)
    {
      moleculeC = c->getRandomAdjoiningCompVoxel(compA, SearchVacant);
    }
  else
    {
      moleculeC = c->getRandomPopulatableVoxel(SearchVacant);
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

Voxel* SpatiocyteNextReactionProcess::reactvAvBC(Species* c)
{
  moleculeC = NULL;
  Comp* aComp2D(getComp2D(c));
  if(aComp2D)
    {
      moleculeC = C->getRandomAdjoiningCompVoxel(aComp2D, SearchVacant);
    }
  else
    {
      moleculeC = C->getRandomPopulatableVoxel(SearchVacant);
    }
  return moleculeC;
}


double SpatiocyteNextReactionProcess::getPropensityZerothOrder() 
{
  return p;
}

double SpatiocyteNextReactionProcess::getPropensityFirstOrder() 
{
  double sizeA(substrateA->getValue());
  if(sizeA < -coefficientA)
    {
      sizeA = 0;
    }
  return p*sizeA;
}

double SpatiocyteNextReactionProcess::getPropensityFirstOrderMultiAC() 
{
  double sizeA(substrateA->getValue());
  if(sizeA < -coefficientA)
    {
      sizeA = 0;
    }
  return p*sizeA;
}

double SpatiocyteNextReactionProcess::getPropensityFirstOrderReactAB() 
{
  double sizeA(theCoordsA.size());
  if(sizeA < -coefficientA)
    {
      sizeA = 0;
    }
  return p*sizeA;
}

double SpatiocyteNextReactionProcess::getPropensityFirstOrderDeoligomerize() 
{
  return p*A->getBoundCnt(theDeoligomerIndex);
}

double SpatiocyteNextReactionProcess::getPropensitySecondOrderReactABvG() 
{
  double sizeA(theCoordsA.size());
  if(sizeA < -coefficientA)
    {
      sizeA = 0;
    }
  else
    {
      sizeA = pow(sizeA, -coefficientA);
    }
  double sizeG(variableG->getValue());
  if(sizeG < -coefficientG)
    {
      sizeG = 0;
    }
  else
    {
      sizeG = pow(sizeG, -coefficientG);
    }
  return p*sizeA*sizeG;
}

//Need to solve homodimerization reaction of two substrate species (Size-1):
double SpatiocyteNextReactionProcess::getPropensitySecondOrderHetero() 
{
  //Don't do the instructions below, if it doesn't work without them,
  //use interrupedEndDiffusion or interrupedAddMolecule or 
  //interruptedRemoveMolecule
  /*
  if(A)
    {
      A->updateMoleculeSize();
    }
  else if(B)
    {
      B->updateMoleculeSize();
    }
    */
  double sizeA(substrateA->getValue());
  double sizeB(substrateB->getValue());
  //Required for HD species when substrate coefficient is < -1
  if(sizeA < -coefficientA)
    {
      sizeA = 0;
    }
  else
    {
      sizeA = pow(sizeA, -coefficientA);
    }
  if(sizeB < -coefficientB)
    {
      sizeB = 0;
    }
  else
    {
      sizeB = pow(sizeB, -coefficientB);
    }
  return p*sizeA*sizeB;
}

double SpatiocyteNextReactionProcess::getPropensitySecondOrderHomo() 
{
  //Don't do the instructions below, if it doesn't work without them,
  //use interrupedEndDiffusion or interrupedAddMolecule or 
  //interruptedRemoveMolecule
  /*
  if(A)
    {
      A->updateMoleculeSize();
    }
    */
  double sizeA(substrateA->getValue());
  //Include coefficientB because A + A -> product, with coefficientA = -1
  //and coefficientB = -1:
  if(sizeA < -coefficientA-coefficientB)
    {
      sizeA = 1;
    }
  //There must be two or more molecules:
  return p*sizeA*(sizeA-1);
}


void SpatiocyteNextReactionProcess::preinitialize()
{
  ReactionProcess::preinitialize();
  if(Deoligomerize)
    {
      setDeoligomerIndex(Deoligomerize);
      for(unsigned i(0); i != Deoligomerize-1; ++i)
        {
          const unsigned aDeoligomerIndex(i+1);
          Process* aProcess(SpatiocyteStepper::createProcess(
             getPropertyInterface().getClassName(),
             String(getFullID().getID()+int2str(aDeoligomerIndex)),
             getStepper()->getModel()->getSystem(getFullID().getSystemPath()),
             getStepper()->getModel()));
          aProcess->setStepper(getStepper());
          aProcess->setPriority(getPriority());
          SpatiocyteNextReactionProcess* aSNRP(
                     dynamic_cast<SpatiocyteNextReactionProcess*>(aProcess));
          if(theRates.size())
            {
              aSNRP->setk(theRates[i]);
            }
          else
            {
              aSNRP->setk(k/aDeoligomerIndex);
            }
          aSNRP->setSearchVacant(SearchVacant);
          aSNRP->setDeoligomerIndex(aDeoligomerIndex);
          aSNRP->setImplicitUnbind(ImplicitUnbind);
          aSNRP->setVariableReferences(theVariableReferenceVector);
          aProcess->preinitialize();
        }
      if(theRates.size())
        {
          setk(theRates[Deoligomerize-1]);
        }
      else
        {
          setk(k/Deoligomerize);
        }
    }
}

void SpatiocyteNextReactionProcess::setDeoligomerIndex(const unsigned anIndex)
{
  theDeoligomerIndex = anIndex;
}

void SpatiocyteNextReactionProcess::setVariableReferences(const 
                          VariableReferenceVector& aVariableReferenceVector)
{
  theVariableReferenceVector = aVariableReferenceVector;
}

void SpatiocyteNextReactionProcess::initializeSecond()
{
  ReactionProcess::initializeSecond();
  if(A)
    {
      substrateA = A->getVariable();
    }
  else if(variableA)
    {
      substrateA = variableA;
    }
  if(B)
    {
      substrateB = B->getVariable();
    }
  else if(variableB)
    {
      substrateB = variableB;
    }
  //if second order, with A and B substrates, both of them
  //must be immobile, or A must be multiscale:
  if(A && B)
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
  if(A && C)
    {
      if(Deoligomerize)
        {
          A->setIsDeoligomerize(C, Deoligomerize);
        }
      //nonHD -> nonHD + nonHD
      //used by multiscale dissociation reaction
      if(!B && D)
        {
          C->setProductPair(D);
          D->setProductPair(C);
        }
    }
  setPropensityMethod();
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

void SpatiocyteNextReactionProcess::initializeBeforePopulate()
{
  if(p != -1)
    {
      return;
    }
  Comp* compA(NULL);
  Comp* compB(NULL);
  Comp* compC(NULL);
  Comp* compD(NULL);
  int secondCoefficient(coefficientB);
  if(A)
    {
      compA = A->getComp();
      compB = compA;
    }
  else if(variableA)
    {
      compA = theSpatiocyteStepper->system2Comp(variableA->getSuperSystem());
      compB = compA;
    }
  if(B)
    {
      compB = B->getComp();
    }
  else if(variableB)
    {
      compB = theSpatiocyteStepper->system2Comp(variableB->getSuperSystem());
    }
  if(variableG)
    {
      compB = theSpatiocyteStepper->system2Comp(variableG->getSuperSystem());
      secondCoefficient = coefficientG;
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
  //Now let's determine the velocity, p (which has the unit 1/s)
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
  else if(theOrder == 1 || (isReactAB && !variableG)) 
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
          THROW_EXCEPTION(ValueError, String(
                            getPropertyInterface().getClassName()) +
                            "[" + getFullID().asString() + 
                            "]: In SpatiocyteNextReactionProcess if two " +
                            "products don't belong to the same compartment, " +
                            "the substrates must also belong to different " +
                            " compartments.");
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
          //unit of k is in (m^3)^(totalCoefficient-1)/s
          //we need to convert k to p which has the unit 1/s
          int totalCoefficient(coefficientA+secondCoefficient);
          p = k/(pow(aVolume, fabs(totalCoefficient)-1));
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
          //unit of k is in (m^2)^(totalCoefficient-1)/s
          //we need to convert k to p which has the unit 1/s
          int totalCoefficient(coefficientA+secondCoefficient);
          p = k/(pow(anArea, fabs(totalCoefficient)-1));
          pFormula << "[k/anArea:" << k << "/" << anArea << "]";
        }
      //If volume + filament <= k(volume)(filament) or
      //   volume + filament <= k(filament)(volume) or
      //   filament + volume <= k(volume)(filament) or
      //   filament + volume <= k(filament)(volume)
      else if((compD && (
        (compC->dimension == 3 && compD->dimension == 1 &&
         compA->dimension == 3 && compB->dimension == 1) ||
        (compC->dimension == 3 && compD->dimension == 1 &&
         compA->dimension == 1 && compB->dimension == 3) ||
        (compC->dimension == 1 && compD->dimension == 3 &&
         compA->dimension == 3 && compB->dimension == 1) ||
        (compC->dimension == 1 && compD->dimension == 3 &&
         compA->dimension == 1 && compB->dimension == 3))) ||
         ((compC->dimension == 3 && compA->dimension == 3
          && compB->dimension == 3) ||
         (compC->dimension == 1 && compA->dimension == 3 
          && compB->dimension == 1) ||
         (compC->dimension == 1 && compA->dimension == 1 
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
          //unit of k is in (m^3)^(totalCoefficient-1)/s
          //we need to convert k to p which has the unit 1/s
          int totalCoefficient(coefficientA+secondCoefficient);
          p = k/(pow(aVolume, fabs(totalCoefficient)-1));
          pFormula << "[k/aVolume:" << k << "/" << aVolume << "]";
        }
      else
        {
          THROW_EXCEPTION(ValueError, String(
                            getPropertyInterface().getClassName()) +
                            "[" + getFullID().asString() + 
                            "]: In SpatiocyteNextReactionProcess: " +
                            "unexpected compartment dimension combinations " +
                            "in pFormula calculation.");
        }
    }
}

void SpatiocyteNextReactionProcess::initializeFourth()
{
  ReactionProcess::initializeFourth();
  /*
  //Done after populating molecules:
  if(isReactAB)
    {
      for(unsigned i(0); i != A->size(); ++i)
        {
          const unsigned coordA(A->getCoord(i));
          const Voxel* molA(A->getMolecule(i));
          for(unsigned j(0); j != molA->diffuseSize; ++j)
            {
              const unsigned adjCoord(molA->adjoiningCoords[j]);
              if(getID((*theLattice)[adjCoord]) == B->getID())
                {
                  theCoordsA.push_back(coordA);
                }
            }
        }
    }
    */
}

/*
getNewIntervalMultiAC()
{
  const double aCurrentTime(getStepper()->getCurrentTime());
  double anInterval(libecs::INF);
  for(unsigned i(0); i != theNextTimes.size(); ++i)
    {
      const double indexInterval(theNextTimes[i]-aCurrentTime);
      if(anInterval > indexInterval)
        {
          anInterval = indexInterval;
          theNextIndex = i;
        }
    }
  return anInterval;
}
*/

/*
double SpatiocyteNextReactionProcess::getIntervalUnbindMultiAB()
{
  if(!A->size() || !p)
    {
      return libecs::INF;
    }
  nextIndexA = 0;
  double fraction(A->getMultiscaleBoundFraction(nextIndexA,
                                            B->getVacantSpecies()->getID())); 
  double rand(theRng->FixedU());
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
      rand = theRng->FixedU();
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
  A->updateMoleculeSize();
  B->updateMoleculeSize();
  if(A->getIsMultiscale())
    {
      return getIntervalUnbindMultiAB();
    }
  updateMoleculesA();
  const double sizeA(moleculesA.size());
  const double rand(theRng->FixedU());
  const double denom((p*sizeA)*(-log(rand)));
  if(denom)
    {
      return 1.0/denom;
    }
  return libecs::INF;
}
*/


void SpatiocyteNextReactionProcess::printParameters()
{
  String aProcess(String(getPropertyInterface().getClassName()) + 
                                      "[" + getFullID().asString() + "]");
  cout << aProcess << " SearchVacant:" << SearchVacant << std::endl;
  if(A)
    {
      cout << "  " << getIDString(A);
    }
  else if(variableA)
    {
      cout << "  " << getIDString(variableA);
    }
  if(B)
    {
      cout << " + " << getIDString(B);
    }
  else if(variableB)
    {
      cout << " + " << getIDString(variableB);
    }
  if(G)
    {
      cout << " + " << getIDString(G);
    }
  else if(variableG)
    {
      cout << " + " << getIDString(variableG);
    }
  if(!A && !variableA)
    {
      if(C)
        {
          cout << "0 -> " << getIDString(C);
        }
      else if(variableC)
        {
          cout << "0 -> " << getIDString(variableC);
        }
    }
  else
    {
      if(C)
        {
          cout << " -> " << getIDString(C);
        }
      else if(variableC)
        {
          cout << " -> " << getIDString(variableC);
        }
    }
  if(D)
    {
      cout << " + " << getIDString(D);
    }
  else if(variableD)
    {
      cout << " + " << getIDString(variableD);
    }
  if(F)
    {
      cout << " + " << getIDString(F);
    }
  else if(variableF)
    {
      cout << " + " << getIDString(variableF);
    }
  double interval(theTime-getStepper()->getCurrentTime());
  double actualInterval(interval);
  double propensity(thePropensity);
  if(theTime == libecs::INF)
    {
      bool a(false);
      bool b(false);
      bool vA(false);
      bool vB(false);
      if(A && !A->getVariable()->getValue())
        {
          A->getVariable()->addValue(1);
          a = true;
        }
      if(B && !B->getVariable()->getValue())
        {
          B->getVariable()->addValue(1);
          b = true;
        }
      if(variableA && !variableA->getValue())
        {
          variableA->addValue(1);
          vA = true;
        }
      if(variableB && !variableB->getValue())
        {
          variableB->addValue(1);
          vB = true;
        }
      interval = getInitInterval();
      propensity = thePropensity; 
      if(a)
        {
          A->getVariable()->addValue(-1);
        }
      if(b)
        {
          B->getVariable()->addValue(-1);
        }
      if(vA)
        {
          variableA->addValue(-1);
        }
      if(vB)
        {
          variableB->addValue(-1);
        }
      actualInterval = getInitInterval();
    }
  cout << " k:" << k << " p = " << pFormula.str() << " = " << p
    << " predicted nextTime:" << interval << " actual nextTime:" << 
    actualInterval << " propensity:" << propensity;
  if(theDeoligomerIndex)
    {
      cout << " DeoligomerIndex:" << theDeoligomerIndex;
    }
  cout << std::endl;
}

double SpatiocyteNextReactionProcess::getNewInterval()
{
  if(isMultiAC)
    {
      if(theNextTimes.size())
        {
          return theNextTimes[theNextIndex]-getStepper()->getCurrentTime();
        }
      return libecs::INF;
      //return theInterval;
    }
  if(getNewPropensity())
    {
      return -log(theRng->FixedU())/thePropensity;
    }
  return libecs::INF;
}

double SpatiocyteNextReactionProcess::getInterval(double aCurrentTime)
{
  if(isMultiAC)
    {
      if(theNextTimes.size())
        {
          return theNextTimes[theNextIndex]-aCurrentTime;
        }
      return libecs::INF;
    }
  if(theTime == libecs::INF)
    {
      return getNewInterval();
    }
  const double oldPropensity(thePropensity);
  if(getNewPropensity())
    {
      return oldPropensity/thePropensity*(theTime-aCurrentTime);
    }
  return libecs::INF;
}

double SpatiocyteNextReactionProcess::getNewPropensity()
{
  updateSubstrateSize();
  thePropensity = (this->*thePropensityMethod)();
  return thePropensity;
}

double SpatiocyteNextReactionProcess::getPropensity() const
{
  return thePropensity;
}

//Find out if this process is interrupted by the aProcess
//by checking if any of the modified variables of aProcess is a
//substrate of this process:
bool SpatiocyteNextReactionProcess::isDependentOn(const Process* aProcess) const
{
  //Check if any variable with netCoefficient != 0 is a substrate
  //of this process:
  for(VariableReferenceVector::const_iterator
      i(theVariableReferenceVector.begin());
      i != theVariableReferenceVector.end(); ++i)
    {
      if((*i).getCoefficient() < 0)
        {
          if(getVariableNetCoefficient(aProcess, (*i).getVariable()))
            {
              return true;
            }
        }
    }
  return false;
}

bool SpatiocyteNextReactionProcess::isDependentOnPre(const Process* aProcess)
{
  return false;
}

bool SpatiocyteNextReactionProcess::isDependentOnPost(const Process* aProcess)
{
  return false;
}

bool SpatiocyteNextReactionProcess::isDependentOnEndDiffusion(Species* aSpecies)
{
  if(isMultiAC)
    {
      if(A == aSpecies)
        {
          return true;
        }
    }
  else if(isReactAB)
    {
      if(aSpecies->getDiffusionCoefficient() && 
         (aSpecies->getVacantSpecies() == A ||
          aSpecies->getVacantSpecies() == B))
        {
          return true;
        }
    }
  return false;
}

//Cannot put isReactiveVacant here because SpatiocyteSpecies does not
//call this interrupt for vacant species:
bool SpatiocyteNextReactionProcess::isDependentOnRemoveMolecule(Species*
                                                                aSpecies)
{
  if(isMultiAC)
    {
      if(A == aSpecies)
        {
          return true;
        }
    }
  else if(isReactAB)
    {
      if(A == aSpecies || B == aSpecies || 
         aSpecies->getVacantSpecies() == A ||
         aSpecies->getVacantSpecies() == B)
        {
          return true;
        }
    }
  return false;
}

//Cannot put isReactiveVacant here because SpatiocyteSpecies does not
//call this interrupt for vacant species:
bool SpatiocyteNextReactionProcess::isDependentOnAddMolecule(Species* aSpecies)
{
  if(isMultiAC)
    {
      if(A == aSpecies)
        {
          return true;
        }
    }
  else if(isReactAB)
    {
      if(A == aSpecies || B == aSpecies || 
         aSpecies->getVacantSpecies() == A ||
         aSpecies->getVacantSpecies() == B)
        {
          return true;
        }
    }
  return false;
}

void SpatiocyteNextReactionProcess::interruptedEndDiffusion(Species* aSpecies)
{
  if(isMultiAC)
    {
      const double aCurrentTime(getStepper()->getCurrentTime());
      theInterval = libecs::INF;
      std::vector<unsigned> inMultiCnts;
      inMultiCnts.resize(theSpecies.size());
      for(unsigned i(0); i != A->size(); ++i)
        {
          A->getInMultiCnts(i, inMultiCnts);
          unsigned cnt(inMultiCnts[E->getID()]);
          inMultiCnts[E->getID()] = 0;
          double fraction(cnt*theRates[1]);
          if(H)
            {
              const unsigned cntH(inMultiCnts[H->getID()]);
              inMultiCnts[H->getID()] = 0;
              cnt += cntH;
              fraction += cntH*theRates[2];
            }
          fraction += (A->getMultiCoordSize()-cnt)*theRates[0];
          fraction /= A->getMultiCoordSize();
          const double aPropensity(p*fraction);
          const double anInterval(thePropensities[i]/aPropensity*
                                 (theNextTimes[i]-aCurrentTime));
          thePropensities[i] = aPropensity;
          theNextTimes[i] = anInterval+aCurrentTime;
          if(anInterval < theInterval)
            {
              theInterval = anInterval;
              theNextIndex = i;
            }
        }
      theSpatiocyteStepper->addInterruptedProcess(
                                      dynamic_cast<SpatiocyteProcess*>(this));
    }
  else if(isReactAB)
    {
      //do the following 
      //if(aSpecies->getDiffusionCoefficient() && 
      //   (aSpecies->getVacantSpecies() == A ||
      //    aSpecies->getVacantSpecies() == B))
      const unsigned prevSize(theCoordsA.size());
      theCoordsA.resize(0);
      for(unsigned i(0); i != A->size(); ++i)
        {
          const unsigned coordA(A->getCoord(i));
          const Voxel* molA(A->getMolecule(i));
          for(unsigned j(0); j != molA->diffuseSize; ++j)
            {
              const unsigned adjCoord(molA->adjoiningCoords[j]);
              if(getID((*theLattice)[adjCoord]) == B->getID())
                {
                  theCoordsA.push_back(coordA);
                }
            }
        }
      if(prevSize != theCoordsA.size())
        {
          theSpatiocyteStepper->addInterruptedProcess(
                                      dynamic_cast<SpatiocyteProcess*>(this));
        }
    }
}

//Cannot put isReactiveVacant here because SpatiocyteSpecies does not
//call this interrupt for vacant species:
//This method is called after the molecule pointed by the index is added
//by aSpecies.
void SpatiocyteNextReactionProcess::interruptedAddMolecule(Species* aSpecies,
                                                           const unsigned index)
{
  if(isMultiAC)
    {
      const double aCurrentTime(getStepper()->getCurrentTime());
      theInterval = theTime-aCurrentTime;
      std::vector<unsigned> inMultiCnts;
      inMultiCnts.resize(theSpecies.size());
      A->getInMultiCnts(index, inMultiCnts);
      unsigned cnt(inMultiCnts[E->getID()]);
      double fraction(cnt*theRates[1]);
      if(H)
        {
          const unsigned cntH(inMultiCnts[H->getID()]);
          cnt += cntH;
          fraction += cntH*theRates[2];
        }
      fraction += (A->getMultiCoordSize()-cnt)*theRates[0];
      fraction /= A->getMultiCoordSize();
      const double aPropensity(p*fraction);
      const double anInterval(-log(theRng->FixedU())/aPropensity);
      thePropensities.push_back(aPropensity);
      theNextTimes.push_back(anInterval+aCurrentTime);
      if(anInterval < theInterval)
        {
          theInterval = anInterval;
          theNextIndex = index;
          theSpatiocyteStepper->addInterruptedProcess(
                                      dynamic_cast<SpatiocyteProcess*>(this));
        }
    }
  else if(isReactAB)
    {
      const unsigned prevSize(theCoordsA.size());
      if(aSpecies == A)
        {
          addCoordsA(A, B, aSpecies->getCoord(index));
        }
      else if(aSpecies == B)
        { 
          addCoordsA(B, A, aSpecies->getCoord(index));
        }
      //Cannot use "else if" here because B could be the vacant species of A
      //and we need to update it accordingly:
      if(aSpecies->getVacantSpecies() == A)
        {
          removeCoordsA(aSpecies->getCoord(index));
        }
      else if(aSpecies->getVacantSpecies() == B)
        {
          removeAdjCoordsA(aSpecies->getCoord(index));
        }
      if(prevSize != theCoordsA.size())
        {
          theSpatiocyteStepper->addInterruptedProcess(
                                      dynamic_cast<SpatiocyteProcess*>(this));
        }
    }
}

//Cannot put isReactiveVacant here because SpatiocyteSpecies does not
//call this interrupt for vacant species:
//This method is called before the molecule pointed by the index is removed
//by aSpecies.
void SpatiocyteNextReactionProcess::interruptedRemoveMolecule(Species* aSpecies,
                                                           const unsigned index)
{
  if(isMultiAC)
    {
      thePropensities[index] = thePropensities.back();
      thePropensities.pop_back();
      theNextTimes[index] = theNextTimes.back();
      theNextTimes.pop_back();
      if(index == theNextIndex)
        {
          theInterval = libecs::INF;
          const double aCurrentTime(getStepper()->getCurrentTime());
          for(unsigned i(0); i != theNextTimes.size(); ++i)
            {
              if(theNextTimes[i]-aCurrentTime < theInterval)
                {
                  theInterval = theNextTimes[i]-aCurrentTime;
                  theNextIndex = i;
                }
            }
          theSpatiocyteStepper->addInterruptedProcess(
                                      dynamic_cast<SpatiocyteProcess*>(this));
        }
      else if(theNextIndex == theNextTimes.size())
        {
          theNextIndex = index;
        }
    }
  else if(isReactAB)
    {
      const unsigned prevSize(theCoordsA.size());
      if(aSpecies == A)
        {
          removeCoordsA(aSpecies->getCoord(index));
        }
      else if(aSpecies == B)
        {
          removeAdjCoordsA(aSpecies->getCoord(index));
        }
      //Cannot use "else if" here because B could be the vacant species of A
      //and we need to update it accordingly:
      if(aSpecies->getVacantSpecies() == A)
        {
          addCoordsA(A, B, aSpecies->getCoord(index));
        }
      else if(aSpecies->getVacantSpecies() == B)
        { 
          addCoordsA(B, A, aSpecies->getCoord(index));
        }
      if(prevSize != theCoordsA.size())
        {
          theSpatiocyteStepper->addInterruptedProcess(
                                      dynamic_cast<SpatiocyteProcess*>(this));
        }
    }
}

void SpatiocyteNextReactionProcess::interruptedPre(ReactionProcess* aProcess)
{
}

void SpatiocyteNextReactionProcess::interruptedPost(ReactionProcess* aProcess)
{
}

void SpatiocyteNextReactionProcess::removeCoordsA(const unsigned coordA)
{
  unsigned i(0);
  //Use the while loop twice here because there could be multiple instances
  //of coordA in the list theCoordsA
  while(i < theCoordsA.size())
    {
      while(theCoordsA[i] == coordA && i < theCoordsA.size())
        {
          theCoordsA[i] = theCoordsA.back();
          theCoordsA.pop_back();
        }
      ++i;
    }
}

void SpatiocyteNextReactionProcess::removeAdjCoordsA(const unsigned coordB)
{
  const Voxel& molB((*theLattice)[coordB]);
  for(unsigned i(0); i != molB.diffuseSize; ++i)
    {
      const unsigned adjCoord(molB.adjoiningCoords[i]);
      if(getID((*theLattice)[adjCoord]) == A->getID())
        {
          removeSingleCoordsA(adjCoord);
        }
    }
}

void SpatiocyteNextReactionProcess::removeSingleCoordsA(const unsigned coordA)
{
  for(unsigned i(0); i != theCoordsA.size(); ++i)
    {
      if(theCoordsA[i] == coordA)
        {
          theCoordsA[i] = theCoordsA.back();
          theCoordsA.pop_back();
          return;
        }
    }
}

void SpatiocyteNextReactionProcess::addCoordsA(Species* a, Species* b,
                                               const unsigned coord)
{
  const Voxel& mol((*theLattice)[coord]);
  for(unsigned i(0); i != mol.diffuseSize; ++i)
    {
      const unsigned adjCoord(mol.adjoiningCoords[i]);
      if(getID((*theLattice)[adjCoord]) == b->getID())
        {
          if(a == A)
            {
              theCoordsA.push_back(coord);
            }
          else
            {
              theCoordsA.push_back(adjCoord);
            }
        }
    }
}

void SpatiocyteNextReactionProcess::setPropensityMethod()
{
  if(theOrder == 0) // no substrate
    {
      thePropensityMethod = &SpatiocyteNextReactionProcess::
        getPropensityZerothOrder;
    }
  else if(theOrder == 1) // one substrate, first order.
    {
      if(theDeoligomerIndex)
        {
          thePropensityMethod = &SpatiocyteNextReactionProcess::
            getPropensityFirstOrderDeoligomerize;
        }
      else if(A && A->getIsMultiscale() && C && !variableC && E)
        {
          isMultiAC = true;
          thePropensityMethod = &SpatiocyteNextReactionProcess::
            getPropensityFirstOrderMultiAC;
        }
      else
        {
          thePropensityMethod = &SpatiocyteNextReactionProcess::
            getPropensityFirstOrder;
        }
    }
  else
    { 
      //nonHD.nonHD -> products
      if(A && B)
        {
          isReactAB = true;
          if(variableG)
            {
              thePropensityMethod = &SpatiocyteNextReactionProcess::
                getPropensitySecondOrderReactABvG;
            }
          else
            {
              thePropensityMethod = &SpatiocyteNextReactionProcess::
                getPropensityFirstOrderReactAB;
            }
        }
      //Two unique substrate species:
      //A + B -> products:
      else if(getZeroVariableReferenceOffset() >= 2 && variableA != variableB)
        {
          thePropensityMethod = &SpatiocyteNextReactionProcess::
            getPropensitySecondOrderHetero;
        }
      //One substrate species, second order
      //A + A -> products:
      else
        {
          thePropensityMethod = &SpatiocyteNextReactionProcess::
            getPropensitySecondOrderHomo;
        }
    }
}

void SpatiocyteNextReactionProcess::checkExternStepperInterrupted()
{
  Model::StepperMap aStepperMap(getModel()->getStepperMap());  
  for(Model::StepperMap::const_iterator i(aStepperMap.begin());
      i != aStepperMap.end(); ++i )
    {
      if(i->second != getStepper())
        {
          std::vector<Process*> aProcessVector(i->second->getProcessVector());
          for(std::vector<Process*>::const_iterator j(aProcessVector.begin());
              j != aProcessVector.end(); ++j)
            {
              if(isDependentOn(*j))
                {
                  isExternInterrupted = true;
                }
            }
        }
    }
}

}
