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
		    "]: This process cannot have a HD substrate species: " + 
        getIDString(variableA));
    }
  if(variableB)
    {
      THROW_EXCEPTION(ValueError, String(
        getPropertyInterface().getClassName()) + " [" + getFullID().asString() +
		    "]: This process cannot have a HD substrate species: " + 
        getIDString(variableB));
    }
}

void DiffusionInfluencedReactionProcess::initializeSecond()
{
  ReactionProcess::initializeSecond(); 
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
  setReactMethod();
}

void DiffusionInfluencedReactionProcess::removeMolecule(Species* aSpecies, 
                                                  Voxel* mol,
                                                  const unsigned index) const
{
  if(A != B)
    {
      aSpecies->removeMolecule(index);
    }
  else
    {
      //If A == B, indexB is no longer valid after molA is removed,
      //so need to use the current index to remove molB:
      aSpecies->removeMolecule(mol);
    }
}

void DiffusionInfluencedReactionProcess::removeMolecule(Species* substrate, 
                                                        Voxel* mol,
                                                        unsigned index,
                                                        Species* product) const
{
  if(A == B)
    { 
      index = substrate->getIndex(mol);
    }
  product->addMolecule(mol, substrate->getTag(index));
  substrate->softRemoveMolecule(index);
}

Voxel* DiffusionInfluencedReactionProcess::getPopulatableVoxel(
                                                             Species* aSpecies,
                                                             Voxel* molA,
                                                             Voxel* molB)
{
  Voxel* mol(aSpecies->getRandomAdjoiningVoxel(molA, SearchVacant));
  if(!mol)
    {
      mol = aSpecies->getRandomAdjoiningVoxel(molB, SearchVacant);
    }
  return mol;
}

Voxel* DiffusionInfluencedReactionProcess::getPopulatableVoxel(
                                                             Species* aSpecies,
                                                             Voxel* molA,
                                                             Voxel* molB,
                                                             Voxel* molC)
{
  Voxel* mol(aSpecies->getRandomAdjoiningVoxel(molA, molC, SearchVacant));
  if(!mol)
    {
      mol = aSpecies->getRandomAdjoiningVoxel(molB, molC, SearchVacant);
    }
  return mol;
}

//A + B -> variableC + [D <- molA]
bool DiffusionInfluencedReactionProcess::reactVarC_AtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  variableC->addValue(1);
  D->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  removeMolecule(B, molB, indexB);
  return true;
}

//A + B -> variableC + [D <- molB]
bool DiffusionInfluencedReactionProcess::reactVarC_BtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  variableC->addValue(1);
  D->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  removeMolecule(A, molA, indexA);
  return true;
}

//A + B -> variableC + [D <- molN]
bool DiffusionInfluencedReactionProcess::reactVarC_NtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{ 
  Voxel* mol(getPopulatableVoxel(D, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      variableC->addValue(1);
      //TODO: need to use the correct tag here:
      D->addMolecule(mol, A->getTag(indexA));
      A->removeMolecule(indexA);
      removeMolecule(B, molB, indexB);
      return true;
    }
  return false;
}

//A + B -> variableC + [D == molA]
bool DiffusionInfluencedReactionProcess::reactVarC_AeqD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  variableC->addValue(1);
  B->removeMolecule(indexB);
  return true;
}

//A + B -> variableC + [D == molB]
bool DiffusionInfluencedReactionProcess::reactVarC_BeqD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  variableC->addValue(1);
  A->removeMolecule(indexA);
  return true;
}

//A + B -> variableD + [C <- molA]
bool DiffusionInfluencedReactionProcess::reactVarD_AtoC(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  //cout << getIDString(variableD) << " beg:" << variableD->getValue() <<
  //  std::endl;
  interruptProcessesPre();
  variableD->addValue(1);
  C->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  removeMolecule(B, molB, indexB);
  //cout << getIDString(variableD) << " end:" << variableD->getValue() <<
  //  std::endl;
  return true;
}

//A + B -> variableD + [C <- molB]
bool DiffusionInfluencedReactionProcess::reactVarD_BtoC(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  variableD->addValue(1);
  C->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  removeMolecule(A, molA, indexA);
  return true;
}

//A + B -> variableD + [C <- molN]
bool DiffusionInfluencedReactionProcess::reactVarD_NtoC(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{ 
  Voxel* mol(getPopulatableVoxel(C, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      variableD->addValue(1);
      //TODO: need to use the correct tag here:
      C->addMolecule(mol, A->getTag(indexA));
      A->removeMolecule(indexA);
      removeMolecule(B, molB, indexB);
      return true;
    }
  return false;
}


//A + B -> variableD + [C == A]
bool DiffusionInfluencedReactionProcess::reactVarD_AeqC(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  //cout << "first:" << getIDString() << " A:" << A->getIDString() << " B:" << B->getIDString() << " C:" << C->getIDString() << " vD:" << getIDString(variableD) << std::endl;
  interruptProcessesPre();
  //cout << "dir1" << std::endl;
  variableD->addValue(1);
  //cout << "dir2" << std::endl;
  B->removeMolecule(indexB);
  //cout << "dir3" << std::endl;
  return true;
}

//A + B -> variableD + [C == B]
bool DiffusionInfluencedReactionProcess::reactVarD_BeqC(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  variableD->addValue(1);
  A->removeMolecule(indexA);
  return true;
}

//A + B -> variableC + variableD
bool DiffusionInfluencedReactionProcess::reactVarC_VarD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  variableC->addValue(1);
  variableD->addValue(1);
  A->removeMolecule(indexA);
  removeMolecule(B, molB, indexB);
  return true;
}

//A + B -> variableC
bool DiffusionInfluencedReactionProcess::reactVarC(Voxel* molA,
                                                   Voxel* molB,
                                                   const unsigned indexA,
                                                   const unsigned indexB)
{
  interruptProcessesPre();
  variableC->addValue(1);
  A->removeMolecule(indexA);
  removeMolecule(B, molB, indexB);
  return true;
}

//A + B -> [A == C] + [D <- molB]
bool DiffusionInfluencedReactionProcess::reactAeqC_BtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  D->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  return true;
}

//A + B -> [A == C] + [D <- molN]
bool DiffusionInfluencedReactionProcess::reactAeqC_NtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  Voxel* mol(getPopulatableVoxel(D, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      D->addMolecule(mol, B->getTag(indexB));
      B->removeMolecule(indexB);
      return true;
    }
  return false;
}

//A + B -> [B == C] + [D <- molA]
bool DiffusionInfluencedReactionProcess::reactBeqC_AtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  D->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  return true;
}

//A + B -> [B == C] + [D <- molN]
bool DiffusionInfluencedReactionProcess::reactBeqC_NtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  Voxel* mol(getPopulatableVoxel(D, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      D->addMolecule(mol, A->getTag(indexA));
      A->removeMolecule(indexA);
      return true;
    }
  return false;
}

//A + B -> [C <- molB] + [A == D]
bool DiffusionInfluencedReactionProcess::reactBtoC_AeqD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  C->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  return true;
}

//A + B -> [C <- molN] + [A == D]
bool DiffusionInfluencedReactionProcess::reactNtoC_AeqD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  Voxel* mol(getPopulatableVoxel(C, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      C->addMolecule(mol, B->getTag(indexB));
      B->removeMolecule(indexB);
      return true;
    }
  return false;
}

//A + B -> [C <- molA] + [B == D]
bool DiffusionInfluencedReactionProcess::reactAtoC_BeqD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  C->addMolecule(molA);
  A->softRemoveMolecule(indexA);
  return true;
}

//A + B -> [C <- molA] + [B == D] + [F <- molN]
bool DiffusionInfluencedReactionProcess::reactAtoC_BeqD_NtoF(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  moleculeF = F->getRandomAdjoiningVoxel(molA, SearchVacant);
  if(moleculeF == NULL)
    {
      moleculeF = F->getRandomAdjoiningVoxel(molB, SearchVacant);
      if(moleculeF == NULL)
        {
          return false;
        }
    }
  interruptProcessesPre();
  C->addMolecule(molA);
  A->softRemoveMolecule(indexA);
  F->addMolecule(moleculeF);
  return true;
}


//A + B -> [C <- molA] + [tagC <- tagA] + [B == D]
bool DiffusionInfluencedReactionProcess::reactAtoC_BeqD_tagAtoC(Voxel* molA,
                                                                Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  C->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  return true;
}

//A + B -> [C <- molA] + [E <- compN]
//zero-coefficient E
//we create a molecule E at random location in the compartment to avoid
//rebinding effect, useful to maintain the concentration of a substrate species
//even after the reaction:
bool DiffusionInfluencedReactionProcess::reactAtoC_compNtoE(Voxel* molA,
                                                            Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  Voxel* mol(E->getRandomPopulatableVoxel(1));
  if(mol)
    { 
      interruptProcessesPre();
      E->addMolecule(mol);
      C->addMolecule(molA, A->getTag(indexA));
      A->softRemoveMolecule(indexA);
      removeMolecule(B, molB, indexB);
      return true;
    }
  return false;
}

//A + B -> [C <- molN] + [B == D]
bool DiffusionInfluencedReactionProcess::reactNtoC_BeqD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  Voxel* mol(getPopulatableVoxel(C, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      C->addMolecule(mol, A->getTag(indexA));
      A->removeMolecule(indexA);
      return true;
    }
  return false;
}


//A + B -> [C <- molA] + [D <- molB]
bool DiffusionInfluencedReactionProcess::reactAtoC_BtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  C->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  removeMolecule(B, molB, indexB, D);
  return true;
}

//A + B -> [C <- molR] + [D <- molB]
bool DiffusionInfluencedReactionProcess::reactRtoC_BtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  Species* vac(C->getVacantSpecies());
  unsigned indexV(vac->getRandomValidIndex());
  if(indexV != vac->size())
    {
      A->removeMolecule(indexA);
      C->addMolecule(vac->getMolecule(indexV));
      removeMolecule(B, molB, indexB, D);
      return true;
    }
  return false;
}


//A + B -> [C <- molR] + [D <- molN]
bool DiffusionInfluencedReactionProcess::reactRtoC_NtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  Species* vac(C->getVacantSpecies());
  unsigned indexV(vac->getRandomValidIndex());
  if(indexV != vac->size())
    {
      Voxel* molD(getPopulatableVoxel(D, molA, molB));
      if(molD)
        {
          Voxel* molC(vac->getMolecule(indexV)); 
          if(molC != molD)
            {
              D->addMolecule(molD, B->getTag(indexB));
              A->removeMolecule(indexA);
              C->addMolecule(molC);
              removeMolecule(B, molB, indexB);
              return true;
            }
        }
    }
  return false;
}


//A + B => B + A
//A + B -> [A swap B]
bool DiffusionInfluencedReactionProcess::swapAB(Voxel* molA,
                                                Voxel* molB,
                                                const unsigned indexA,
                                                const unsigned indexB)
{
  interruptProcessesPre();
  Tag tagB(A->getTag(indexA));
  Tag tagA(B->getTag(indexB));
  unsigned oriA(tagA.boundCnt);
  unsigned oriB(tagB.boundCnt);
  unsigned boundCnt(tagA.boundCnt);
  tagA.boundCnt = tagB.boundCnt;
  tagB.boundCnt = boundCnt;
  A->softReplaceMolecule(indexA, molB, tagA, B);
  B->softReplaceMolecule(indexB, molA, tagB, A);
  return true;
}


//A + B => C + A
//A + B -> [C <- molA] + [move A to D]
bool DiffusionInfluencedReactionProcess::moveAtoD_reactBtoC(Voxel* molA,
                                                          Voxel* molB,
                                                          const unsigned indexA,
                                                          const unsigned indexB)
{
  interruptProcessesPre();
  Tag tagD(A->getTag(indexA));
  tagD.boundCnt = B->getTag(indexB).boundCnt;
  Tag tagC(A->getTag(indexA));
  D->softReplaceMolecule(indexA, molB, tagD, B);
  C->addMolecule(molA, tagC);
  B->softRemoveMolecule(indexB);
  return true;
}

//A + B => B + D
//A + B -> [move B to C] + [D <- molB]
bool DiffusionInfluencedReactionProcess::reactAtoD_moveBtoC(Voxel* molA,
                                                          Voxel* molB,
                                                          const unsigned indexA,
                                                          const unsigned indexB)
{
  interruptProcessesPre();
  Tag tagC(B->getTag(indexB));
  tagC.boundCnt = A->getTag(indexA).boundCnt;
  Tag tagD(B->getTag(indexB));
  C->softReplaceMolecule(indexB, molA, tagC, A);
  D->addMolecule(molB, tagD);
  A->softRemoveMolecule(indexA);
  return true;
}


//A + B -> [C <- molA] + [D <- molN]
bool DiffusionInfluencedReactionProcess::reactAtoC_NtoD(
                                                  Voxel* molA, Voxel* molB,
                                                  const unsigned indexA,
                                                  const unsigned indexB)
{
  Voxel* mol(getPopulatableVoxel(D, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      D->addMolecule(mol, B->getTag(indexB));
      C->addMolecule(molA, A->getTag(indexA));
      A->softRemoveMolecule(indexA);
      removeMolecule(B, molB, indexB);
      return true;
    }
  return false;
}

//A + B -> [C <- molB] + [D <- molA]
bool DiffusionInfluencedReactionProcess::reactBtoC_AtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  interruptProcessesPre();
  C->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  removeMolecule(A, molA, indexA, D);
  return true;
}

//A + B -> [C <- molB] + [D <- molN]
bool DiffusionInfluencedReactionProcess::reactBtoC_NtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  Voxel* mol(getPopulatableVoxel(D, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      D->addMolecule(mol, A->getTag(indexA));
      C->addMolecule(molB, B->getTag(indexB));
      B->softRemoveMolecule(indexB);
      removeMolecule(A, molA, indexA);
      return true;
    }
  return false;
}

//A + B -> [C <- molN] + [D <- molN]
bool DiffusionInfluencedReactionProcess::reactNtoC_NtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  Voxel* molC(getPopulatableVoxel(C, molA, molB));
  if(molC)
    {
      Voxel* molD(getPopulatableVoxel(C, molA, molB, molC));
      if(molD)
        {
          interruptProcessesPre();
          C->addMolecule(molC, A->getTag(indexA));
          D->addMolecule(molD, B->getTag(indexB));
          A->removeMolecule(indexA);
          removeMolecule(B, molB, indexB);
          return true;
        }
    }
  return false;
}

//A + B -> [C <- molN] + [D <- molB]
bool DiffusionInfluencedReactionProcess::reactNtoC_BtoD(Voxel* molA,
                                                        Voxel* molB,
                                                        const unsigned indexA,
                                                        const unsigned indexB)
{
  Voxel* mol(getPopulatableVoxel(C, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      C->addMolecule(mol, A->getTag(indexA));
      D->addMolecule(molB, B->getTag(indexB));
      A->removeMolecule(indexA);
      removeMolecule(B, molB, indexB);
      return true;
    }
  return false;
}

//A + B -> [A == C]
bool DiffusionInfluencedReactionProcess::reactAeqC(Voxel* molA, Voxel* molB,
                                                   const unsigned indexA,
                                                   const unsigned indexB)
{
  interruptProcessesPre();
  B->removeMolecule(indexB);
  return true;
}

//A + B -> [B == C]
bool DiffusionInfluencedReactionProcess::reactBeqC(Voxel* molA, Voxel* molB,
                                                   const unsigned indexA,
                                                   const unsigned indexB)
{
  interruptProcessesPre();
  A->removeMolecule(indexA);
  return true;
}

//A + B -> [C <- molA]
bool DiffusionInfluencedReactionProcess::reactAtoC(Voxel* molA, Voxel* molB,
                                                   const unsigned indexA,
                                                   const unsigned indexB)
{
  interruptProcessesPre();
  C->addMolecule(molA, A->getTag(indexA));
  A->softRemoveMolecule(indexA);
  removeMolecule(B, molB, indexB);
  return true;
}

//A + B -> [C <- molB]
bool DiffusionInfluencedReactionProcess::reactBtoC(Voxel* molA, Voxel* molB,
                                                   const unsigned indexA,
                                                   const unsigned indexB)
{
  interruptProcessesPre();
  C->addMolecule(molB);
  B->softRemoveMolecule(indexB);
  removeMolecule(A, molA, indexA);
  return true;
}

//A + B -> [C <- molB] + [tagC <- tagA]
bool DiffusionInfluencedReactionProcess::reactBtoC_tagAtoC(Voxel* molA,
                                                           Voxel* molB,
                                                   const unsigned indexA,
                                                   const unsigned indexB)
{
  interruptProcessesPre();
  C->addMolecule(molB, A->getTag(indexA));
  B->softRemoveMolecule(indexB);
  removeMolecule(A, molA, indexA);
  return true;
}

//A + B -> [C <- molB] + [tagC <- tagB]
bool DiffusionInfluencedReactionProcess::reactBtoC_tagBtoC(Voxel* molA,
                                                           Voxel* molB,
                                                   const unsigned indexA,
                                                   const unsigned indexB)
{
  interruptProcessesPre();
  C->addMolecule(molB, B->getTag(indexB));
  B->softRemoveMolecule(indexB);
  removeMolecule(A, molA, indexA);
  return true;
}

//A + B -> [C <- molN]
bool DiffusionInfluencedReactionProcess::reactNtoC(Voxel* molA, Voxel* molB,
                                                   const unsigned indexA,
                                                   const unsigned indexB)
{
  Voxel* mol(getPopulatableVoxel(C, molA, molB));
  if(mol)
    {
      interruptProcessesPre();
      C->addMolecule(mol, A->getTag(indexA));
      A->removeMolecule(indexA);
      removeMolecule(B, molB, indexB);
      return true;
    }
  return false;
}

void DiffusionInfluencedReactionProcess::setReactMethod()
{
  if(RandomC)
    {
      ForcedSequence = true;
    }
  if(ForcedSequence)
    {
      setForcedSequenceReactMethod();
    }
  else
    {
      setFreeSequenceReactMethod();
    }
}

void DiffusionInfluencedReactionProcess::setForcedSequenceReactMethod()
{
  if(!A || !B || !C || !D)
    {
      THROW_EXCEPTION(ValueError,
                   String(getPropertyInterface().getClassName()) +
                  "[" + getFullID().asString() + "]: For DIRP ForcedSequence" +
                  "or RandomC, A, B, C and D must all be nonHD species.");
    }
  setGeneralForcedSequenceReactMethod();
  /*
  if(C && D)
    {
      if(A == D)
        {
          //A + B => B + A
          if(B == C)
            {
              if(A->getIsDiffusing() || B->getIsDiffusing())
                {
                  //A + B -> [A swap B]
                  setGeneralForcedSequenceReactMethod();
                  //reactM = &DiffusionInfluencedReactionProcess::swapAB;
                }
              else
                {
                  setGeneralForcedSequenceReactMethod();
                }
            }
          //A + B => C + A
          else if(A->getIsDiffusing()) 
            {
              //A + B -> [C <- molA] + [move A to D]
              setGeneralForcedSequenceReactMethod();
              //reactM = &DiffusionInfluencedReactionProcess::moveAtoD_reactBtoC;
            }
          else
            {
              setGeneralForcedSequenceReactMethod();
            }/mo
        }
      //A + B => B + D
      else if(B == C)
        {
          if(B->getIsDiffusing())
            {
               //A + B -> [move B to C] + [D <- molB]
               setGeneralForcedSequenceReactMethod();
               //reactM = &DiffusionInfluencedReactionProcess::reactAtoD_moveBtoC;
            }
          else
            {
              setGeneralForcedSequenceReactMethod();
            }
        }
      else
        {
          setGeneralForcedSequenceReactMethod();
        }
    }
    */
}

void DiffusionInfluencedReactionProcess::setGeneralForcedSequenceReactMethod()
{
  if(C && D)
    {
      if(RandomC)
        {
          if(B->isReplaceable(D))
            { 
              //A + B -> [C <- molR] + [D <- molB]
              reactM = &DiffusionInfluencedReactionProcess::reactRtoC_BtoD;
            }
          else
            {
              //A + B -> [C <- molR] + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactRtoC_NtoD;
            }

        }
      else if(A->isReplaceable(C))
        {
          if(B->isReplaceable(D))
            {
              if(B == D && F)
                {
                  //A + B -> [C <- molA] + [B == D] + [D <- molN]
                  reactM = 
                    &DiffusionInfluencedReactionProcess::reactAtoC_BeqD_NtoF;
                }
              else
                {
                  //A + B -> [C <- molA] + [D <- molB]
                  reactM = &DiffusionInfluencedReactionProcess::reactAtoC_BtoD;
                }
            }
          else
            {
              //A + B -> [C <- molA] + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactAtoC_NtoD;
            }
        }
      else
        {
          if(B->isReplaceable(D))
            {
              //A + B -> [C <- molN] + [D <- molB]
              reactM = &DiffusionInfluencedReactionProcess::reactNtoC_BtoD;
            }
          else
            {
              //A + B -> [C <- molN] + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactNtoC_NtoD;
            }
        }
    }
}

void DiffusionInfluencedReactionProcess::setFreeSequenceReactMethod()
{
  if(variableC && D)
    {
      if(A == D)
        {
          //A + B -> variableC + [A == D]
          reactM = &DiffusionInfluencedReactionProcess::reactVarC_AeqD;
        }
      else if(B == D)
        {
          //A + B -> variableC + [B == D]
          reactM = &DiffusionInfluencedReactionProcess::reactVarC_BeqD;
        }
      else
        { 
          if(A->isReplaceable(D))
            {
              //A + B -> variableC + [D <- molA]
              reactM = &DiffusionInfluencedReactionProcess::reactVarC_AtoD;
            }
          else if(B->isReplaceable(D))
            {
              //A + B -> variableC + [D <- molB]
              reactM = &DiffusionInfluencedReactionProcess::reactVarC_BtoD;
            }
          else
            {
              //A + B -> variableC + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactVarC_NtoD;
            }
        }
    }
  else if(variableD && C)
    {
      if(A == C)
        {
          //cout << "reactVarD_AeqC" << std::endl;
          //A + B -> variableD + [A == C]
          reactM = &DiffusionInfluencedReactionProcess::reactVarD_AeqC;
        }
      else if(B == C)
        {
          //A + B -> variableD + [B == C]
          reactM = &DiffusionInfluencedReactionProcess::reactVarD_BeqC;
        }
      else
        { 
          if(A->isReplaceable(C))
            {
              //A + B -> variableD + [C <- molA]
              reactM = &DiffusionInfluencedReactionProcess::reactVarD_AtoC;
            }
          else if(B->isReplaceable(C))
            {
              //A + B -> variableD + [C <- molB]
              reactM = &DiffusionInfluencedReactionProcess::reactVarD_BtoC;
            }
          else
            {
              //A + B -> variableD + [C <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactVarD_NtoC;
            }
        }
    }
  else if(variableC)
    {
      if(variableD)
        {
          //A + B -> variableC + variableD
          reactM = &DiffusionInfluencedReactionProcess::reactVarC_VarD;
        }
      else
        {
          //A + B -> variableC
          reactM = &DiffusionInfluencedReactionProcess::reactVarC;
        }
    }
  //A + B -> C + D
  else if(D)
    {
      if(A == C && B == D)
        {
          //A + B -> [A == C] + [B == D]
          reactM = &DiffusionInfluencedReactionProcess::reactNone;
        }
      else if(B == C && A == D)
        {
          //A + B -> [B == C] + [A == D]
          reactM = &DiffusionInfluencedReactionProcess::reactNone;
        }
      else if(A == C)
        {
          if(B->isReplaceable(D))
            {
              //cout << "reactAeqC_BtoD:" << getIDString() << std::endl;
              //A + B -> [A == C] + [D <- molB]
              reactM = &DiffusionInfluencedReactionProcess::reactAeqC_BtoD;
            }
          else
            {
              //A + B -> [A == C] + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactAeqC_NtoD;
            }
        }
      else if(B == C)
        {
          if(A->isReplaceable(D))
            {
              //cout << "reactBeqC_AtoD:" << getIDString() << std::endl;
              //A + B -> [B == C] + [D <- molA]
              reactM = &DiffusionInfluencedReactionProcess::reactBeqC_AtoD;
            }
          else
            {
              //A + B -> [B == C] + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactBeqC_NtoD;
            }
        }
      else if(A == D)
        {
          if(B->isReplaceable(C))
            {
              //cout << "reactBtoC_AeqD:" << getIDString() << std::endl;
              //A + B -> [C <- molB] + [A == D]
              reactM = &DiffusionInfluencedReactionProcess::reactBtoC_AeqD;
            }
          else
            {
              //cout << "reactNtoC_AeqD:" << getIDString() << std::endl;
              //A + B -> [C <- molN] + [A == D]
              reactM = &DiffusionInfluencedReactionProcess::reactNtoC_AeqD;
            }
        }
      else if(B == D)
        {
          if(A->isReplaceable(C))
            {
              if(C->getIsTagged() && A->getIsTagged())
                {
              //cout << "reactAtoC_BeqD_tagAtoC:" << getIDString() << std::endl;
                  //A + B -> [C <- molA] + [tagC <- tagA] + [B == D]
                  reactM = 
                    &DiffusionInfluencedReactionProcess::reactAtoC_BeqD_tagAtoC;
                }
              else
                {
                  if(F)
                    {
                      //cout << "reactAtoC_BeqD_NtoF:" << getIDString() <<
                      //  std::endl;
                      //A + B -> [C <- molA] + [B == D] + [D <- molN]
                      reactM = 
                        &DiffusionInfluencedReactionProcess::reactAtoC_BeqD_NtoF;
                    }
                  else
                    {
                      //cout << "reactAtoC_BeqD:" << getIDString() <<
                      //std::endl;
                      //A + B -> [C <- molA] + [B == D]
                      reactM = 
                        &DiffusionInfluencedReactionProcess::reactAtoC_BeqD;
                    }
                }
            }
          else
            {
              //A + B -> [C <- molN] + [B == D]
              reactM = &DiffusionInfluencedReactionProcess::reactNtoC_BeqD;
            }
        }
      else
        {
          if(A->isReplaceable(C))
            {
              if(B->isReplaceable(D))
                {
                  //A + B -> [C <- molA] + [D <- molB]
                  reactM = &DiffusionInfluencedReactionProcess::reactAtoC_BtoD;
                }
              else
                {
                  //A + B -> [C <- molA] + [D <- molN]
                  reactM = &DiffusionInfluencedReactionProcess::reactAtoC_NtoD;
                }
            }
          else if(B->isReplaceable(C))
            {
              if(A->isReplaceable(D))
                {
                  //A + B -> [C <- molB] + [D <- molA]
                  reactM = &DiffusionInfluencedReactionProcess::reactBtoC_AtoD;
                }
              else
                {
                  //A + B -> [C <- molB] + [D <- molN]
                  reactM = &DiffusionInfluencedReactionProcess::reactBtoC_NtoD;
                }
            }
          else
            {
              //A + B -> [C <- molN] + [D <- molN]
              reactM = &DiffusionInfluencedReactionProcess::reactNtoC_NtoD;
            }
        }
    }
  //A + B -> C + E(0 coefficient, random comp voxel)
  else if(E)
    {
      if(A == C)
        {
          throwException("reactAeqC_E");
        }
      else if(B == C)
        {
          throwException("reactBeqC_E");
        }
      else if(A->isReplaceable(C))
        {
          reactM = &DiffusionInfluencedReactionProcess::reactAtoC_compNtoE;
        }
      else
        {
          throwException("reactBtoC_E");
        }
    }
  else
    {
      if(A == C)
        {
          //A + B -> [A == C]
          reactM = &DiffusionInfluencedReactionProcess::reactAeqC;
        }
      else if(B == C)
        {
          //A + B -> [B == C]
          reactM = &DiffusionInfluencedReactionProcess::reactBeqC;
        }
      else if(A->isReplaceable(C))
        {
          //A + B -> [C <- molA]
          reactM = &DiffusionInfluencedReactionProcess::reactAtoC;
        }
      else if(B->isReplaceable(C))
        {
          if(C->getIsTagged())
            {
              if(B->getIsTagged())
                {
                  //A + B -> [C <- molB] + [tagC <- tagB]
                  reactM = 
                    &DiffusionInfluencedReactionProcess::reactBtoC_tagBtoC;
                }
              else
                {
                  //A + B -> [C <- molB] + [tagC <- tagA]
                  reactM =
                    &DiffusionInfluencedReactionProcess::reactBtoC_tagAtoC;
                }
            }
          else
            {
              //A + B -> [C <- molB]
              reactM = &DiffusionInfluencedReactionProcess::reactBtoC;
            }
        }
      else
        {
          //A + B -> [C <- molN]
          reactM = &DiffusionInfluencedReactionProcess::reactNtoC;
        }
    }
}

void DiffusionInfluencedReactionProcess::throwException(String aString)
{
  THROW_EXCEPTION(ValueError, String(getPropertyInterface().getClassName()) +
                  "[" + getFullID().asString() + "]: " + aString + " is not " +
                  "yet implemented.");
}

//positive-coefficient F
void DiffusionInfluencedReactionProcess::addMoleculeF()
{
  if(!F)
    {
      return;
    }
  moleculeF = F->getRandomAdjoiningVoxel(moleculeC, SearchVacant);
  if(moleculeF == NULL)
    {
      moleculeF = F->getRandomAdjoiningVoxel(moleculeD, SearchVacant);
      if(moleculeF == NULL)
        {
          return;
        }
    }
  F->addMolecule(moleculeF);
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
  else if(A->getDimension() == 3 && B->getIsCompVacant() && 
          B->getDimension() == 2)
    {
      double nv(A->getComp()->vacantSpecies->size());
      double ns(B->size());
      cout << "1st ns:" << ns << " " << getIDString() << std::endl;
      /*
      if(B->getComp()->interfaceID != theSpecies.size())
        {
          ns = theSpecies[B->getComp()->interfaceID]->size();
        }
        */
      double S(B->getComp()->specArea);
      double V(A->getComp()->actualVolume);
      cout << "1: nv:" << nv << " ns:" << ns << " S:" << S << " V:" << V << std::endl;
      if(p == -1)
        {
          //p = 24*k*r_v/((6+3*sqrt(3)+2*sqrt(6))*D_A); //averaged hcp surface
          //p = 4*k*sqrt(3)*r_v/(3*sqrt(2)*D_A); //perfect hcp surface
          p = 8*k*nv*r_v*r_v*S/(3*ns*V*D_A); //more accurate for full lattice
        }
      else
        {
          //k = p*((6+3*sqrt(3)+2*sqrt(6))*D_A)/(24*r_v); //averaged hcp surf
          k = p*3*ns*V*D_A/(8*nv*r_v*r_v*S); //more accurate for full lattice
        }
    }
  else if(A->getIsCompVacant() && A->getDimension() == 2 &&
          B->getDimension() == 3)
    {
      double nv(B->getComp()->vacantSpecies->size());
      double ns(A->size());
      cout << "2nd ns:" << ns << " " << getIDString() << std::endl;
      /*
      if(A->getComp()->interfaceID != theSpecies.size())
        {
          ns = theSpecies[A->getComp()->interfaceID]->size();
        }
        */
      double S(A->getComp()->specArea);
      double V(B->getComp()->actualVolume);
      cout << "1: nv:" << nv << " ns:" << ns << " S:" << S << " V:" << V << std::endl;
      if(p == -1)
        {
          //p = 24*k*r_v/((6+3*sqrt(3)+2*sqrt(6))*D_B); //averaged hcp surface
          //p = 4*k*sqrt(3)*r_v/(3*sqrt(2)*D_B); //perfect hcp surface
          p = 8*k*nv*r_v*r_v*S/(3*ns*V*D_B); //more accurate for full lattice
        }
      else
        {
          //k = p*((6+3*sqrt(3)+2*sqrt(6))*D_B)/(24*r_v); //averaged hcp surf
          k = p*3*ns*V*D_B/(8*nv*r_v*r_v*S); //more accurate for full lattice
        }
    }
  else if(A->getDimension() == 3 && B->getIsCompVacant() && 
          B->getDimension() == 1)
    {
      double nv(A->getComp()->vacantSpecies->compVoxelSize());
      double nl(B->compVoxelSize());
      cout << "31 1st nl:" << nl << " " << getIDString() << std::endl;
      /*
      if(B->getComp()->interfaceID != theSpecies.size())
        {
          nl = theSpecies[B->getComp()->interfaceID]->size();
        }
        */
      double L(B->getComp()->specLength);
      double V(A->getComp()->actualVolume);
      cout << "1: nv:" << nv << " nl:" << nl << " L:" << L << " V:" <<
        V << std::endl;
      if(p == -1)
        {
          if(B->getComp()->interfaceID != theSpecies.size())
            {
              p = k/D_A;
              //p = k*L/(D_A*nl);
              Species* interface(theSpecies[B->getComp()->interfaceID]);
              cout << "process:" << getIDString() << " " << getIDString(interface) << " interface size:" << interface->size() << std::endl;
              for(unsigned i(0); i != interface->size(); ++i)
                {
                  Voxel* mol(interface->getMolecule(i));
                  for(unsigned j(0); j != mol->adjoiningSize-mol->diffuseSize;
                      ++j)
                    {
                      double val(interface->getInterfaceConst(mol, j));
                      if(val*p > 1)
                        {
                          cout << i << " " << j << " " << val*p << std::endl;
                        }
                    }
                }
            }
          else
            { 
              p = 4*k*nv*r_v*r_v*L/(5*nl*V*D_A);
            }
        }
      else
        {
          // [m^s/s] = [m^3][m^2/s]/([m][m][m]) = [m^2/s]
          k = p*5*nl*V*D_A/(4*nv*r_v*r_v*L);
          if(B->getComp()->interfaceID != theSpecies.size())
            {
              Species* interface(theSpecies[B->getComp()->interfaceID]);
              interface->setUnityInterfaceConsts();
            }
        }
    }
  else if(A->getIsCompVacant() && A->getDimension() == 1 &&
          B->getDimension() == 3)
    {
      double nv(B->getComp()->vacantSpecies->compVoxelSize());
      double nl(A->compVoxelSize());
      cout << "13 2nd nl:" << nl << " " << getIDString() << std::endl;
      /*
      if(A->getComp()->interfaceID != theSpecies.size())
        {
          nl = theSpecies[A->getComp()->interfaceID]->size();
        }
        */
      double L(A->getComp()->specLength);
      double V(B->getComp()->actualVolume);
      cout << "1: nv:" << nv << " nl:" << nl << " L:" << L << " V:" << 
        V << std::endl;
      if(p == -1)
        {
          if(A->getComp()->interfaceID != theSpecies.size())
            {
              p = k/D_B;
              //p = k*L/(D_B*nl);
              cout << "process:" << getIDString() << std::endl;
              Species* interface(theSpecies[A->getComp()->interfaceID]);
              for(unsigned i(0); i != interface->size(); ++i)
                {
                  Voxel* mol(interface->getMolecule(i));
                  for(unsigned j(0); j != mol->adjoiningSize-mol->diffuseSize;
                      ++j)
                    {
                      double val(interface->getInterfaceConst(mol, j));
                      if(val*p > 1)
                        {
                          cout << i << " " << j << " " << val*p << std::endl;
                        }
                    }
                }
            }
          else
            {
              p = 4*k*nv*r_v*r_v*L/(5*nl*V*D_B);
            }
        }
      else
        {
          // [m^s/s] = [m^3][m^2/s]/([m][m][m]) = [m^2/s]
          k = p*5*nl*V*D_B/(4*nv*r_v*r_v*L);
          if(A->getComp()->interfaceID != theSpecies.size())
            {
              Species* interface(theSpecies[A->getComp()->interfaceID]);
              interface->setUnityInterfaceConsts();
            }
        }
    }
  else if(A->getDimension() == 3 && B->getDimension() == 2)
    {
      /*
      //Need to use nInterface/nVacant ratio to improve accuracy when a 
      //molecule in an off-lattice compartment is connected by multiple
      //interface voxels.
      double i_v(1);
      if(B->getComp()->interfaceID != theSpecies.size())
        {
          double nInterface(theSpecies[B->getComp()->interfaceID]->size());
          double nVacant(B->getVacantSpecies()->compVoxelSize());
          i_v = nInterface/nVacant;
        }
        */
      if(p == -1)
        {
          p = sqrt(2)*k/(3*D_A*r_v);
        }
      else
        {
          k = p*(3*D_A*r_v)/sqrt(2); //[m^3/s]
        }
    }
  else if(A->getDimension() == 2 && B->getDimension() == 3)
    {
      /*
      //Need to use nInterface/nVacant ratio to improve accuracy when a 
      //molecule in an off-lattice compartment is connected by multiple
      //interface voxels.
      double i_v(1);
      if(A->getComp()->interfaceID != theSpecies.size())
        {
          double nInterface(theSpecies[A->getComp()->interfaceID]->size());
          double nVacant(A->getVacantSpecies()->compVoxelSize());
          i_v = nInterface/nVacant;
        }
        */
      if(p == -1)
        {
          p = sqrt(2)*k/(3*D_B*r_v);
        }
      else
        {
          k = p*(3*D_B*r_v)/sqrt(2); //[m^3/s]
        }
    }
  else if(A->getDimension() == 3 && B->getDimension() == 1)
    {
      double nv(A->getComp()->vacantSpecies->compVoxelSize());
      double V(A->getComp()->actualVolume);
      /*
      //Need to use nInterface/nVacant ratio to improve accuracy when a 
      //molecule in an off-lattice compartment is connected by multiple
      //interface voxels.
      double i_v(1);
      if(B->getComp()->interfaceID != theSpecies.size())
        {
          double nInterface(theSpecies[B->getComp()->interfaceID]->size());
          double nVacant(B->getVacantSpecies()->compVoxelSize());
          i_v = nInterface/nVacant;
        }
        */
      if(p == -1)
        {
          if(B->getComp()->interfaceID != theSpecies.size())
            {
              p = k/D_A;
              cout << "process:" << getIDString() << std::endl;
              Species* interface(theSpecies[B->getComp()->interfaceID]);
              for(unsigned i(0); i != interface->size(); ++i)
                {
                  Voxel* mol(interface->getMolecule(i));
                  for(unsigned j(0); j != mol->adjoiningSize-mol->diffuseSize;
                      ++j)
                    {
                      double val(interface->getInterfaceConst(mol, j));
                      if(val*p > 1)
                        {
                          cout << i << " " << j << " " << val*p << std::endl;
                        }
                    }
                }
            }
          else
            {
              p = 4*k*nv*r_v*r_v/(5*V*D_A);
            }
        }
      else
        {
          k = p*5*V*D_A/(4*nv*r_v*r_v); //[m^3/s]
          if(B->getComp()->interfaceID != theSpecies.size())
            {
              Species* interface(theSpecies[B->getComp()->interfaceID]);
              interface->setUnityInterfaceConsts();
            }
        }
    }
  else if(A->getDimension() == 1 && B->getDimension() == 3)
    {
      double nv(B->getComp()->vacantSpecies->compVoxelSize());
      double V(B->getComp()->actualVolume);
      /*
      //Need to use nInterface/nVacant ratio to improve accuracy when a 
      //molecule in an off-lattice compartment is connected by multiple
      //interface voxels.
      double i_v(1);
      if(A->getComp()->interfaceID != theSpecies.size())
        {
          double nInterface(theSpecies[A->getComp()->interfaceID]->size());
          double nVacant(A->getVacantSpecies()->compVoxelSize());
          i_v = nInterface/nVacant;
        }
        */
      if(p == -1)
        {
          if(A->getComp()->interfaceID != theSpecies.size())
            {
              p = k/D_B;
              cout << "process:" << getIDString() << std::endl;
              Species* interface(theSpecies[A->getComp()->interfaceID]);
              for(unsigned i(0); i != interface->size(); ++i)
                {
                  Voxel* mol(interface->getMolecule(i));
                  for(unsigned j(0); j != mol->adjoiningSize-mol->diffuseSize;
                      ++j)
                    {
                      double val(interface->getInterfaceConst(mol, j));
                      if(val*p > 1)
                        {
                          cout << i << " " << j << " " << val*p << std::endl;
                        }
                    }
                }
            }
          else
            {
              p = 4*k*nv*r_v*r_v/(5*V*D_B);
            }
        }
      else
        {
          k = p*5*V*D_B/(4*nv*r_v*r_v); //[m^3/s]
          if(A->getComp()->interfaceID != theSpecies.size())
            {
              Species* interface(theSpecies[A->getComp()->interfaceID]);
              interface->setUnityInterfaceConsts();
            }
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
  cout << aProcess << std::endl;
  cout << "  " << getIDString(A) << "[dim:" << A->getDimension() << "]"  <<
    " + " <<  getIDString(B) << "[dim:" << B->getDimension() << "]" << " -> ";
  if(C)
    {
      cout << getIDString(C) << "[dim:" << C->getDimension() << "]";
    }
  else
    {
      cout << getIDString(variableC);
    }
  if(D)
    {
      cout << " + " << getIDString(D) << "[dim:" << D->getDimension() << "]";
    }
  else if(variableD)
    {
      cout << " + " << getIDString(variableD);
    }
  cout << ": k=" << k << ", p=" << p << 
    ", p_A=" << A->getReactionProbability(B->getID()) <<
    ", p_B=" << B->getReactionProbability(A->getID()) << std::endl; 
}

}
