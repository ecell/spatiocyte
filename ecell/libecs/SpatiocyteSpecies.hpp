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


#ifndef __SpatiocyteSpecies_hpp
#define __SpatiocyteSpecies_hpp

#include <sstream>
#include <libecs/Variable.hpp>
#include <libecs/SpatiocyteCommon.hpp>
#include <libecs/SpatiocyteStepper.hpp>
#include <libecs/SpatiocyteNextReactionProcess.hpp>
#include <libecs/DiffusionInfluencedReactionProcess.hpp>
#include <libecs/MoleculePopulateProcess.hpp>
#include <libecs/SpatiocyteVector.hpp>

namespace libecs
{

using namespace libecs;

// The size of Coord must be 128 bytes to avoid cacheline splits
// The Core 2 has 64-byte cacheline


static String int2str(int anInt)
{
  std::stringstream aStream;
  aStream << anInt;
  return aStream.str();
}

static String double2str(double aVal)
{
  std::stringstream aStream;
  aStream << aVal;
  return aStream.str();
}

/*
 * isVacant: vacant species definition: a species on which other species can
 * diffuse on or occupy. One major optimization in speed for vacant species is
 * that theMolecules list is not updated when other species diffuse on it. 
 * There are four possible types of vacant species:
 * 1. !isDiffusiveVacant && !isReactiveVacant: a vacant species
 *    whose theMolecules list and theMoleculeSize are never updated. This is
 *    the most basic version. This vacant species never diffuses, or reacts
 *    using SNRP. 
 * 2. !isDiffusiveVacant && isReactiveVacant: a vacant species which is a
 *    substrate of a SNRP reaction. In this case theMoleculeSize is updated
 *    before the step interval of SNRP is calculated, and theMolecules list is
 *    updated before the SNRP reaction (fire) is executed to get a valid
 *    molecule list. Set by SNRP.
 * 3. isDiffusiveVacant && !isReactiveVacant: a vacant species which also
 *    diffuses. In this case, theMolecules list and theMoleculeSize are
 *    updated just before it is diffused. Set by DiffusionProcess.   
 * 4. isDiffusiveVacant && isReactiveVacant: a vacant species which reacts
 *    using SNRP and also diffuses.
 * 5. isCompVacant: the VACANT species declared for each compartment. It can
 *    also be isDiffusiveVacant and isReactiveVacant. Set during compartment
 *    registration. It also persistently stores all the compartment voxels.
 *    Referred to as theVacantSpecies. For isCompVacant, initially there
 *    are no molecules in its list. All voxels are stored in theCompVoxels. Only
 *    if it is updated before being called by VisualizationLogProcess, SNRP or
 *    DiffusionProcess, it will be theMolecules will be populated with the
 *    CompCoords.
 * 6. isVacant {isCompVacant; isDiffusiveVacant; isReactiveVacant): the
 *    general name used to identify either isCompVacant, isDiffusiveVacant or
 *    isReactiveVacant to reduce comparison operations.
 */

class Species
{
public:
  Species(SpatiocyteStepper* aStepper, Variable* aVariable, int anID, 
          int anInitCoordSize, RandomLib::Random& aRng,
          double voxelRadius, std::vector<Voxel>& aLattice,
          std::vector<Species*>& aSpeciesList):
    isCentered(false),
    isCompVacant(false),
    isDeoligomerize(false),
    isDiffusing(false),
    isDiffusiveVacant(false),
    isFixedAdjoins(false),
    isGaussianPopulation(false),
    isInContact(false),
    isInterface(false),
    isMultiMultiReactive(false),
    isMultiscale(false),
    isMultiscaleComp(false),
    isOffLattice(false),
    isOnMultiscale(false),
    isOrigins(false),
    isPeriodic(false),
    isPolymer(false),
    isReactiveVacant(false),
    isRegularLattice(false),
    isSubunitInitialized(false),
    isTag(false),
    isTagged(false),
    isVacant(false),
    theID(anID),
    theInitCoordSize(anInitCoordSize),
    theMoleculeSize(0),
    theRotateSize(1),
    theSpeciesCollisionCnt(0),
    D(0),
    nDiffuseRadius(0.5),
    theDiffuseRadius(voxelRadius),
    theDiffusionInterval(libecs::INF),
    theMoleculeRadius(voxelRadius),
    theVoxelRadius(voxelRadius),
    theWalkProbability(1),
    theRng(aRng),
    theTrailSpecies(NULL),
    theDeoligomerizedProduct(NULL),
    thePopulateProcess(NULL),
    theStepper(aStepper),
    theVariable(aVariable),
    cout(std::cout),
    theCompVoxels(&theMolecules),
    theLattice(aLattice),
    theSpecies(aSpeciesList) {}
  ~Species() {}
  void initialize(const int anAdjoiningCoordSize, const unsigned aNullCoord,
                  const unsigned aNullID)
    {
      const unsigned speciesSize(theSpecies.size());
      theAdjoiningCoordSize = anAdjoiningCoordSize;
      theNullCoord = aNullCoord;
      theNullID = aNullID;
      theStride = UINT_MAX/speciesSize;
      theReactionProbabilities.resize(speciesSize, 0);
      theDiffusionInfluencedReactions.resize(speciesSize, NULL);
      isFinalizeReactions.resize(speciesSize, false);
      isMultiscaleBinderID.resize(speciesSize, false);
      isMultiMultiReactantID.resize(speciesSize, false);
      isMultiscaleBoundID.resize(speciesSize, false);
      isProductPair.resize(speciesSize, false);
      theMultiscaleUnbindIDs.resize(speciesSize);
      if(theComp)
        {
          setVacantSpecies(theComp->vacantSpecies);
        }
      theVacantIdx = theStride*theVacantID;
      theNullTag.speciesID = theNullID;
      theNullTag.rotIndex = 0;
      theNullTag.multiIdx = 0;
      theNullTag.boundCnt = 0;
      cout.setLevel(theStepper->getDebugLevel());
    }
  void setDiffusionInfluencedReaction(DiffusionInfluencedReactionProcess*
                                      aReaction, int anID, double aProbability)
    {
      if(theDiffusionInfluencedReactions[anID] != NULL &&
         aReaction->getIDString() != 
         theDiffusionInfluencedReactions[anID]->getIDString())
        {
          THROW_EXCEPTION(ValueError, getIDString() +
             ": has duplicate DiffusionInfluencedReactionProcesses (1):" +
             aReaction->getIDString() + " (2):" +
             theDiffusionInfluencedReactions[anID]->getIDString() +
             ". Remove one of them before proceeding.");
        }
      theDiffusionInfluencedReactions[anID] = aReaction;
      theReactionProbabilities[anID] = aProbability;
    }
  void setDiffusionInfluencedReactantPair(Species* aSpecies)
    {
      theDiffusionInfluencedReactantPairs.push_back(aSpecies);
    }
  void setPopulateProcess(MoleculePopulateProcess* aProcess, double aDist)
    {
      if(aDist)
        {
          isGaussianPopulation = true;
        }
      thePopulateProcess = aProcess;
    }
  void setProductPair(Species* aSpecies)
    {
      isProductPair[aSpecies->getID()] = true;
    }
  bool isTrailSpecies(Species* aSpecies)
    {
      return (theTrailSpecies == aSpecies);
    }
  bool getIsInterface()
    {
      return isInterface;
    }
  bool getIsDeoligomerize()
    {
      return isDeoligomerize;
    }
  bool getIsOnMultiscale()
    {
      return isOnMultiscale;
    }
  bool getIsGaussianPopulation()
    {
      return isGaussianPopulation;
    }
  int getPopulatePriority()
    {
      return thePopulateProcess->getPriority();
    }
  void populateCompGaussian()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateGaussian(this);
        }
      else if(theMoleculeSize)
        {
          cout << "Warning: Species " << getIDString() <<
            " not populated." << std::endl;
        }
    }
  void addTaggedSpecies(Species* aSpecies)
    {
      isTag = true;
      //If one of the tagged species is off-lattice then
      //make the tag species off-lattice. The getPoint method
      //will convert the coord of lattice species to point
      //when logging:
      if(aSpecies->getIsOffLattice())
        {
          isOffLattice = true;
        }
      theTaggedSpeciesList.push_back(aSpecies);
    }
  void addTagSpecies(Species* aSpecies)
    {
      isTagged = true;
      theTagSpeciesList.push_back(aSpecies);
      aSpecies->addTaggedSpecies(this);
    }
  void setIsDeoligomerize(Species* aSpecies, const unsigned aSize)
    {
      isDeoligomerize = true;
      theBoundCnts.resize(aSize+1);
      theDeoligomerizedProduct = aSpecies;
    }
  void setIsPeriodic()
    {
      isPeriodic = true;
    }
  void setIsInterface()
    {
      isInterface = true;
    }
  void setIsRegularLattice(unsigned aDiffuseSize)
    {
      isRegularLattice = true;
      theDiffuseSize = aDiffuseSize;
    }
  void setWalkPropensity(const double aPropensity)
    {
      theWalkPropensity = aPropensity;
      if(aPropensity > 0)
        {
          isDiffusing = true;
        }
    }
  bool getIsRegularLattice()
    {
      return isRegularLattice;
    }
  bool getIsTagged()
    {
      return isTagged;
    }
  bool getIsTag()
    {
      return isTag;
    }
  bool getIsPopulateSpecies()
    {
      return (thePopulateProcess != NULL);
    }
  void populateCompUniformDense(unsigned* voxelIDs, unsigned* aCount)
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformDense(this, voxelIDs, aCount);
        }
      else if(theMoleculeSize)
        {
          cout << "Species:" << theVariable->getFullID().asString() <<
            " not populated." << std::endl;
        }
    }
  void populateCompUniformSparse()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformSparse(this);
        }
      else if(theMoleculeSize)
        {
          cout << "Species:" << theVariable->getFullID().asString() <<
            " not populated." << std::endl;
        }
    }
  void populateUniformOnDiffusiveVacant()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformOnDiffusiveVacant(this);
        }
      else if(theMoleculeSize)
        {
          cout << "Species:" << theVariable->getFullID().asString() <<
            " not populated." << std::endl;
        }
    }
  void populateUniformOnMultiscale()
    {
      if(thePopulateProcess)
        {
          thePopulateProcess->populateUniformOnMultiscale(this);
        }
      else if(theMoleculeSize)
        {
          cout << "Species:" << theVariable->getFullID().asString() <<
            " not populated." << std::endl;
        }
    }
  Variable* getVariable() const
    {
      return theVariable;
    }
  std::vector<unsigned> getSourceCoords()
    {
      std::vector<unsigned> aCoords;
      /*
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          std::vector<unsigned>& 
            aSourceCoords(theMolecules[i]->subunit->sourceCoords);
          for(unsigned j(0); j != aSourceCoords.size(); ++j)
            {
              if(aSourceCoords[j] != theNullCoord)
                {
                  aCoords.push_back(aSourceCoords[j]);
                }
            }
        }
        */
      return aCoords;
    }
  std::vector<unsigned> getTargetCoords()
    {
      std::vector<unsigned> aCoords;
      /*
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          std::vector<unsigned>& 
            aTargetCoords(theMolecules[i]->subunit->targetCoords);
          for(unsigned j(0); j != aTargetCoords.size(); ++j)
            {
              if(aTargetCoords[j] != theNullCoord)
                {
                  aCoords.push_back(aTargetCoords[j]);
                }
            }
        }
        */
      return aCoords;
    }
  std::vector<unsigned> getSharedCoords()
    {
      std::vector<unsigned> aCoords;
      /*
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          std::vector<unsigned>& 
            aSharedLipids(theMolecules[i]->subunit->sharedLipids);
          for(unsigned j(0); j != aSharedLipids.size(); ++j)
            {
              if(aSharedLipids[j] != theNullCoord)
                {
                  aCoords.push_back(aSharedLipids[j]);
                }
            }
        }
        */
      return aCoords;
    }
  unsigned size() const
    {
      return theMoleculeSize;
    }
  Voxel* getMolecule(const unsigned anIndex)
    {
      return theMolecules[anIndex];
    }
  Point getPoint(const unsigned anIndex) const
    {
      if(isOffLattice)
        {
          if(theMolecules[anIndex]->point)
            {
              return *theMolecules[anIndex]->point;
            }
          return theStepper->coord2point(getCoord(anIndex));
        }
      return theStepper->coord2point(getCoord(anIndex));
    }
  Point coord2point(const unsigned aCoord) const
    {
      if(theLattice[aCoord].point)
        {
          return *theLattice[aCoord].point;
        }
      return theStepper->coord2point(aCoord);
    }
  unsigned short getID() const
    {
      return theID;
    }
  //Call updateMolecules() before calling this function if
  //this species is a Tag (GFP) species, whose tag and molecules
  //are always not up to date.
  double getSquaredDisplacement(const unsigned index)
    {
      Point aCurrentPoint;
      if(isOrigins)
        {
          aCurrentPoint = getPeriodicPoint(index);
        }
      else
        {
          aCurrentPoint = theStepper->getPeriodicPoint(getCoord(index),
                                                       theDimension,
                                                   &theTags[index].origin);
        }
      const double nDisplacement(distance(theTags[index].origin.point,
                                          aCurrentPoint)*
                                 theStepper->getVoxelRadius()*2);
      return nDisplacement*nDisplacement;
    }
  double getMeanSquaredDisplacement()
    {
      //For GFP species, whose molecules are not always up to date:
      updateMolecules();
      if(!theMoleculeSize)
        {
          return 0;
        }
      double aTotalSquaredDisplacement(0);
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          aTotalSquaredDisplacement += getSquaredDisplacement(i);
        }
      return aTotalSquaredDisplacement/theMoleculeSize;
    }
  //Only works for 2D periodic diffusion now:
  Point getPeriodicPoint(const unsigned index)
    {
      Origin& anOrigin(theTags[index].origin);
      Point aPoint(getPoint(index));
      disp_(aPoint, theWidthVector, 
            anOrigin.row*lipRows*nDiffuseRadius*sqrt(3));
      disp_(aPoint, theLengthVector, anOrigin.col*lipCols*nDiffuseRadius*2);
      return aPoint;
    }
  void setIsSubunitInitialized()
    {
      isSubunitInitialized = true;
    }
  void setIsMultiscale()
    {
      isMultiscale = true;
      isMultiscaleComp = true;
    }
  bool getIsMultiscale()
    {
      return isMultiscale;
    }
  bool getIsMultiscaleComp()
    {
      return isMultiscaleComp;
    }
  void setIsCompVacant()
    {
      isCompVacant = true;
      isVacant = true;
      setVacantSpecies(this);
    }
  void setIsInContact()
    {
      isInContact = true;
    }
  void setIsCentered()
    {
      isCentered = true;
    }
  void setIsPopulated()
    {
      theInitCoordSize = theMoleculeSize;
      getVariable()->setValue(theMoleculeSize);
    }
  void shuffle()
    {
      std::random_shuffle(theMolecules.begin(), theMolecules.end());
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          theMolecules[i]->idx = i+theStride*theID;
        }
    }
  void initCollisionCnt()
    {
      theSpeciesCollisionCnt = 0;
      collisionCnts.resize(theMoleculeSize);
      for(std::vector<unsigned>::iterator 
          i(collisionCnts.begin()); i != collisionCnts.end(); ++i)
        {
          *i = 0;
        }
    }
  void finalizeSpecies()
    {
      for(unsigned i(0); i != theDiffusionInfluencedReactions.size(); ++i)
        {
          if(theDiffusionInfluencedReactions[i] &&
             theDiffusionInfluencedReactions[i]->getCollision())
            {
              initCollisionCnt();
              //For non-diffusing target species that does not have
              //DIRP list:
              theSpecies[i]->initCollisionCnt();
              break;
            }
        }
      //need to shuffle molecules of the compVacant species if it has
      //diffusing vacant species to avoid bias when random walking:
      if(isCompVacant)
        {
          for(unsigned i(0); i != theComp->species.size(); ++i)
            {
              if(theComp->species[i]->getIsDiffusiveVacant())
                {
                  shuffle();
                  break;
                }
            }
          //Is GFP tagged:
          if(isTagged)
            {
              for(unsigned i(0); i != theMoleculeSize; ++i)
                {
                  resetTagOrigin(i);
                }
            }
        }
      if(!isVacant && !getIsPopulated())
        {
          THROW_EXCEPTION(ValueError, getIDString() +
             ": has a non-zero value, theMoleculeSize:" + 
             int2str(theMoleculeSize) + ", variable value:" + 
             int2str(getVariable()->getValue()) +
             ", init size:" + int2str(theInitCoordSize) + 
             ", but is not populated.");
        }
      /*
      std::cout << getIDString() << std::endl;
      std::cout << "\tDiffusionInfluenced:" << std::endl;
      for(unsigned i(0); i != theDiffusionInfluencedReactions.size(); ++i)
        {
          if(theDiffusionInfluencedReactions[i])
            {
              std::cout << "\t\t" << theDiffusionInfluencedReactions[i
                ]->getIDString() << std::endl;
            }
        }
      std::cout << "\tAddMolecule:" << std::endl;
      for(unsigned i(0); i != theInterruptedProcessesAddMolecule.size(); ++i)
        {
          std::cout << "\t\t" << theInterruptedProcessesAddMolecule[i
            ]->getIDString() << std::endl;
        }
      std::cout << "\tRemoveMolecule:" << std::endl;
      for(unsigned i(0); i != theInterruptedProcessesRemoveMolecule.size(); ++i)
        {
          std::cout << "\t\t" << theInterruptedProcessesRemoveMolecule[i
            ]->getIDString() << std::endl;
        }
      std::cout << "\tEndDiffusion:" << std::endl;
      for(unsigned i(0); i != theInterruptedProcessesEndDiffusion.size(); ++i)
        {
          std::cout << "\t\t" << theInterruptedProcessesEndDiffusion[i
            ]->getIDString() << std::endl;
        }
        */
    }
  unsigned getCollisionCnt(unsigned anIndex)
    {
      return collisionCnts[anIndex];
    }
  unsigned getSpeciesCollisionCnt()
    {
      return theSpeciesCollisionCnt;
    }
  void setIsDiffusiveVacant()
    {
      isDiffusiveVacant = true;
      isVacant = true;
    }
  void setIsReactiveVacant()
    {
      isReactiveVacant = true;
      isVacant = true;
    }
  void setIsOffLattice()
    {
      isOffLattice = true;
    }
  void resetFixedAdjoins()
    {
      isFixedAdjoins = false;
      for(unsigned i(0); i != theComp->species.size(); ++i)
        {
          theComp->species[i]->setIsFixedAdjoins(false);
        }
    }
  void setIsFixedAdjoins(bool state)
    {
      isFixedAdjoins = state;
    }
  void setIsPolymer(std::vector<double> bendAngles, int aDirectionality)
    {
      theBendAngles.resize(0);
      thePolymerDirectionality = aDirectionality;
      isPolymer = true;
      for(std::vector<double>::const_iterator i(bendAngles.begin()); 
          i != bendAngles.end(); ++i)
        {
          theBendAngles.push_back(*i);
          if(thePolymerDirectionality != 0)
            {
              theBendAngles.push_back(*i+M_PI);
            }
        }
    }
  void setDiffusionCoefficient(double aCoefficient)
    {
      D = aCoefficient;
      if(D > 0)
        {
          isDiffusing = true;
        }
    }
  double getDiffusionCoefficient() const
    {
      return D;
    }
  double getWalkProbability() const
    {
      return theWalkProbability;
    }
  bool getIsPolymer() const
    {
      return isPolymer;
    }
  bool getIsOffLattice()
    {
      return isOffLattice;
    }
  bool getIsSubunitInitialized() const
    {
      return isSubunitInitialized;
    }
  bool getIsDiffusing() const
    {
      return isDiffusing;
    }
  bool getIsCompVacant() const
    {
      return isCompVacant;
    }
  bool getIsVacant() const
    {
      return isVacant;
    }
  bool getIsDiffusiveVacant()
    {
      return isDiffusiveVacant;
    }
  bool getIsReactiveVacant()
    {
      return isReactiveVacant;
    }
  bool getIsInContact() const
    {
      return isInContact;
    }
  bool getIsCentered() const
    {
      return isCentered;
    }
  bool getIsPopulated() const
    {
      return theMoleculeSize == theInitCoordSize;
    }
  double getDiffusionInterval() const
    {
      return theDiffusionInterval;
    }
  unsigned getDimension()
    {
      return theDimension;
    }
  void setDimension(unsigned aDimension)
    {
      theDimension = aDimension;
      if(theDimension == 3)
        {
          isFixedAdjoins = true;
        }
    }
  void resetFinalizeReactions()
    {
      for(unsigned i(0); i != isFinalizeReactions.size(); ++i)
        {
          isFinalizeReactions[i] = false;
        }
    }
  void finalizeReactions()
    {
      for(unsigned i(0); i != isFinalizeReactions.size(); ++i)
        {
          if(isFinalizeReactions[i])
            {
              theDiffusionInfluencedReactions[i]->finalizeReaction();
            }
        }
      for(unsigned i(0); i != theInterruptedProcessesEndDiffusion.size(); ++i)
        {
          theInterruptedProcessesEndDiffusion[i]->interruptedEndDiffusion(this);
        }
    }
  void addCollisionCnt()
    {
      ++theSpeciesCollisionCnt;
    }
  void addCollision(Voxel* aVoxel)
    {
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          if(aVoxel == theMolecules[i])
            {
              ++collisionCnts[i];
              return;
            }
        }
      cout << "error in species add collision" << std::endl;
    }
  unsigned getID(const Voxel* aVoxel) const
    {
      return aVoxel->idx/theStride;
    }
  unsigned getID(const Voxel& aVoxel) const
    {
      return aVoxel.idx/theStride;
    }
  void pushInterfaceConsts(std::vector<double>& interfaceConsts)
    {
      theInterfaceConsts.push_back(interfaceConsts);
    }
  void setUnityInterfaceConsts()
    {
      for(unsigned i(0); i != theInterfaceConsts.size(); ++i)
        {
          for(unsigned j(0); j != theInterfaceConsts[i].size(); ++j)
            {
              theInterfaceConsts[i][j] = 1;
            }
        }
    }
  double getInterfaceConst(Voxel* aVoxel, unsigned index)
    {
      return theInterfaceConsts[getIndex(aVoxel)][index];
    }
  void walk()
    {
      double ic(1);
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoiningCoordSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          if(!isFixedAdjoins)
            {
              size = source->diffuseSize;
            }
          Voxel* target(&theLattice[source->adjoiningCoords[
                        theRng.Integer(size)]]);
          if(getID(target) == theVacantID)
            {
              if(theWalkProbability == 1 || theRng.Fixed() < theWalkProbability)
                {
                  source->idx = target->idx;
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                }
            }
          else
            {
              if(theSpecies[getID(target)]->getIsInterface())
                {
                  //Some interface voxels do not have pointers to the
                  //off lattice subunits, so their adjoiningSize == diffuseSize:
                  if(target->adjoiningSize == target->diffuseSize)
                    {
                      continue;
                    }
                  unsigned coord(theRng.Integer(target->adjoiningSize-
                                                target->diffuseSize));
                  //For interfaceSpecies, theReactionProbabilities[tarID]
                  //contains the value k/DA, so 
                  //p = (k/DA)*ic
                  //  = (k/DA)*8*nv*rv*rv/(V*bf)
                  ic = theSpecies[getID(target)]->getInterfaceConst(target,
                                                                    coord);
                  coord = target->adjoiningCoords[coord+target->diffuseSize];
                  target = &theLattice[coord];
                }
              const unsigned tarID(getID(target));
              DiffusionInfluencedReactionProcess* 
                aReaction(theDiffusionInfluencedReactions[tarID]);
              if(aReaction)
                {
                  if(aReaction->getCollision() == 3)
                    { 
                      ++theSpeciesCollisionCnt;
                      Species* targetSpecies(theSpecies[tarID]);
                      targetSpecies->addCollisionCnt();
                      return;
                    }
                  //If it meets the reaction probability:
                  if(theReactionProbabilities[tarID] == 1 ||
                     theRng.Fixed() < theReactionProbabilities[tarID]*ic)
                    { 
                      if(aReaction->getCollision() && 
                         aReaction->getCollision() != 3)
                        { 
                          ++collisionCnts[i];
                          Species* targetSpecies(theSpecies[tarID]);
                          targetSpecies->addCollision(target);
                          if(aReaction->getCollision() != 2)
                            {
                              return;
                            }
                        }
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(source, target, i, target->idx%theStride, tarID);
                      //If the reaction is successful, the last molecule of this
                      //species will replace the pointer of i, so we need to 
                      //decrement i to perform the diffusion on it. However, if
                      //theMoleculeSize didn't decrease, that means the
                      //currently walked molecule was a product of this
                      //reaction and so we don't need to walk it again by
                      //decrementing i.
                      if(theMoleculeSize < aMoleculeSize)
                        {
                          --i;
                        }
                    }
                }
            }
        }
    }
  //Just perform DIRP at every diffusion step interval without walking
  void walkReact()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoiningCoordSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          if(!isFixedAdjoins)
            {
              size = source->diffuseSize;
            }
          Voxel* target(&theLattice[source->adjoiningCoords[
                        theRng.Integer(size)]]);
          if(theSpecies[getID(target)]->getIsInterface())
            {
              //Some interface voxels do not have pointers to the
              //off lattice subunits, so their adjoiningSize == diffuseSize:
              if(target->adjoiningSize == target->diffuseSize)
                {
                  continue;
                }
              unsigned coord(theRng.Integer(target->adjoiningSize-
                                            target->diffuseSize));
              coord = target->adjoiningCoords[coord+target->diffuseSize];
              target = &theLattice[coord];
            }
          const unsigned tarID(getID(target));
          DiffusionInfluencedReactionProcess* 
            aReaction(theDiffusionInfluencedReactions[tarID]);
          if(aReaction)
            {
              if(aReaction->getCollision() == 3)
                { 
                  ++theSpeciesCollisionCnt;
                  Species* targetSpecies(theSpecies[tarID]);
                  targetSpecies->addCollisionCnt();
                  return;
                }
              //If it meets the reaction probability:
              if(theReactionProbabilities[tarID] == 1 ||
                 theRng.Fixed() < theReactionProbabilities[tarID])
                { 
                  if(aReaction->getCollision() &&
                     aReaction->getCollision() != 3)
                    { 
                      ++collisionCnts[i];
                      Species* targetSpecies(theSpecies[tarID]);
                      targetSpecies->addCollision(target);
                      if(aReaction->getCollision() != 2)
                        {
                          return;
                        }
                    }
                  unsigned aMoleculeSize(theMoleculeSize);
                  react(source, target, i, target->idx%theStride, tarID);
                  //If the reaction is successful, the last molecule of this
                  //species will replace the pointer of i, so we need to 
                  //decrement i to perform the diffusion on it. However, if
                  //theMoleculeSize didn't decrease, that means the
                  //currently walked molecule was a product of this
                  //reaction and so we don't need to walk it again by
                  //decrementing i.
                  if(theMoleculeSize < aMoleculeSize)
                    {
                      i -= aMoleculeSize-theMoleculeSize;
                    }
                }
            }
        }
    }
  void walkVacant()
    {
      updateVacantMolecules();
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          int size;
          if(isFixedAdjoins)
            {
              size = theAdjoiningCoordSize;
            }
          else
            {
              size = source->diffuseSize;
            }
          Voxel* target(&theLattice[source->adjoiningCoords[
                        theRng.Integer(size)]]);
          if(getID(target) == theVacantID)
            {
              if(theWalkProbability == 1 || theRng.Fixed() < theWalkProbability)
                {
                  source->idx = target->idx;
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                }
            }
        }
    }
  void walkTrail()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoiningCoordSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          if(!isFixedAdjoins)
            {
              size = source->diffuseSize;
            }
          Voxel* target(&theLattice[source->adjoiningCoords[
                        theRng.Integer(size)]]);
          if(getID(target) == theVacantID)
            {
              if(theWalkProbability == 1 || theRng.Fixed() < theWalkProbability)
                {
                  theTrailSpecies->addMolecule(source);
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                }
            }
          else if(getID(target) == theTrailSpecies->getID())
            {
              if(theWalkProbability == 1 || theRng.Fixed() < theWalkProbability)
                {
                  //Can be optimized by directly replacing source with target
                  //voxel of theTrailSpecies:
                  theTrailSpecies->addMolecule(source);
                  theTrailSpecies->softRemoveMolecule(target);
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                }
            }
          else
            {
              if(theSpecies[getID(target)]->getIsInterface())
                {
                  //Some interface voxels do not have pointers to the
                  //off lattice subunits, so their adjoiningSize == diffuseSize:
                  if(target->adjoiningSize == target->diffuseSize)
                    {
                      continue;
                    }
                  unsigned coord(theRng.Integer(target->adjoiningSize-
                                                target->diffuseSize));
                  coord = target->adjoiningCoords[coord+target->diffuseSize];
                  target = &theLattice[coord];
                }
              const unsigned tarID(getID(target));
              DiffusionInfluencedReactionProcess* 
                aReaction(theDiffusionInfluencedReactions[tarID]);
              if(aReaction)
                {
                  //If it meets the reaction probability:
                  if(theReactionProbabilities[tarID] == 1 ||
                     theRng.Fixed() < theReactionProbabilities[tarID])
                    { 
                      if(aReaction->getCollision())
                        { 
                          ++collisionCnts[i];
                          Species* targetSpecies(theSpecies[tarID]);
                          targetSpecies->addCollision(target);
                          if(aReaction->getCollision() != 2)
                            {
                              return;
                            }
                        }
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(source, target, i, target->idx%theStride, tarID);
                      //If the reaction is successful, the last molecule of this
                      //species will replace the pointer of i, so we need to 
                      //decrement i to perform the diffusion on it. However, if
                      //theMoleculeSize didn't decrease, that means the
                      //currently walked molecule was a product of this
                      //reaction and so we don't need to walk it again by
                      //decrementing i.
                      if(theMoleculeSize < aMoleculeSize)
                        {
                          --i;
                        }
                    }
                }
            }
        }
    }
  void walkTrailNeighbor()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoiningCoordSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          if(!isFixedAdjoins)
            {
              size = source->diffuseSize;
            }
          Voxel* target(&theLattice[source->adjoiningCoords[
                        theRng.Integer(size)]]);
          if(getID(target) == theVacantID)
            {
              if(theWalkProbability == 1 || theRng.Fixed() < theWalkProbability)
                {
                  theTrailSpecies->addMolecule(source);
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                  for(unsigned j(target->diffuseSize); j != target->trailSize;
                      ++j)
                    {
                      Voxel* trail(&theLattice[target->adjoiningCoords[j]]);
                      if(getID(trail) == theVacantID)
                        {
                          theTrailSpecies->addMolecule(trail);
                        }
                    }
                }
            }
          else if(getID(target) == theTrailSpecies->getID())
            {
              if(theWalkProbability == 1 || theRng.Fixed() < theWalkProbability)
                {
                  //Can be optimized by directly replacing source with target
                  //voxel of theTrailSpecies:
                  theTrailSpecies->addMolecule(source);
                  theTrailSpecies->softRemoveMolecule(target);
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                  for(unsigned j(target->diffuseSize); j != target->trailSize;
                      ++j)
                    {
                      Voxel* trail(&theLattice[target->adjoiningCoords[j]]);
                      if(getID(trail) == theVacantID)
                        {
                          theTrailSpecies->addMolecule(trail);
                        }
                    }
                }
            }
          else
            {
              if(theSpecies[getID(target)]->getIsInterface())
                {
                  //Some interface voxels do not have pointers to the
                  //off lattice subunits, so their adjoiningSize == diffuseSize:
                  if(target->adjoiningSize == target->diffuseSize)
                    {
                      continue;
                    }
                  unsigned coord(theRng.Integer(target->adjoiningSize-
                                                target->diffuseSize));
                  coord = target->adjoiningCoords[coord+target->diffuseSize];
                  target = &theLattice[coord];
                }
              const unsigned tarID(getID(target));
              DiffusionInfluencedReactionProcess* 
                aReaction(theDiffusionInfluencedReactions[tarID]);
              if(aReaction)
                {
                  //If it meets the reaction probability:
                  if(theReactionProbabilities[tarID] == 1 ||
                     theRng.Fixed() < theReactionProbabilities[tarID])
                    { 
                      if(aReaction->getCollision())
                        { 
                          ++collisionCnts[i];
                          Species* targetSpecies(theSpecies[tarID]);
                          targetSpecies->addCollision(target);
                          if(aReaction->getCollision() != 2)
                            {
                              return;
                            }
                        }
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(source, target, i, target->idx%theStride, tarID);
                      //If the reaction is successful, the last molecule of this
                      //species will replace the pointer of i, so we need to 
                      //decrement i to perform the diffusion on it. However, if
                      //theMoleculeSize didn't decrease, that means the
                      //currently walked molecule was a product of this
                      //reaction and so we don't need to walk it again by
                      //decrementing i.
                      if(theMoleculeSize < aMoleculeSize)
                        {
                          --i;
                        }
                    }
                }
            }
        }
    }
  void walkMultiscalePropensity()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          Voxel* target(&theLattice[source->adjoiningCoords[
                        theRng.Integer(source->diffuseSize)]]);
          if(getID(target) == theVacantID)
            {
              if(!isIntersectMultiscale(source->coord, target->coord) &&
                 isMultiscaleWalkPropensity(source->coord, target->coord))
                {
                  removeMultiscaleMolecule(source, theTags[i].rotIndex);
                  addMultiscaleMolecule(target, i);
                  source->idx = target->idx;
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                }
            }
        }
    }
  void walkMultiscale()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          Voxel* target(&theLattice[source->adjoiningCoords[
                        theRng.Integer(source->diffuseSize)]]);
          if(getID(target) == theVacantID)
            {
              if(!isIntersectMultiscale(source->coord, target->coord))
                {
                  removeMultiscaleMolecule(source, theTags[i].rotIndex);
                  addMultiscaleMolecule(target, i);
                  source->idx = target->idx;
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                }
            }
        }
    }
  //This is used for perfect flat surface compartments only, so that we can
  //get correct 2D diffusion behavior. Also much faster than walk() when
  //performing periodic boundary diffusion (eg. walkRegular:5s vs walk:8.7s).
  void walkRegular()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoiningCoordSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned tarIndex(theRng.Integer(theDiffuseSize)); 
          const unsigned row((source->coord-lipStartCoord)/lipCols);
          int coordA(source->coord-lipStartCoord+
                     theAdjoinOffsets[row%2][tarIndex]);
          if(!isInLattice(coordA, theRowOffsets[tarIndex]+row))
            {
              continue;
            }
          Voxel* target(&theLattice[coordA+lipStartCoord]);
          if(getID(target) == theVacantID)
            {
              if(theWalkProbability == 1 || theRng.Fixed() < theWalkProbability)
                {
                  source->idx = target->idx;
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                }
            }
          else
            {
              if(theSpecies[getID(target)]->getIsInterface())
                {
                  //Some interface voxels do not have pointers to the
                  //off lattice subunits, so their adjoiningSize == diffuseSize:
                  if(target->adjoiningSize == target->diffuseSize)
                    {
                      continue;
                    }
                  unsigned coord(theRng.Integer(target->adjoiningSize-
                                                target->diffuseSize));
                  coord = target->adjoiningCoords[coord+target->diffuseSize];
                  target = &theLattice[coord];
                }
              const unsigned tarID(getID(target));
              DiffusionInfluencedReactionProcess* 
                aReaction(theDiffusionInfluencedReactions[tarID]);
              if(aReaction)
                {
                  if(aReaction->getCollision() == 3)
                    { 
                      ++theSpeciesCollisionCnt;
                      Species* targetSpecies(theSpecies[tarID]);
                      targetSpecies->addCollisionCnt();
                      return;
                    }
                  //If it meets the reaction probability:
                  if(theReactionProbabilities[tarID] == 1 ||
                     theRng.Fixed() < theReactionProbabilities[tarID])
                    { 
                      if(aReaction->getCollision() &&
                         aReaction->getCollision() != 3)
                        { 
                          ++collisionCnts[i];
                          Species* targetSpecies(theSpecies[tarID]);
                          targetSpecies->addCollision(target);
                          if(aReaction->getCollision() != 2)
                            {
                              return;
                            }
                        }
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(source, target, i, target->idx%theStride, tarID);
                      //If the reaction is successful, the last molecule of this
                      //species will replace the pointer of i, so we need to 
                      //decrement i to perform the diffusion on it. However, if
                      //theMoleculeSize didn't decrease, that means the
                      //currently walked molecule was a product of this
                      //reaction and so we don't need to walk it again by
                      //decrementing i.
                      if(theMoleculeSize < aMoleculeSize)
                        {
                          --i;
                        }
                    }
                }
            }
        }
    }
  void walkReactRegular()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoiningCoordSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned tarIndex(theRng.Integer(theDiffuseSize)); 
          const unsigned row((source->coord-lipStartCoord)/lipCols);
          int coordA(source->coord-lipStartCoord+
                     theAdjoinOffsets[row%2][tarIndex]);
          if(!isInLattice(coordA, theRowOffsets[tarIndex]+row))
            {
              continue;
            }
          Voxel* target(&theLattice[coordA+lipStartCoord]);
          if(theSpecies[getID(target)]->getIsInterface())
            {
              //Some interface voxels do not have pointers to the
              //off lattice subunits, so their adjoiningSize == diffuseSize:
              if(target->adjoiningSize == target->diffuseSize)
                {
                  continue;
                }
              unsigned coord(theRng.Integer(target->adjoiningSize-
                                            target->diffuseSize));
              coord = target->adjoiningCoords[coord+target->diffuseSize];
              target = &theLattice[coord];
            }
          const unsigned tarID(getID(target));

          DiffusionInfluencedReactionProcess* 
            aReaction(theDiffusionInfluencedReactions[tarID]);
          if(aReaction)
            {
              if(aReaction->getCollision() == 3)
                { 
                  ++theSpeciesCollisionCnt;
                  Species* targetSpecies(theSpecies[tarID]);
                  targetSpecies->addCollisionCnt();
                  return;
                }
              //If it meets the reaction probability:
              if(theReactionProbabilities[tarID] == 1 ||
                 theRng.Fixed() < theReactionProbabilities[tarID])
                { 
                  if(aReaction->getCollision() &&
                     aReaction->getCollision() != 3)
                    { 
                      ++collisionCnts[i];
                      Species* targetSpecies(theSpecies[tarID]);
                      targetSpecies->addCollision(target);
                      if(aReaction->getCollision() != 2)
                        {
                          return;
                        }
                    }
                  unsigned aMoleculeSize(theMoleculeSize);
                  react(source, target, i, target->idx%theStride, tarID);
                  //If the reaction is successful, the last molecule of this
                  //species will replace the pointer of i, so we need to 
                  //decrement i to perform the diffusion on it. However, if
                  //theMoleculeSize didn't decrease, that means the
                  //currently walked molecule was a product of this
                  //reaction and so we don't need to walk it again by
                  //decrementing i.
                  if(theMoleculeSize < aMoleculeSize)
                    {
                      --i;
                    }
                }
            }
        }
    }
  void walkReactRegularOrigins()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoiningCoordSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned tarIndex(theRng.Integer(theDiffuseSize)); 
          const unsigned row((source->coord-lipStartCoord)/lipCols);
          int coordA(source->coord-lipStartCoord+
                     theAdjoinOffsets[row%2][tarIndex]);
          Origin anOrigin(theTags[i].origin);
          if(!isInLatticeOrigins(coordA, theRowOffsets[tarIndex]+row,
                                 anOrigin))
            {
              continue;
            }
          Voxel* target(&theLattice[coordA+lipStartCoord]);
          if(theSpecies[getID(target)]->getIsInterface())
            {
              //Some interface voxels do not have pointers to the
              //off lattice subunits, so their adjoiningSize == diffuseSize:
              if(target->adjoiningSize == target->diffuseSize)
                {
                  continue;
                }
              unsigned coord(theRng.Integer(target->adjoiningSize-
                                            target->diffuseSize));
              coord = target->adjoiningCoords[coord+target->diffuseSize];
              target = &theLattice[coord];
            }
          const unsigned tarID(getID(target));
          DiffusionInfluencedReactionProcess* 
            aReaction(theDiffusionInfluencedReactions[tarID]);
          if(aReaction)
            {
              if(aReaction->getCollision() == 3)
                { 
                  ++theSpeciesCollisionCnt;
                  Species* targetSpecies(theSpecies[tarID]);
                  targetSpecies->addCollisionCnt();
                  return;
                }
              //If it meets the reaction probability:
              if(theReactionProbabilities[tarID] == 1 ||
                 theRng.Fixed() < theReactionProbabilities[tarID])
                { 
                  if(aReaction->getCollision() &&
                     aReaction->getCollision() != 3)
                    { 
                      ++collisionCnts[i];
                      Species* targetSpecies(theSpecies[tarID]);
                      targetSpecies->addCollision(target);
                      if(aReaction->getCollision() != 2)
                        {
                          return;
                        }
                    }
                  unsigned aMoleculeSize(theMoleculeSize);
                  react(source, target, i, target->idx%theStride, tarID);
                  //If the reaction is successful, the last molecule of this
                  //species will replace the pointer of i, so we need to 
                  //decrement i to perform the diffusion on it. However, if
                  //theMoleculeSize didn't decrease, that means the
                  //currently walked molecule was a product of this
                  //reaction and so we don't need to walk it again by
                  //decrementing i.
                  if(theMoleculeSize < aMoleculeSize)
                    {
                      --i;
                    }
                }
            }
        }
    }
  void walkRegularOrigins()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      unsigned size(theAdjoiningCoordSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned tarIndex(theRng.Integer(theDiffuseSize)); 
          const unsigned row((source->coord-lipStartCoord)/lipCols);
          int coordA(source->coord-lipStartCoord+
                     theAdjoinOffsets[row%2][tarIndex]);
          Origin anOrigin(theTags[i].origin);
          if(!isInLatticeOrigins(coordA, theRowOffsets[tarIndex]+row,
                                 anOrigin))
            {
              continue;
            }
          Voxel* target(&theLattice[coordA+lipStartCoord]);
          if(getID(target) == theVacantID)
            {
              if(theWalkProbability == 1 || theRng.Fixed() < theWalkProbability)
                {
                  source->idx = target->idx;
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                  theTags[i].origin = anOrigin;
                }
            }
          else
            {
              if(theSpecies[getID(target)]->getIsInterface())
                {
                  //Some interface voxels do not have pointers to the
                  //off lattice subunits, so their adjoiningSize == diffuseSize:
                  if(target->adjoiningSize == target->diffuseSize)
                    {
                      continue;
                    }
                  unsigned coord(theRng.Integer(target->adjoiningSize-
                                                target->diffuseSize));
                  coord = target->adjoiningCoords[coord+target->diffuseSize];
                  target = &theLattice[coord];
                }
              const unsigned tarID(getID(target));
              DiffusionInfluencedReactionProcess* 
                aReaction(theDiffusionInfluencedReactions[tarID]);
              if(aReaction)
                {
                  if(aReaction->getCollision() == 3)
                    { 
                      ++theSpeciesCollisionCnt;
                      Species* targetSpecies(theSpecies[tarID]);
                      targetSpecies->addCollisionCnt();
                      return;
                    }
                  //If it meets the reaction probability:
                  if(theReactionProbabilities[tarID] == 1 ||
                     theRng.Fixed() < theReactionProbabilities[tarID])
                    { 
                      if(aReaction->getCollision() &&
                         aReaction->getCollision() != 3)
                        { 
                          ++collisionCnts[i];
                          Species* targetSpecies(theSpecies[tarID]);
                          targetSpecies->addCollision(target);
                          if(aReaction->getCollision() != 2)
                            {
                              return;
                            }
                        }
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(source, target, i, target->idx%theStride, tarID);
                      //If the reaction is successful, the last molecule of this
                      //species will replace the pointer of i, so we need to 
                      //decrement i to perform the diffusion on it. However, if
                      //theMoleculeSize didn't decrease, that means the
                      //currently walked molecule was a product of this
                      //reaction and so we don't need to walk it again by
                      //decrementing i.
                      if(theMoleculeSize < aMoleculeSize)
                        {
                          --i;
                        }
                    }
                }
            }
        }
    }
  void walkOnMultiscaleRegular()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned tarIndex(theRng.Integer(theDiffuseSize)); 
          const unsigned row((source->coord-lipStartCoord)/lipCols);
          int coordA(source->coord-lipStartCoord+
                     theAdjoinOffsets[row%2][tarIndex]);
          if(!isInLattice(coordA, theRowOffsets[tarIndex]+row))
            {
              continue;
            }
          Voxel* target(&theLattice[coordA+lipStartCoord]);
          if(getID(target) == theVacantID)
            {
              if(theWalkProbability == 1 || theRng.Fixed() < theWalkProbability)
                {
                  source->idx = theTags[i].multiIdx;
                  if(theTags[i].multiIdx != target->idx)
                    {
                      theTags[i].multiIdx = target->idx;
                    }
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                }
            }
          else
            {
              const unsigned tarID(getID(target));
              if(theDiffusionInfluencedReactions[tarID])
                {
                  //If it meets the reaction probability:
                  if(theReactionProbabilities[tarID] == 1 ||
                     theRng.Fixed() < theReactionProbabilities[tarID])
                    {
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(source, target, i, target->idx%theStride, tarID);
                      //If the reaction is successful, the last molecule of this
                      //species will replace the pointer of i, so we need to 
                      //decrement i to perform the diffusion on it. However, if
                      //theMoleculeSize didn't decrease, that means the
                      //currently walked molecule was a product of this
                      //reaction and so we don't need to walk it again by
                      //decrementing i.
                      if(theMoleculeSize < aMoleculeSize)
                        {
                          --i;
                        }
                    }
                }
            }
        }
    }
  void rotateMultiscaleRegular()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const int tarIndex(theRng.Integer(2)); 
          const unsigned coordA(source->coord-vacStartCoord);
          const int rowA(coordA/lipCols);
          if(!isIntersectMultiscaleRegular(i, coordA, rowA,
                         theRotOffsets[rowA%2][theTags[i].rotIndex][tarIndex]))
            { 
              unsigned srcIndex(0);
              unsigned srcRotIndex(theTags[i].rotIndex+1);
              if(!tarIndex)
                {
                  srcIndex = 1;
                  if(!theTags[i].rotIndex)
                    {
                      srcRotIndex = theRotateSize-1;
                    }
                  else
                    {
                      srcRotIndex = theTags[i].rotIndex-1; 
                    }
                }
              else if(srcRotIndex == theRotateSize)
                {
                  srcRotIndex = 0;
                }
              moveMultiscaleMoleculeRegular(coordA, rowA, 
                    theRotOffsets[rowA%2][theTags[i].rotIndex][tarIndex],
                    theRotOffsets[rowA%2][srcRotIndex][srcIndex], i);
              theTags[i].rotIndex = srcRotIndex;
            }
        }
    }
  void rotateMultiscalePropensityRegular()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const int tarIndex(theRng.Integer(2)); 
          const unsigned coordA(source->coord-vacStartCoord);
          const int rowA(coordA/lipCols); 
          if(!isIntersectMultiscaleRegular(i, coordA, rowA,
                 theRotOffsets[rowA%2][theTags[i].rotIndex][tarIndex]))
            {
              unsigned srcIndex(0);
              unsigned srcRotIndex(theTags[i].rotIndex+1);
              if(!tarIndex)
                {
                  srcIndex = 1;
                  if(!theTags[i].rotIndex)
                    {
                      srcRotIndex = theRotateSize-1;
                    }
                  else
                    {
                      srcRotIndex = theTags[i].rotIndex-1; 
                    }
                }
              else if(srcRotIndex == theRotateSize)
                {
                  srcRotIndex = 0;
                }
              if(isMultiscaleWalkPropensityRegular(coordA, rowA,
                     theRotOffsets[rowA%2][theTags[i].rotIndex][tarIndex],
                     theRotOffsets[rowA%2][srcRotIndex][srcIndex]))
                {
                  moveMultiscaleMoleculeRegular(coordA, rowA, 
                        theRotOffsets[rowA%2][theTags[i].rotIndex][tarIndex],
                        theRotOffsets[rowA%2][srcRotIndex][srcIndex], i);
                  theTags[i].rotIndex = srcRotIndex;
                }
            }
        }
    }
  void walkMultiscalePropensityRegular()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned srcCoord(source->coord-vacStartCoord);
          const unsigned tarIndex(theRng.Integer(theDiffuseSize)); 
          const unsigned row(srcCoord/lipCols);
          int tarCoord(srcCoord+theAdjoinOffsets[row%2][tarIndex]);
          if(!isInLattice(tarCoord, theRowOffsets[tarIndex]+row))
            {
              continue;
            }
          Voxel* target(&theLattice[tarCoord+vacStartCoord]);
          if(getID(target) == theVacantID)
            {
              if(isMoveMultiscaleWalkPropensityRegular(i, srcCoord, row,
                     theTarOffsets[row%2][theTags[i].rotIndex][tarIndex],
                     theSrcOffsets[row%2][theTags[i].rotIndex][tarIndex], i))
                {
                  source->idx = target->idx;
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                }
            }
        }
    }
  void walkMultiscalePropensityRegularOrigins()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned srcCoord(source->coord-vacStartCoord);
          const unsigned tarIndex(theRng.Integer(theDiffuseSize)); 
          const unsigned row(srcCoord/lipCols);
          int tarCoord(srcCoord+theAdjoinOffsets[row%2][tarIndex]);
          Origin anOrigin(theTags[i].origin);
          if(!isInLatticeOrigins(tarCoord, theRowOffsets[tarIndex]+row,
                                 anOrigin))
            {
              continue;
            }
          Voxel* target(&theLattice[tarCoord+vacStartCoord]);
          if(getID(target) == theVacantID)
            {
              if(isMoveMultiscaleWalkPropensityRegular(i, srcCoord, row,
                     theTarOffsets[row%2][theTags[i].rotIndex][tarIndex],
                     theSrcOffsets[row%2][theTags[i].rotIndex][tarIndex], i))
                {
                  source->idx = target->idx;
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                  theTags[i].origin = anOrigin;
                }
            }
        }
    }
  void walkMultiscaleRegular()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned srcCoord(source->coord-vacStartCoord);
          const unsigned tarIndex(theRng.Integer(theDiffuseSize)); 
          const unsigned row(srcCoord/lipCols);
          int tarCoord(srcCoord+theAdjoinOffsets[row%2][tarIndex]);
          if(!isInLattice(tarCoord, theRowOffsets[tarIndex]+row))
            {
              continue;
            }
          Voxel* target(&theLattice[tarCoord+vacStartCoord]);
          if(!isIntersectMultiscaleRegular(i, srcCoord, row,
                   theTarOffsets[row%2][theTags[i].rotIndex][tarIndex]))
            {
              moveMultiscaleMoleculeRegular(srcCoord, row, 
                 theTarOffsets[row%2][theTags[i].rotIndex][tarIndex],
                 theSrcOffsets[row%2][theTags[i].rotIndex][tarIndex], i);
              source->idx = target->idx;
              target->idx = i+theStride*theID;
              theMolecules[i] = target;
            }
        }
    }
  void walkMultiscaleRegularOrigins()
    {
      const unsigned beginMoleculeSize(theMoleculeSize);
      for(unsigned i(0); i < beginMoleculeSize && i < theMoleculeSize; ++i)
        {
          Voxel* source(theMolecules[i]);
          const unsigned srcCoord(source->coord-vacStartCoord);
          const unsigned tarIndex(theRng.Integer(theDiffuseSize)); 
          const unsigned row(srcCoord/lipCols);
          int tarCoord(srcCoord+theAdjoinOffsets[row%2][tarIndex]);
          Origin anOrigin(theTags[i].origin);
          if(!isInLatticeOrigins(tarCoord, theRowOffsets[tarIndex]+row,
                                 anOrigin))
            {
              continue;
            }
          Voxel* target(&theLattice[tarCoord+vacStartCoord]);
          if(getID(target) == theVacantID)
            {
              if(!isIntersectMultiscaleRegular(i, srcCoord, row,
                   theTarOffsets[row%2][theTags[i].rotIndex][tarIndex]))
                {
                  moveMultiscaleMoleculeRegular(srcCoord, row, 
                     theTarOffsets[row%2][theTags[i].rotIndex][tarIndex],
                     theSrcOffsets[row%2][theTags[i].rotIndex][tarIndex], i);
                  source->idx = target->idx;
                  target->idx = i+theStride*theID;
                  theMolecules[i] = target;
                  theTags[i].origin = anOrigin;
                }
            }
        }
    }
  void react(Voxel* src, Voxel* tar, const unsigned srcIndex,
             const unsigned tarIndex, const unsigned tarID)
    {
      DiffusionInfluencedReactionProcess* aReaction(
                      theDiffusionInfluencedReactions[tarID]);
      if(aReaction->getA() == this)
        { 
          aReaction->react(src, tar, srcIndex, tarIndex);
        }
      else
        {
          aReaction->react(tar, src, tarIndex, srcIndex);
        }
      isFinalizeReactions[tarID] = true;
    }
  void setComp(Comp* aComp)
    {
      theComp = aComp;
    }
  Comp* getComp()
    {
      return theComp;
    }
  void setVariable(Variable* aVariable)
    {
      theVariable = aVariable;
    }
  unsigned getCoord(const unsigned anIndex) const
    {
      return theMolecules[anIndex]->coord;
    }
  void removeSurfacePlanes()
    {
      int newCoordSize(0);
      for(unsigned i(0); i < theMoleculeSize; ++i) 
        {
          unsigned aCoord(getCoord(i));
          if(theStepper->isRemovableEdgeCoord(aCoord, theComp))
            {
              Comp* aSuperComp(
                 theStepper->system2Comp(theComp->system->getSuperSystem())); 
              aSuperComp->vacantSpecies->addCompVoxel(aCoord);
            }
          else 
            { 
              theMolecules[newCoordSize] = &theLattice[aCoord];
              ++newCoordSize; 
            }
        }
      theMoleculeSize = newCoordSize;
      //Must resize, otherwise compVoxelSize will be inaccurate:
      theMolecules.resize(theMoleculeSize);
      theVariable->setValue(theMoleculeSize);
    }
  void removePeriodicEdgeVoxels()
    {
      int newCoordSize(0);
      for(unsigned i(0); i < theMoleculeSize; ++i) 
        {
          unsigned aCoord(getCoord(i));
          if(theStepper->isPeriodicEdgeCoord(aCoord, theComp))
            {
              theLattice[aCoord].idx = theLattice[theNullCoord].idx;
            }
          else 
            { 
              theMolecules[newCoordSize] = &theLattice[aCoord];
              ++newCoordSize; 
            }
        }
      theMoleculeSize = newCoordSize;
      //Must resize, otherwise compVoxelSize will be inaccurate:
      theMolecules.resize(theMoleculeSize);
      theVariable->setValue(theMoleculeSize);
    }
  void updateSpecies()
    {
      if(isCompVacant && (isDiffusiveVacant || isReactiveVacant))
        {
          theCompVoxels = new std::vector<Voxel*>;
          for(unsigned i(0); i != theMoleculeSize; ++i)
            { 
              theCompVoxels->push_back(theMolecules[i]);
            }
        }
    }
  //If it isReactiveVacant it will only be called by SNRP when it is substrate
  //If it isDiffusiveVacant it will only be called by DiffusionProcess before
  //being diffused. So we need to only check if it isVacant:
  void updateMolecules()
    {
      if(isDiffusiveVacant || isReactiveVacant)
        {
          updateVacantMolecules();
        }
      else if(isTag)
        {
          updateTagMolecules();
        }
    }
  //If it isReactiveVacant it will only be called by SNRP when it is substrate:
  void updateMoleculeSize()
    {
      if(isDiffusiveVacant || isReactiveVacant)
        {
          updateVacantCoordSize();
        }
    }
  void allocateTags(const unsigned size)
    {
      theMoleculeSize = size;
      theMolecules.resize(size);
      theTags.resize(size);
      updateTagMolecules();
    }
  void updateTagMolecules()
    {
      for(unsigned i(0); i != theTaggedSpeciesList.size(); ++i)
        {
          Species* aSpecies(theTaggedSpeciesList[i]);
          for(unsigned j(0); j != aSpecies->size(); ++j)
            {
              Tag& aTag(aSpecies->getTag(j));
              if(aTag.speciesID == theID)
                {
                  theMolecules[aTag.molID] = aSpecies->getMolecule(j);
                  theTags[aTag.molID] = aTag;
                }
            }
        }
      /* Debug
      unsigned cnt(0);
      unsigned size(0);
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          theMolecules[i] = NULL;
        }
      for(unsigned i(0); i != theTaggedSpeciesList.size(); ++i)
        {
          Species* aSpecies(theTaggedSpeciesList[i]);
          size += aSpecies->size();
          for(unsigned j(0); j != aSpecies->size(); ++j)
            {
              Tag& aTag(aSpecies->getTag(j));
              if(aTag.speciesID == theID)
                {
                  ++cnt;
                  if(theMolecules[aTag.molID] != NULL)
                    {
                      std::cout << "not null" << std::endl;
                    }
                  theMolecules[aTag.molID] = aSpecies->getMolecule(j);
                  theTags[aTag.molID] = aTag;
                }
              else
                {
                  std::cout << "not tagged:" << aSpecies->getIDString() << std::endl;
                }
            }
        }
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          if(theMolecules[i] == NULL)
            {
              std::cout << "nul mol" << std::endl;
            }
        }
      if(cnt != theMoleculeSize || size != cnt)
        {
          std::cout << getIDString() << std::endl;
          std::cout << "cnt:" << cnt << " msize:" << theMoleculeSize << " size:" << size << std::endl;
        }
        */
    }
  void updateMoleculeList()
    {
      unsigned aSize(0);
      for(unsigned i(0); i != theMoleculeSize; ++i)
        { 
          if(getID(theMolecules[i]) == theID)
            {
              theMolecules[aSize++] = theMolecules[i];
            }
        }
      theMoleculeSize = aSize;
      theVariable->setValue(aSize);
    }
  //Even if it is a isCompVacant, this method will be called by
  //VisualizationLogProcess (either Reactive or Diffusive),
  //or SNRP if it is Reactive, or DiffusionProcess if it is Diffusive:
  void updateVacantMolecules()
    {
      theMoleculeSize = 0;
      int aSize(theVacantSpecies->compVoxelSize());
      for(int i(0); i != aSize; ++i)
        { 
          Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
          if(getID(aVoxel) == theID)
            {
              ++theMoleculeSize;
              if(theMoleculeSize > theMolecules.size())
                {
                  theMolecules.push_back(aVoxel);
                }
              else
                {
                  theMolecules[theMoleculeSize-1] = aVoxel;
                }
            }
        }
      theVariable->setValue(theMoleculeSize);
    }
  void updateVacantCoordSize()
    {
      theMoleculeSize = 0;
      int aSize(theVacantSpecies->compVoxelSize());
      for(int i(0); i != aSize; ++i)
        { 
          Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
          if(getID(aVoxel) == theID)
            {
              ++theMoleculeSize;
            }
        }
      if(theMoleculeSize > theMolecules.size())
        {
          theMolecules.resize(theMoleculeSize);
        }
      theVariable->setValue(theMoleculeSize);
    }
  Tag& getTag(const unsigned anIndex)
    {
      if(isVacant)
        {
          return theNullTag;
        }
      return theTags[anIndex];
    }
  Tag& getTag(Voxel* aVoxel)
    {
      return getTag(getIndex(aVoxel));
    }
  Species* getMultiscaleVacantSpecies()
    {
      return theMultiscaleVacantSpecies;
    }
  void getInMultiCnts(const unsigned anIndex, std::vector<unsigned>& cnts)
    {
      std::vector<unsigned> coords;
      getMultiCoords(anIndex, coords);
      for(unsigned i(0); i != coords.size(); ++i)
        {
          unsigned id(theLattice[coords[i]].idx/theStride);
          if(id != getID())
            {
              ++cnts[id];
            }
        }
    }
  void addMolecule(Voxel* aVoxel)
    {
      addMolecule(aVoxel, theNullTag);
    }
  void addMoleculeExMulti(Voxel* aVoxel, const unsigned boundCnt)
    {
      //use direct since we don't want to count bounds:
      addMoleculeDirect(aVoxel);
      if(isDeoligomerize)
        {
          theTags[theMoleculeSize-1].boundCnt = boundCnt;
          theBoundCnts[boundCnt]++;
        }
    }
  void addMoleculeInMulti(Voxel* aVoxel, const unsigned multiIdx,
                          const unsigned boundCnt)
    {
      //use direct since we don't want to count bounds:
      addMoleculeDirect(aVoxel);
      theTags[theMoleculeSize-1].multiIdx = multiIdx;
      if(isDeoligomerize)
        {
          theTags[theMoleculeSize-1].boundCnt = boundCnt;
          theBoundCnts[boundCnt]++;
        }
    }
  void addMoleculeInMulti(Voxel* aVoxel, const unsigned multiIdx)
    {
      addMoleculeTagless(aVoxel);
      theTags[theMoleculeSize-1].multiIdx = multiIdx;
    }
  void softReplaceMolecule(const unsigned index, Voxel* aVoxel, const Tag& aTag,
                           Species* aSpecies)
    {
      theMolecules[index] = aVoxel;
      theTags[index] = aTag;
      aVoxel->idx = index+theStride*theID;
      /*
      if(isDeoligomerize && !aSpecies->getIsDeoligomerize())
        {
          addBounds(aVoxel);
        }
      else if(!isDeoligomerize && aSpecies->getIsDeoligomerize())
        {
          removeAdjoinBounds(aVoxel);
        }
        */
    }
  void addMoleculeDirect(Voxel* aVoxel)
    {
      aVoxel->idx = theMoleculeSize+theStride*theID;
      softAddMolecule(aVoxel);
    }
  void softAddMolecule(Voxel* aVoxel)
    {
      ++theMoleculeSize; 
      if(theMoleculeSize > theMolecules.size())
        {
          theMolecules.resize(theMoleculeSize);
          theTags.push_back(theNullTag);
        }
      else
        {
          theTags[theMoleculeSize-1] = theNullTag;
        }
      theMolecules[theMoleculeSize-1] = aVoxel;
      theVariable->setValue(theMoleculeSize);
    }
  void addMoleculeTagless(Voxel* aVoxel)
    {
      addMoleculeDirect(aVoxel);
      if(isDeoligomerize)
        {
          addBounds(aVoxel);
        }
      //In future do this for even vacant species:
      //Only call addMolecule or removeMolecule when you know that
      //at the end of the reaction that molecule will not be overwritten
      //by other species.
      for(unsigned i(0); i != theInterruptedProcessesAddMolecule.size(); ++i)
        {
          theInterruptedProcessesAddMolecule[i
            ]->interruptedAddMolecule(this, theMoleculeSize-1);
        }
    }
  void addMolecule(Voxel* aVoxel, Tag& aTag)
    {
      if(!isVacant)
        {
          if(isMultiscale)
            {
              addMultiscaleMolecule(aVoxel, theMoleculeSize);
            }
          addMoleculeTagged(aVoxel, aTag);
        }
      else
        {
          aVoxel->idx = theID*theStride;
        }
    }
  void setIsTagged()
    {
      //Is GFP tagged:
      isTagged = true;
    }
  void addMoleculeTagged(Voxel* aVoxel, Tag& aTag)
    {
      addMoleculeTagless(aVoxel);
      //Is GFP tagged:
      if(isTagged)
        {
          //If it is theNullTag, the molecule of this species is being tagged
          //for the first time:
          if(aTag.speciesID == theNullID)
            {
              resetTagOrigin(theMoleculeSize-1);
            }
          //Previous tag exists in another species and is being transferred to
          //the molecule of this species:
          else
            {
              theTags[theMoleculeSize-1].origin = aTag.origin;
              theTags[theMoleculeSize-1].speciesID = aTag.speciesID;
              theTags[theMoleculeSize-1].molID = aTag.molID;
            }
        }
    }
  bool isIntersectMultiscaleRegular(unsigned& srcIndex,
                                    const unsigned coordA, 
                                    const unsigned rowA,
                                    const std::vector<int>& anOffsets)
    {
      bool isIntersect(false);
      for(unsigned i(0); i != anOffsets.size(); ++i)
        {
          const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coordB(coordA+anOffsets[i]);
          if(isInLattice(coordB, offsetRow))
            {
              unsigned idx(theLattice[coordB+lipStartCoord].idx);
              const unsigned tarID(idx/theStride);
              if(!isMultiMultiReactive &&
                 (theSpecies[tarID]->getIsMultiscale() || 
                  theSpecies[tarID]->getIsOnMultiscale()))
                {
                  return true;
                }
              else if(isMultiMultiReactantID[tarID])
                {
                  //If it meets the reaction probability:
                  if(theReactionProbabilities[tarID] == 1 ||
                     theRng.Fixed() < theReactionProbabilities[tarID])
                    { 
                      Voxel* target(&theLattice[coordB+lipStartCoord]);
                      if(!theSpecies[tarID]->getIsMultiscale())
                        {
                          idx = theSpecies[tarID
                            ]->getTag(idx%theStride).multiIdx;
                        }
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(target, target, srcIndex, idx%theStride,
                            idx/theStride);
                      if(theMoleculeSize < aMoleculeSize)
                        {
                          --srcIndex;
                        }
                      return true;
                    }
                  isIntersect = true;
                }
              else if(theSpecies[tarID]->getIsMultiscale() ||
                      theSpecies[tarID]->getIsOnMultiscale())
                {
                  isIntersect = true;
                }
            }
        }
      return isIntersect;
    }
  bool isMoveMultiscaleWalkPropensityRegular(unsigned& srcIndex,
                                             const unsigned coordA, 
                                             const unsigned rowA,
                                             const std::vector<int>& tarOffsets,
                                             const std::vector<int>& srcOffsets,
                                             const unsigned index)
    {
      bool isIntersect(false);
      theIdxList.resize(0);
      theBindList.resize(0);
      theVacantList.resize(0);
      theUnbindList.resize(0);
      int tarCnt(0);
      int srcCnt(0);
      //count tar
      for(unsigned i(0); i != tarOffsets.size(); ++i)
        {
          const int offsetRow((tarOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coordB(coordA+tarOffsets[i]);
          if(isInLattice(coordB, offsetRow))
            {
              const unsigned coord(coordB+lipStartCoord);
              unsigned idx(theLattice[coord].idx);
              const unsigned tarID(idx/theStride);
              //If it doesn't perform MultiMulti reaction and
              //intersects with another Multi:
              if(!isMultiMultiReactive &&
                 (theSpecies[tarID]->getIsMultiscale() ||
                  theSpecies[tarID]->getIsOnMultiscale()))
                {
                  return false;
                }
              else if(isMultiMultiReactantID[tarID])
                {
                  //If it meets the reaction probability:
                  if(theReactionProbabilities[tarID] == 1 ||
                     theRng.Fixed() < theReactionProbabilities[tarID])
                    { 
                      Voxel* target(&theLattice[coord]);
                      if(!theSpecies[tarID]->getIsMultiscale())
                        {
                          idx = theSpecies[tarID
                            ]->getTag(idx%theStride).multiIdx;
                        }
                      unsigned aMoleculeSize(theMoleculeSize);
                      react(target, target, srcIndex, idx%theStride,
                            idx/theStride);
                      if(theMoleculeSize < aMoleculeSize)
                        {
                          --srcIndex;
                        }
                      return false;
                    }
                  isIntersect = true;
                }
              //If it performs MultiMulti reaction and
              //intersects with another non-reacting Multi:
              else if(theSpecies[tarID]->getIsMultiscale() ||
                      theSpecies[tarID]->getIsOnMultiscale())
                {
                  isIntersect = true;
                }
              //If it doesn't intersect another Multi:
              else
                {
                  if(isMultiscaleBinderID[tarID])
                    {
                      ++tarCnt;
                    }
                  if(getID(theLattice[coord]) ==
                     theMultiscaleVacantSpecies->getID())
                    {
                      theIdxList.push_back(coord);
                    }
                  else
                    {
                      theBindList.push_back(coord);
                    }
                }
            }
        }
      if(isIntersect)
        {
          return false;
        }
      //count src
      for(unsigned i(0); i != srcOffsets.size(); ++i)
        {
          const int offsetRow((srcOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coordB(coordA+srcOffsets[i]);
          if(isInLattice(coordB, offsetRow))
            {
              const unsigned coord(coordB+lipStartCoord);
              if(isMultiscaleBoundID[getID(theLattice[coord])])
                {
                  ++srcCnt;
                }
              if(getID(theLattice[coord]) == theID)
                {
                  theVacantList.push_back(coord);
                }
              else
                {
                  theUnbindList.push_back(coord);
                }
            }
        }
      if(theRng.Fixed() < exp((tarCnt-srcCnt)/theWalkPropensity))
        {
          updateMoveMultiscale(index);
          return true;
        }
      return false;
    }
  void updateMoveMultiscale(const unsigned index)
    {
      const unsigned idx(index+theID*theStride);
      for(unsigned i(0); i != theIdxList.size(); ++i)
        {
          theLattice[theIdxList[i]].idx = idx;
        }
      for(unsigned i(0); i != theBindList.size(); ++i)
        {
          bindMultiscale(&theLattice[theBindList[i]], idx);
        }
      for(unsigned i(0); i != theVacantList.size(); ++i)
        {
          theLattice[theVacantList[i]].idx = 
            theMultiscaleVacantSpecies->getID()*theStride;
        }
      for(unsigned i(0); i != theUnbindList.size(); ++i)
        {
          unbindMultiscale(&theLattice[theUnbindList[i]]);
        }
    }
  bool isNonIntersectMultiscaleWalkPropensityRegular(const unsigned coordA, 
                                                     const unsigned rowA,
                                           const std::vector<int>& tarOffsets,
                                           const std::vector<int>& srcOffsets)
    {
      int tarCnt(0);
      int srcCnt(0);
      //count tar
      for(unsigned i(0); i != tarOffsets.size(); ++i)
        {
          const int offsetRow((tarOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coordB(coordA+tarOffsets[i]);
          if(isInLattice(coordB, offsetRow))
            {
              const unsigned anID(getID(theLattice[coordB+lipStartCoord]));
              if(isMultiscaleBinderID[anID])
                {
                  ++tarCnt;
                }
              else if(theSpecies[anID]->getIsMultiscale() || 
                      theSpecies[anID]->getIsOnMultiscale())
                {
                  return false;
                }
            }
        }
      //count src
      for(unsigned i(0); i != srcOffsets.size(); ++i)
        {
          const int offsetRow((srcOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coordB(coordA+srcOffsets[i]);
          if(isInLattice(coordB, offsetRow))
            {
              if(isMultiscaleBoundID[getID(theLattice[coordB+lipStartCoord])])
                {
                  ++srcCnt;
                }
            }
        }
      if(theRng.Fixed() < exp((tarCnt-srcCnt)/theWalkPropensity))
        {
          return true;
        }
      return false;
    }
  bool isMultiscaleWalkPropensityRegular(const unsigned coordA, 
                                         const unsigned rowA,
                                         const std::vector<int>& tarOffsets,
                                         const std::vector<int>& srcOffsets)
    {
      int tarCnt(0);
      int srcCnt(0);
      //count tar
      for(unsigned i(0); i != tarOffsets.size(); ++i)
        {
          const int offsetRow((tarOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coordB(coordA+tarOffsets[i]);
          if(isInLattice(coordB, offsetRow))
            {
              if(isMultiscaleBinderID[getID(theLattice[coordB+lipStartCoord])])
                {
                  ++tarCnt;
                }
            }
        }
      //count src
      for(unsigned i(0); i != srcOffsets.size(); ++i)
        {
          const int offsetRow((srcOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coordB(coordA+srcOffsets[i]);
          if(isInLattice(coordB, offsetRow))
            {
              if(isMultiscaleBoundID[getID(theLattice[coordB+lipStartCoord])])
                {
                  ++srcCnt;
                }
            }
        }
      if(theRng.Fixed() < exp((tarCnt-srcCnt)/theWalkPropensity))
        {
          return true;
        }
      return false;
    }
  bool isMultiscaleWalkPropensity(const unsigned srcCoord,
                                  const unsigned tarCoord)
    {
      unsigned srcCnt(0);
      unsigned tarCnt(0);
      const unsigned coordA(srcCoord-vacStartCoord);
      const unsigned coordB(tarCoord-vacStartCoord);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coord(theIntersectLipids[coordA][i]+lipStartCoord);
          if(isMultiscaleBoundID[getID(theLattice[coord])])
            {
              ++srcCnt;
            }
        }
      for(unsigned i(0); i != theIntersectLipids[coordB].size(); ++i)
        {
          unsigned coord(theIntersectLipids[coordB][i]+lipStartCoord);
          if(isMultiscaleBinderID[getID(theLattice[coord])])
            {
              ++tarCnt;
            }
        }
      if(theRng.Fixed() < exp((tarCnt-srcCnt)/theWalkPropensity))
        {
          return true;
        }
      return false;
    }
  void addMultiscaleMolecule(Voxel* aVoxel, const unsigned index)
    {
      const unsigned coordA(aVoxel->coord-vacStartCoord);
      const unsigned idx(index+theID*theStride);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2][0]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  const unsigned coord(coordB+lipStartCoord);
                  if(getID(theLattice[coord]) == 
                     theMultiscaleVacantSpecies->getID())
                    {
                      theLattice[coord].idx = idx;
                    }
                  else
                    {
                      bindMultiscale(&theLattice[coord], idx);
                    }
                }
            }
        }
      else
        {
          for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
            {
              unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
              if(getID(theLattice[coordB]) ==
                 theMultiscaleVacantSpecies->getID())
                {
                  theLattice[coordB].idx = idx;
                }
              else
                { 
                  bindMultiscale(&theLattice[coordB], idx);
                }
            }
        }
    }
  void moveMultiscaleMoleculeRegular(const unsigned coordA, 
                                     const unsigned rowA,
                                     const std::vector<int>& tarOffsets,
                                     const std::vector<int>& srcOffsets,
                                     const unsigned index)
    {
      const unsigned idx(index+theID*theStride);
      //Add tar
      for(unsigned i(0); i != tarOffsets.size(); ++i)
        {
          const int offsetRow((tarOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coordB(coordA+tarOffsets[i]);
          if(isInLattice(coordB, offsetRow))
            {
              const unsigned coord(coordB+lipStartCoord);
              if(getID(theLattice[coord]) ==
                 theMultiscaleVacantSpecies->getID())
                {
                  theLattice[coord].idx = idx;
                }
              else
                {
                  bindMultiscale(&theLattice[coord], idx);
                }
            }
        }
      //Remove src
      for(unsigned i(0); i != srcOffsets.size(); ++i)
        {
          const int offsetRow((srcOffsets[i]+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coordB(coordA+srcOffsets[i]);
          if(isInLattice(coordB, offsetRow))
            {
              const unsigned coord(coordB+lipStartCoord);
              if(getID(theLattice[coord]) == theID)
                {
                  theLattice[coord].idx = 
                    theMultiscaleVacantSpecies->getID()*theStride;
                }
              else
                {
                  unbindMultiscale(&theLattice[coord]);
                }
            }
        }
    }
  void updateBoundMultiIdx(Voxel* aVoxel, const unsigned rotIndex)
    {
      const unsigned coordA(aVoxel->coord-vacStartCoord);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2][rotIndex]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  const unsigned coord(coordB+lipStartCoord);
                  const unsigned anID(getID(theLattice[coord]));
                  if(anID == theID)
                    {
                      theLattice[coord].idx = aVoxel->idx;
                    }
                  else
                    {
                      theSpecies[anID]->getTag(&theLattice[coord]).multiIdx = 
                        aVoxel->idx;
                    }
                }
            }
        }
      else
        {
        }
    }
  void removeMultiscaleMolecule(Voxel* aVoxel, const unsigned rotIndex)
    {
      const unsigned coordA(aVoxel->coord-vacStartCoord);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2][rotIndex]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  const unsigned coord(coordB+lipStartCoord);
                  if(getID(theLattice[coord]) == theID)
                    {
                      theLattice[coord].idx = 
                        theMultiscaleVacantSpecies->getID()*theStride;
                    }
                  else
                    {
                      unbindMultiscale(&theLattice[coord]);
                    }
                }
            }
        }
      else
        {
          for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
            {
              unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
              if(getID(theLattice[coordB]) == theID)
                {
                  theLattice[coordB].idx =
                    theMultiscaleVacantSpecies->getID()*theStride;
                }
              else
                {
                  unbindMultiscale(&theLattice[coordB]);
                }
            }
        }
    }
  bool isMultiIntersectCoord(const unsigned srcCoord)
    {
      const unsigned coordA(srcCoord-vacStartCoord);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2][0]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  Species* aSpecies(theSpecies[
                                    getID(theLattice[coordB+lipStartCoord])]);
                  if(aSpecies->getIsMultiscale() ||
                     aSpecies->getIsOnMultiscale())
                    {
                      return true;
                    }
                }
            }
        }
      else
        {
          for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
            {
              unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
              unsigned anID(getID(theLattice[coordB]));
              if(theSpecies[anID]->getIsMultiscale() ||
                 theSpecies[anID]->getIsOnMultiscale())
                {
                  return true;
                }
            }
        }
      return false;
    }
  bool isMultiIntersectCoord(const unsigned srcCoord,
                             std::vector<unsigned>& multiCompCnts)
    {
      std::vector<unsigned> cnts;
      cnts.resize(theSpecies.size(), 0);
      const unsigned coordA(srcCoord-vacStartCoord);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2][0]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  Species* aSpecies(theSpecies[
                                    getID(theLattice[coordB+lipStartCoord])]);
                  if(aSpecies->getIsMultiscale() ||
                     aSpecies->getIsOnMultiscale())
                    {
                      return true;
                    }
                  else
                    {
                      ++cnts[aSpecies->getID()];
                    }
                }
            }
        }
      else
        {
          //do for non regular lattice
        }
      multiCompCnts = cnts;
      return false;
    }
  bool isMultiIntersectCoord(const unsigned coordA, const unsigned ignoreIdx)
    {
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2][0]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  const unsigned idx(theLattice[coordB+lipStartCoord].idx);
                  Species* aSpecies(theSpecies[idx/theStride]);
                  if((aSpecies->getIsMultiscale() && idx != ignoreIdx) ||
                     (aSpecies->getIsOnMultiscale() && 
                      aSpecies->getTag(idx%theStride).multiIdx != ignoreIdx))
                    {
                      return true;
                    }
                }
            }
        }
      else
        {
          //do for non regular lattice
        }
      return false;
    }
  bool isMultiIntersectCoord(const unsigned coordA, const unsigned ignoreIdx,
                             const unsigned ignoreIdx2)
    {
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2][0]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coordB(coordA+anOffsets[i]);
              if(isInLattice(coordB, offsetRow+rowA))
                {
                  const unsigned idx(theLattice[coordB+lipStartCoord].idx);
                  Species* aSpecies(theSpecies[idx/theStride]);
                  if((aSpecies->getIsMultiscale() && idx != ignoreIdx &&
                      idx != ignoreIdx2) ||
                     (aSpecies->getIsOnMultiscale() && 
                      aSpecies->getTag(idx%theStride).multiIdx != ignoreIdx &&
                      aSpecies->getTag(idx%theStride).multiIdx != ignoreIdx2))
                    {
                      return true;
                    }
                }
            }
        }
      else
        {
          //do for non regular lattice
        }
      return false;
    }
  unsigned getLipStartCoord() const
    {
      return lipStartCoord;
    }
  unsigned getVacStartCoord() const
    {
      return vacStartCoord;
    }
  bool isInLatticeOrigins(int& coord, const int tarRow, Origin& anOrigin) const
    {
      if(isPeriodic)
        {
          if(coord < 0)
            {
              if((coord+1)/lipCols+1 > -tarRow)
                {
                  --anOrigin.col;
                }
              else if((coord+1)/lipCols+1 < -tarRow)
                {
                  ++anOrigin.col;
                }
              coord += lipCols*(tarRow-(coord+1)/lipCols+1);
            }
          else
            {
              if(coord/lipCols < tarRow)
                {
                  --anOrigin.col;
                }
              else if(coord/lipCols > tarRow)
                {
                  ++anOrigin.col;
                }
              coord += lipCols*(tarRow-(coord/lipCols));
            }
          if(coord < 0)
            {
              coord += lipRows*lipCols;
              --anOrigin.row;
            }
          else if(coord >= lipRows*lipCols)
            {
              coord -= lipRows*lipCols;
              ++anOrigin.row;
            }
        }
      else if(coord/lipCols != tarRow || coord < 0 ||
              coord >= lipRows*lipCols)
        {
          return false;
        }
      return true;
    }
  bool isInLattice(int& coord, const int rowOffset) const
    {
      if(isPeriodic)
        {
          if(coord < 0)
            {
              coord += lipCols*(rowOffset-(coord+1)/lipCols+1);
            }
          else
            {
              coord += lipCols*(rowOffset-(coord/lipCols));
            }
          if(coord < 0)
            {
              coord += lipRows*lipCols;
            }
          else if(coord >= lipRows*lipCols)
            {
              coord -= lipRows*lipCols;
            }
        }
      else if(coord/lipCols != rowOffset || coord < 0 ||
              coord >= lipRows*lipCols)
        {
          return false;
        }
      return true;
    }
  bool isIntersectMultiscale(const unsigned srcCoord, const unsigned tarCoord)
    {
      bool isIntersect(false);
      /*
      const unsigned coordA(srcCoord-vacStartCoord);
      std::vector<unsigned> temp;
      temp.resize(theIntersectLipids[coordA].size());
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
          temp[i] = theLattice[coordB].id;
          theLattice[coordB].id = theSpeciesSize;
        }
      isIntersect = isIntersectMultiscale(tarCoord);
      for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
        {
          unsigned coordB(theIntersectLipids[coordA][i]+lipStartCoord);
          theLattice[coordB].id = temp[i];
        }
        */
      return isIntersect;
    }
  void bindMultiscale(Voxel* aVoxel, const unsigned multiIdx)
    { 
      unsigned anID(getID(aVoxel));
      theDiffusionInfluencedReactions[anID]->bind(aVoxel, multiIdx);
      isFinalizeReactions[anID] = true;
    }
  void unbindMultiscale(Voxel* aVoxel)
    {
      unsigned anID(theMultiscaleUnbindIDs[getID(aVoxel)]);
      theDiffusionInfluencedReactions[anID]->unbind(aVoxel);
      isFinalizeReactions[anID] = true;
    }
  void addCompVoxel(const unsigned aCoord)
    {
      theLattice[aCoord].idx = theMoleculeSize+theStride*theID;
      ++theMoleculeSize; 
      theCompVoxels->push_back(&theLattice[aCoord]);
      theVariable->setValue(theMoleculeSize);
    }
  //must use anIndex to get the voxel. Using coord and getIndex may not work.
  void removeCompVoxel(const unsigned anIndex)
    {
      --theMoleculeSize;
      (*theCompVoxels)[anIndex]->idx = theVacantID*theStride;
      (*theCompVoxels)[anIndex] = (*theCompVoxels)[theMoleculeSize];
      theCompVoxels->pop_back();
      theVariable->setValue(theMoleculeSize);
    }
  void clearCompVoxels()
    {
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          (*theCompVoxels)[i]->idx = theVacantID*theStride;
        }
      theMoleculeSize = 0;
      theCompVoxels->resize(0);
      theTags.resize(0);
      theVariable->setValue(0);
    }
  void clearEmptyCompVoxels()
    {
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          if(getID((*theCompVoxels)[i]) == getID())
            {
              (*theCompVoxels)[i]->idx = theVacantID*theStride;
            }
        }
      theMoleculeSize = 0;
      theCompVoxels->resize(0);
      theTags.resize(0);
      theVariable->setValue(0);
    }
  String getIDString(unsigned anID) const
    {
      return theSpecies[anID]->getIDString();
    }
  String getIDString(Species* aSpecies) const
    {
      return aSpecies->getIDString();
    }
  String getIDString() const
    {
      const Variable* aVariable(getVariable());
      if(aVariable)
        {
          return "["+aVariable->getSystemPath().asString()+":"+
            aVariable->getID()+"]["+int2str(getID())+"]";
        }
      else if(getID() == theNullID)
        {
          return "[theNullID]["+int2str(getID())+"]";
        }
      return "[unknown]";
    }
  unsigned compVoxelSize()
    {
      return theCompVoxels->size();
    }
  Voxel* getCompVoxel(unsigned index)
    {
      return (*theCompVoxels)[index];
    }
  const unsigned getIndex(const Voxel* aVoxel) const
    {
      return aVoxel->idx%theStride;
    }
  //It is soft remove because the id of the molecule is not changed:
  void softRemoveMolecule(Voxel* aVoxel)
    {
      softRemoveMolecule(getIndex(aVoxel));
    }
  void removeMolecule(Voxel* aVoxel)
    {
      removeMolecule(getIndex(aVoxel));
    }
  void removeMolecule(unsigned anIndex)
    {
      if(!isVacant)
        {
          theMolecules[anIndex]->idx = theVacantID*theStride;
          softRemoveMolecule(anIndex);
        }
    }
  void softRemoveMolecule(unsigned anIndex)
    {
      if(isDeoligomerize)
        {
          removeBounds(anIndex);
          return;
        }
      if(!isVacant)
        {
          if(isMultiscale)
            {
              removeMultiscaleMolecule(theMolecules[anIndex],
                                       theTags[anIndex].rotIndex);
            }
          removeMoleculeDirect(anIndex);
        }
    }
  void addBound(const unsigned index)
    {
      theBoundCnts[theTags[index].boundCnt++]--;
      theBoundCnts[theTags[index].boundCnt]++;
    }
  void removeBound(const unsigned index)
    {
      Tag& aTag(theTags[index]);
      theBoundCnts[aTag.boundCnt--]--;
      if(!aTag.boundCnt)
        {
          Voxel* aVoxel(theMolecules[index]);
          if(isOnMultiscale)
            {
              theDeoligomerizedProduct->addMoleculeInMulti(aVoxel,
                                                           aTag.multiIdx);
            }
          else
            {
              theDeoligomerizedProduct->addMolecule(aVoxel, aTag);
            }
          removeMoleculeDirect(index);
        }
      else
        {
          theBoundCnts[aTag.boundCnt]++;
        }
    }
  void addBounds(Voxel* aVoxel)
    {
      unsigned& cnt(theTags[getIndex(aVoxel)].boundCnt);
      cnt = 0;
      for(unsigned i(0); i != aVoxel->diffuseSize; ++i)
        {
          unsigned aCoord(aVoxel->adjoiningCoords[i]);
          const unsigned anID(getID(theLattice[aCoord]));
          if(theSpecies[anID]->getIsDeoligomerize())
            {
              cnt++;
              theSpecies[anID]->addBound(theLattice[aCoord].idx%theStride);
            }
        }
      theBoundCnts[cnt]++;
    }
  void removeBounds(const unsigned anIndex)
    {
      Voxel* aVoxel(theMolecules[anIndex]);
      theBoundCnts[theTags[anIndex].boundCnt]--;
      removeMoleculeDirect(anIndex);
      removeAdjoinBounds(aVoxel);
    }
  void removeAdjoinBounds(Voxel* aVoxel)
    {
      for(unsigned i(0); i != aVoxel->diffuseSize; ++i)
        {
          const unsigned aCoord(aVoxel->adjoiningCoords[i]);
          const unsigned anID(getID(theLattice[aCoord]));
          if(theSpecies[anID]->getIsDeoligomerize())
            {
              theSpecies[anID]->removeBound(theLattice[aCoord].idx%theStride);
            }
        }
    }
  unsigned getBoundCnt(const unsigned anIndex)
    {
      return theBoundCnts[anIndex];
    }
  std::vector<unsigned>& getBoundCnts()
    {
      return theBoundCnts;
    }
  void removeMoleculeBoundDirect(unsigned anIndex)
    {
      if(isDeoligomerize)
        {
          theBoundCnts[theTags[anIndex].boundCnt]--;
        }
      --theMoleculeSize;
      if(theMoleculeSize > anIndex)
        {
          theMolecules[anIndex] = theMolecules[theMoleculeSize];
          theMolecules[anIndex]->idx = anIndex+theStride*theID;
          theTags[anIndex] = theTags[theMoleculeSize];
        }
      theVariable->setValue(theMoleculeSize);
    }
  void removeMoleculeDirect(unsigned anIndex)
    {
      for(unsigned i(0); i != theInterruptedProcessesRemoveMolecule.size(); ++i)
        {
          theInterruptedProcessesRemoveMolecule[i
            ]->interruptedRemoveMolecule(this, anIndex);
        }
      --theMoleculeSize;
      if(theMoleculeSize > anIndex)
        {
          theMolecules[anIndex] = theMolecules[theMoleculeSize];
          theMolecules[anIndex]->idx = anIndex+theStride*theID;
          theTags[anIndex] = theTags[theMoleculeSize];
          if(isMultiscale)
            {
              updateBoundMultiIdx(theMolecules[anIndex],
                                  theTags[anIndex].rotIndex);
            }
        }
      theVariable->setValue(theMoleculeSize);
    }
  //Used to remove all molecules and free memory used to store the molecules
  void clearMolecules()
    {
      theMolecules.resize(0);
      theTags.resize(0);
      theMoleculeSize = 0;
      theVariable->setValue(0);
    }
  //Used by the SpatiocyteStepper when resetting an interation, so must
  //clear the whole compartment using theComp->vacantSpecies->getVacantID():
  void removeMolecules()
    {
      if(!isCompVacant && !isInterface)
        {
          while(theMoleculeSize)
            {
              removeMolecule(theMoleculeSize-1);
            }
          theVariable->setValue(theMoleculeSize);
          if(theDeoligomerizedProduct)
            {
              theDeoligomerizedProduct->removeMolecules();
            }
        }
    }
  int getPopulateCoordSize()
    {
      return theInitCoordSize-theMoleculeSize;
    }
  int getInitCoordSize()
    {
      return theInitCoordSize;
    }
  void resetMoleculeOrigins()
    {
      //For GFP species, whose molecules are not always up to date:
      if(isTag)
        {
          for(unsigned i(0); i != theTaggedSpeciesList.size(); ++i)
            {
              theTaggedSpeciesList[i]->resetMoleculeOrigins();
            }
          updateMolecules();
          return;
        }
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          resetTagOrigin(i);
        }
    }
  void resetTagOrigin(const unsigned index)
    {
      Origin& anOrigin(theTags[index].origin);
      anOrigin.point = getPoint(index);
      anOrigin.row = 0;
      anOrigin.layer = 0;
      anOrigin.col = 0;
    }
  void setInitCoordSize(const unsigned& val)
    {
      theInitCoordSize = val;
    }
  void initMoleculeOrigins()
    {
      isOrigins = true;
      resetMoleculeOrigins();
    }
  void removeBoundaryMolecules()
    {
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          if(theStepper->isBoundaryCoord(getCoord(i), theDimension))
            {
              cout << "is still there" << std::endl;
            }
        }
      theVariable->setValue(theMoleculeSize);
    }
  void relocateBoundaryMolecules()
    {
      for(unsigned i(0); i < theMoleculeSize; ++i)
        {
          Origin anOrigin(theTags[i].origin);
          unsigned periodicCoord(theStepper->getPeriodicCoord(
                                                getCoord(i),
                                                theDimension, &anOrigin));
          if(getID(theLattice[periodicCoord]) == theVacantID)
            {
              theMolecules[i]->idx = theVacantID*theStride;
              theMolecules[i] = &theLattice[periodicCoord];
              theMolecules[i]->idx = i+theID*theStride;
              theTags[i].origin = anOrigin;
            }
        }
    }
  unsigned getVacantID() const
    {
      return theVacantID;
    }
  Species* getVacantSpecies() const
    {
      return theVacantSpecies;
    }
  //Called by initializeFirst of CompartmentProcess:
  void setVacantSpecies(Species* aVacantSpecies)
    {
      theVacantSpecies = aVacantSpecies;
      theVacantID = aVacantSpecies->getID();
      theVacantIdx = theStride*theVacantID;
      if(aVacantSpecies->getIsMultiscale())
        {
          isOnMultiscale = true;
          isMultiscaleComp = true;
        }
    }
  void setTrailSpecies(Species* aTrailSpecies)
    {
      theTrailSpecies = aTrailSpecies;
    }
   const std::vector<double>& getBendAngles() const
    {
      return theBendAngles;
    }
  const Point getWestPoint() const
    {
      Point aWestPoint(theComp->centerPoint);
      aWestPoint.x = theComp->centerPoint.x-theComp->lengthX/2+
        theComp->lengthY/2;
      return aWestPoint;
    }
  const Point getEastPoint() const
    {
      Point anEastPoint(theComp->centerPoint); 
      anEastPoint.x = theComp->centerPoint.x+theComp->lengthX/2-
        theComp->lengthY/2;
      return anEastPoint;
    }
  double getCompRadius() const
    {
      return theComp->lengthY/2;
    }
  double getMoleculeRadius() const
    {
      return theMoleculeRadius;
    }
  double getDiffuseRadius() const
    {
      return theDiffuseRadius;
    }
  void setMoleculeRadius(double aRadius)
    {
      theMoleculeRadius = aRadius;
      theDiffuseRadius = aRadius;
      nDiffuseRadius = theDiffuseRadius/(theStepper->getVoxelRadius()*2);
    }
  void setDiffuseRadius(double aRadius)
    {
      theDiffuseRadius = aRadius;
      nDiffuseRadius = theDiffuseRadius/(theStepper->getVoxelRadius()*2);
    }
  Species* getDiffusionInfluencedReactantPair()
    {
      if(theDiffusionInfluencedReactantPairs.empty())
        {
          return NULL;
        }
      return theDiffusionInfluencedReactantPairs[0];
    }
  double getReactionProbability(int anID)
    {
      return theReactionProbabilities[anID];
    }
  double getMaxReactionProbability()
    {
      double maxProbability(0);
      for(std::vector<double>::const_iterator
          i(theReactionProbabilities.begin()); 
          i != theReactionProbabilities.end(); ++i)
        {
          if(maxProbability < *i)
            {
              maxProbability = *i;
            }
        }
      return maxProbability;
    }
  void rescaleReactionProbabilities(double aWalkProbability)
    {
      theWalkProbability = aWalkProbability;
      for(std::vector<double>::iterator i(theReactionProbabilities.begin()); 
          i != theReactionProbabilities.end(); ++i)
        {
          *i = (*i)*aWalkProbability;
        }
    }
  void setDiffusionInterval(double anInterval)
    {
      if(anInterval < theDiffusionInterval)
        {
          theDiffusionInterval = anInterval;
        }
      for(unsigned i(0); i != theTagSpeciesList.size(); ++i)
        {
          theTagSpeciesList[i]->setDiffusionInterval(theDiffusionInterval);
        }
    }
  void addInterruptAddMolecule(SpatiocyteProcess* aProcess)
    {
      if(std::find(theInterruptedProcessesAddMolecule.begin(),
                   theInterruptedProcessesAddMolecule.end(), aProcess) ==
         theInterruptedProcessesAddMolecule.end())
        {
          theInterruptedProcessesAddMolecule.push_back(aProcess);
        }
    }
  void addInterruptRemoveMolecule(SpatiocyteProcess* aProcess)
    {
      if(std::find(theInterruptedProcessesRemoveMolecule.begin(),
                   theInterruptedProcessesRemoveMolecule.end(), aProcess) ==
         theInterruptedProcessesRemoveMolecule.end())
        {
          theInterruptedProcessesRemoveMolecule.push_back(aProcess);
        }
    }
  void addInterruptEndDiffusion(SpatiocyteProcess* aProcess)
    {
      if(std::find(theInterruptedProcessesEndDiffusion.begin(),
                   theInterruptedProcessesEndDiffusion.end(), aProcess) ==
         theInterruptedProcessesEndDiffusion.end())
        {
          theInterruptedProcessesEndDiffusion.push_back(aProcess);
        }
    }
  unsigned getRandomIndex()
    {
      return theRng.Integer(theMoleculeSize);
    }
  unsigned getRandomValidIndex()
    {
      unsigned anIndex(getRandomIndex());
      const unsigned ranIndex(anIndex);
      while(getID(theMolecules[anIndex]) != theID)
        {
          ++anIndex;
          if(anIndex == theMoleculeSize)
            {
              anIndex = 0;
            }
          else if(anIndex == ranIndex)
            {
              return theMoleculeSize;
            }
        }
      return anIndex;
    }
  unsigned getRandomOligomerIndex(const unsigned boundCnt)
    {
      const unsigned start(theRng.Integer(theMoleculeSize));
      for(unsigned i(start); i != theMoleculeSize; ++i)
        {
          if(theTags[i].boundCnt == boundCnt)
            {
              return i;
            }
        }
      for(unsigned i(0); i != start; ++i)
        {
          if(theTags[i].boundCnt == boundCnt)
            {
              return i;
            }
        }
      cout << "shouldn't get here deoligomerize:" << getIDString() <<
        std::endl;
      return 0;
    }
  Voxel* getRandomMolecule()
    {
      return theMolecules[getRandomIndex()];
    }
  int getBendIndex(double aBendAngle)
    {
      for(unsigned i(0); i != theBendAngles.size(); ++i)
        {
          if(theBendAngles[i] == aBendAngle)
            {
              return i;
            }
        }
      return 0;
    }
  Voxel* getRandomAdjoiningVoxel(Voxel* source, int searchVacant)
    {
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(isPopulatable(&theLattice[aCoord]))
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theStepper->id2Comp(getID(theLattice[aCoord])) == theComp)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantVoxel(compCoords);
    } 
  Voxel* getBindingSiteAdjoiningVoxel(Voxel* source, int bindingSite)
    {
      if(bindingSite < source->diffuseSize)
        { 
          unsigned aCoord(source->adjoiningCoords[bindingSite]);
          if(isPopulatable(&theLattice[aCoord]))
            {
              return &theLattice[aCoord];
            }
        }
      return NULL;
    } 
  Voxel* getBindingSiteAdjoiningVoxel(Voxel* source, int bindingSite,
                                      Species* aSpecies)
    {
      if(bindingSite < source->diffuseSize)
        { 
          unsigned aCoord(source->adjoiningCoords[bindingSite]);
          if(getID(&theLattice[aCoord]) == aSpecies->getID())
            {
              return &theLattice[aCoord];
            }
        }
      return NULL;
    } 
  Voxel* getRandomDiffuseVoxel(Voxel* source, Species* aTargetSpecies,
                               int searchVacant)
    {
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source->diffuseSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(getID(theLattice[aCoord]) == aTargetSpecies->getID())
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source->diffuseSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theStepper->id2Comp(getID(theLattice[aCoord])) == theComp)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantVoxel(compCoords, aTargetSpecies);
    } 
  Voxel* getRandomAdjoiningVoxel(Voxel* source,
                                 Species* aTargetSpecies,
                                 int searchVacant)
    {
      std::vector<unsigned> compCoords;
      if(searchVacant)
        {
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(getID(theLattice[aCoord]) == aTargetSpecies->getID())
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theStepper->id2Comp(getID(theLattice[aCoord])) == theComp)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantVoxel(compCoords, aTargetSpecies);
    } 
  Voxel* getMultiRandomAdjoiningVoxel(Species* A, Species* D, Voxel* molA,
                                      const unsigned searchVacant)
    {
      //molA will be replaced with molecule from this species (C) after the
      //reaction:
      const unsigned coordA(molA->coord-vacStartCoord);
      const unsigned idxA(molA->idx);
      const unsigned row(coordA/lipCols);
      if(!isMultiIntersectCoord(coordA, idxA))
        {
          std::vector<int>& anOffsets(theProductPairOffsets[D->getID()][row%2]);
          std::vector<unsigned> doneAdjs;
          unsigned i;
          do
            {
              do
                {
                  i = theRng.Integer(anOffsets.size()); 
                }
              while(searchVacant && std::find(doneAdjs.begin(), doneAdjs.end(),
                                              i) != doneAdjs.end());
              doneAdjs.push_back(i);
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols+row);
              int coordD(coordA+anOffsets[i]);
              if(isInLattice(coordD, offsetRow) &&
                 !D->isMultiIntersectCoord(coordD, idxA))
                {
                  return &theLattice[coordD+vacStartCoord];
                }
            }
          while(searchVacant && doneAdjs.size() < anOffsets.size());
        }
      return NULL;
    }
  unsigned getMultiCoordSize() const
    {
      return theOffsets[0][0].size();
    }
  void getMultiCoords(const unsigned molIndex,
                      std::vector<unsigned>& coords) const
    {
      const unsigned coordA(getCoord(molIndex)-vacStartCoord);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const std::vector<int>& anOffsets(theOffsets[rowA%2][0]);
          for(unsigned i(0); i != anOffsets.size(); ++i)
            {
              const int offsetRow((anOffsets[i]+theRegLatticeCoord)/lipCols-
                                  theRegLatticeCoord/lipCols);
              int coord(coordA+anOffsets[i]);
              if(isInLattice(coord, offsetRow+rowA))
                {
                  coords.push_back(coord+lipStartCoord);
                }
            }
        }
    }
  unsigned getMultiCoord(const unsigned molIndex,
                         const unsigned offsetIndex) const
    { 
      const unsigned coordA(getCoord(molIndex)-vacStartCoord);
      if(isRegularLattice)
        {
          const int rowA(coordA/lipCols);
          const int anOffset(theOffsets[rowA%2][0][offsetIndex]);
          const int offsetRow((anOffset+theRegLatticeCoord)/lipCols-
                              theRegLatticeCoord/lipCols+rowA);
          int coord(coordA+anOffset);
          if(isInLattice(coord, offsetRow))
            {
              return coord;
            }
        }
      return theNullCoord;
    }
  Voxel* getRandomAdjoiningVoxel(Species* A, Species* C, Voxel* molA,
                                 const unsigned searchVacant)
    {
      if(getIsMultiscale())
        {
          return C->getMultiRandomAdjoiningVoxel(A, this, molA, searchVacant);
        }
      return getRandomAdjoiningVoxel(molA, searchVacant);
    }
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Voxel* target,
                                 int searchVacant)
    {
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(&theLattice[aCoord] != target && 
                 isPopulatable(&theLattice[aCoord]))
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theStepper->id2Comp(getID(theLattice[aCoord])) == theComp &&
                 &theLattice[aCoord] != target)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantVoxel(compCoords);
    }
  bool isAdjoinedSpecies(Voxel* source, Species* aTargetSpecies)
    {
      for(unsigned i(0); i != source->adjoiningSize; ++i)
        {
          if(getID(theLattice[source->adjoiningCoords[i]]) == 
             aTargetSpecies->getID())
            {
              return true;
            }
        }
      return false;
    }
  unsigned getAdjoiningMoleculeCnt(Voxel* source, Species* aTargetSpecies)
    {
      unsigned cnt(0);
      for(unsigned i(0); i != source->adjoiningSize; ++i)
        {
          if(getID(theLattice[source->adjoiningCoords[i]]) == 
             aTargetSpecies->getID())
            {
              ++cnt;
            }
        }
      return cnt;
    }
  Voxel* getRandomAdjoiningVoxel(Voxel* source, Voxel* targetA, Voxel* targetB,
                                 int searchVacant)
    {
      std::vector<unsigned> compCoords;
      if(searchVacant)
        { 
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(source != targetA && source != targetB && 
                 isPopulatable(&theLattice[aCoord]))
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      else
        {
          for(unsigned i(0); i != source->adjoiningSize; ++i)
            {
              unsigned aCoord(source->adjoiningCoords[i]);
              if(theStepper->id2Comp(getID(theLattice[aCoord])) == theComp &&
                 source != targetA && source != targetB)
                {
                  compCoords.push_back(aCoord);
                }
            }
        }
      return getRandomVacantVoxel(compCoords);
    }
  Voxel* getRandomVacantVoxel(std::vector<unsigned>& aCoords)
    {
      if(aCoords.size())
        {
          const unsigned r(theRng.Integer(aCoords.size())); 
          const unsigned aCoord(aCoords[r]);
          if(isPopulatable(&theLattice[aCoord]))
            {
              return &theLattice[aCoord];
            }
        }
      return NULL;
    }
  Voxel* getRandomVacantVoxel(std::vector<unsigned>& aCoords,
                              Species* aVacantSpecies)
    {
      if(aCoords.size())
        {
          const unsigned r(theRng.Integer(aCoords.size())); 
          const unsigned aCoord(aCoords[r]);
          if(getID(theLattice[aCoord]) == aVacantSpecies->getID())
            {
              return &theLattice[aCoord];
            }
        }
      return NULL;
    }
  Voxel* getRandomPopulatableMulti(int searchVacant,
                                   std::vector<unsigned>& multiCompCnts)
    {
      const unsigned aSize(theVacantSpecies->compVoxelSize());
      const unsigned r(theRng.Integer(aSize)); 
      if(searchVacant)
        {
          for(unsigned i(r); i != aSize; ++i)
            {
              Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
              if(!isMultiIntersectCoord(aVoxel->coord, multiCompCnts))
                {
                  return aVoxel;
                }
            }
          for(unsigned i(0); i != r; ++i)
            {
              Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
              if(!isMultiIntersectCoord(aVoxel->coord, multiCompCnts))
                {
                  return aVoxel;
                }
            }
        }
      else
        {
          Voxel* aVoxel(theVacantSpecies->getCompVoxel(r));
          if(!isMultiIntersectCoord(aVoxel->coord, multiCompCnts))
            {
              return aVoxel;
            }
        }
      return NULL;
    }
  Voxel* getRandomPopulatableVoxel(int searchVacant)
    {
      const unsigned aSize(theVacantSpecies->compVoxelSize());
      const unsigned r(theRng.Integer(aSize)); 
      if(searchVacant)
        {
          for(unsigned i(r); i != aSize; ++i)
            {
              Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
              if(isPopulatable(aVoxel))
                {
                  return aVoxel;
                }
            }
          for(unsigned i(0); i != r; ++i)
            {
              Voxel* aVoxel(theVacantSpecies->getCompVoxel(i));
              if(isPopulatable(aVoxel))
                {
                  return aVoxel;
                }
            }
        }
      else
        {
          Voxel* aVoxel(theVacantSpecies->getCompVoxel(r));
          if(isPopulatable(aVoxel))
            {
              return aVoxel;
            }
        }
      return NULL;
    }
  Voxel* getRandomAdjoiningCompVoxel(Comp* aComp, int searchVacant)
    {
      int aSize(theVacantSpecies->size());
      const unsigned r(theRng.Integer(aSize)); 
      Voxel* aVoxel(theVacantSpecies->getMolecule(r));
      return getRandomAdjoiningVoxel(aVoxel, searchVacant);
    }
  //We need to updateMolecules to set the valid address of voxels
  //since they may have been changed when theLattice is resized by 
  //processes:
  void updateMoleculePointers()
    {
      for(unsigned i(0); i != theCoords.size(); ++i)
        {
          theMolecules[i] = &theLattice[theCoords[i]];
        }
      theCoords.resize(0);
    }
  void saveCoords()
    {
      theCoords.resize(theMoleculeSize);
      for(unsigned i(0); i != theMoleculeSize; ++i)
        {
          theCoords[i] = getCoord(i);
        }
    }
  void setVacStartCoord(unsigned aCoord, unsigned aVacantRows,
                        unsigned aVacantCols)
    {
      vacStartCoord = aCoord;
      vacRows = aVacantRows;
      vacCols = aVacantCols;
    }
  void setLipStartCoord(unsigned aCoord, unsigned aLipidRows,
                        unsigned aLipidCols)
    {
      lipStartCoord = aCoord;
      lipRows = aLipidRows;
      lipCols = aLipidCols;
    }
  void setIntersectLipids(Species* aLipid, Point& aLipidStart, double aGridSize,
                          unsigned aGridCols, unsigned aGridRows, 
                          std::vector<std::vector<unsigned> >& aGrid,
                          unsigned aVacantRows, unsigned aVacantCols)
    {
      double nDist((aLipid->getMoleculeRadius()+theMoleculeRadius)/
                   (2*theVoxelRadius));
      if(aLipid->getMoleculeRadius() >= theMoleculeRadius)
        {
          nDist = 1.5*theMoleculeRadius/(2*theVoxelRadius);
        }
      //Traverse through the entire compartment voxels:
      unsigned endA(vacStartCoord+theVacantSpecies->size());
      theIntersectLipids.resize(theVacantSpecies->size());
      for(unsigned i(vacStartCoord); i != endA; ++i)
        {
          getIntersectLipids(i, nDist, aLipidStart, aGridSize,
                             aGridCols, aGridRows, aGrid,
                             theIntersectLipids[i-vacStartCoord]);
        }
    }
  void getIntersectLipids(unsigned coordA, double nDist, Point& aLipidStart,
                          double aGridSize, unsigned aGridCols, 
                          unsigned aGridRows, 
                          std::vector<std::vector<unsigned> >& aGrid,
                          std::vector<unsigned>& anIntersectLipids)
    {
      Point& pointA(*theLattice[coordA].point);
      unsigned row((unsigned)((pointA.y-aLipidStart.y)/aGridSize));
      unsigned col((unsigned)((pointA.z-aLipidStart.z)/aGridSize));
      unsigned rowStart(std::max(unsigned(1), row)-1);
      unsigned rowEnd(std::min(unsigned(aGridRows), row+2));
      for(unsigned j(rowStart); j != rowEnd; ++j)
        {
          unsigned colStart(std::max(unsigned(1), col)-1);
          unsigned colEnd(std::min(unsigned(aGridCols), col+2));
          for(unsigned k(colStart); k != colEnd; ++k)
            {
              std::vector<unsigned>& coords(aGrid[k+aGridCols*j]);
              for(unsigned l(0); l != coords.size(); ++l)
                {
                  unsigned m(coords[l]);
                  Point& pointB(*theLattice[m].point);
                  if(distance(pointA, pointB) < nDist)
                    {
                      anIntersectLipids.push_back(m-lipStartCoord);
                    }
                }
            }
        }
    }
  void setAdjoinOffsets(const std::vector<std::vector<int> > anAdjoinOffsets,
                        const std::vector<int> aRowOffsets)
    {
      theAdjoinOffsets = anAdjoinOffsets;
      theRowOffsets = aRowOffsets;
    }
  void setProductPairOffsets(const std::vector<std::vector<int> > 
                             anAdjoinOffsets,
                             const std::vector<int> aRowOffsets,
                             Species* aLipid, const Point& aLipidStart, 
                             const unsigned aVacantRows,
                             const unsigned aVacantCols,
                             const double nLipidRadius,
                             const double aSubunitAngle, Point& aSurfaceNormal,
                             Point& widthVector, Point& lengthVector)
    {
      unsigned rowA(aVacantRows/2);
      unsigned rowB(rowA+1);
      if(rowB >= aVacantRows && rowA > 0)
        {
          rowB = rowA-1;
        }
      else if(!rowA)
        {
          rowB = 0;
        }
      unsigned aCol(aVacantCols/2);
      unsigned coordA(rowA*aVacantCols+aCol);
      unsigned coordB(rowB*aVacantCols+aCol);
      theProductPairOffsets.resize(theSpecies.size());
      std::vector<unsigned> coordsA;
      for(unsigned i(0); i != theOffsets[rowA%2][0].size(); ++i)
        {
          coordsA.push_back(theOffsets[rowA%2][0][i]+coordA);
        }
      std::vector<unsigned> coordsB;
      for(unsigned i(0); i != theOffsets[rowB%2][0].size(); ++i)
        {
          coordsB.push_back(theOffsets[rowB%2][0][i]+coordB);
        }
      for(unsigned i(0); i != theSpecies.size(); ++i)
        {
          if(theProductPairOffsets[i].empty() && isProductPair[i])
            {
              double nDist((theMoleculeRadius+
                            2*theSpecies[i]->getMoleculeRadius())/
                           (2*theVoxelRadius));
              theProductPairOffsets[i].resize(2);
              pushProductPairOffsets(rowA, coordA, coordsA, nDist, aLipidStart,
                                     nLipidRadius, 
                                     theProductPairOffsets[i][rowA%2], i,
                                     widthVector, lengthVector);
              pushProductPairOffsets(rowB, coordB, coordsB, nDist, aLipidStart,
                                     nLipidRadius, 
                                     theProductPairOffsets[i][rowB%2], i,
                                     widthVector, lengthVector);
            }
        }
    }
  void getRowColStartEnd(const unsigned coordA, const double nDist,
                         const Point& aLipidStart, const double nLipidRadius,
                         unsigned& rowStart, unsigned& colStart,
                         unsigned& rowEnd, unsigned& colEnd,
                         Point& widthVector, Point& lengthVector)
    {
      Point& pointA(*theLattice[coordA+vacStartCoord].point);
      double minWidth(point2lineDist(pointA, lengthVector, aLipidStart));
      double minLength(point2lineDist(pointA, widthVector, aLipidStart));
      double maxWidth(minWidth+nDist*2);
      double maxLength(minLength+nDist*2);
      minWidth = std::max(0.0, minWidth-nDist*2);
      minLength = std::max(0.0, minLength-nDist*2);
      rowStart = (unsigned)std::max(minWidth/(nLipidRadius*sqrt(3))-1, 0.0);
      rowEnd = (unsigned)std::min(maxWidth/(nLipidRadius*sqrt(3))+1,
                                         double(lipRows));
      colStart = (unsigned)std::max(minLength/(nLipidRadius*2)-1, 0.0);
      colEnd = (unsigned)std::min(maxLength/(nLipidRadius*2)+1,
                                  double(lipCols));
    }
  void pushProductPairOffsets(const unsigned row, const unsigned coordA,
                              const std::vector<unsigned>& myCoords,
                              const double nDist, const Point& aLipidStart,
                              const double nLipidRadius,
                              std::vector<int>& aProductPairOffsets,
                              const unsigned aSpeciesID,
                              Point& widthVector, Point& lengthVector)
    {
      unsigned rowStart, rowEnd, colStart, colEnd;
      getRowColStartEnd(coordA, nDist, aLipidStart, nLipidRadius, rowStart,
                        colStart, rowEnd, colEnd, widthVector, lengthVector);
      for(unsigned i(rowStart); i != rowEnd; ++i)
        {
          for(unsigned j(colStart); j != colEnd; ++j)
            {
              unsigned pairCoord(i*lipCols+j);
              if(theSpecies[aSpeciesID
                 ]->isValidProductPairCoord(i, pairCoord, myCoords))
                { 
                  aProductPairOffsets.push_back(long(pairCoord)-long(coordA));
                }
            }
        }
    }
  bool isValidProductPairCoord(const unsigned row, const unsigned coord,
                               const std::vector<unsigned>& pairCoords)  
    {
      for(unsigned i(0); i != theOffsets[row%2][0].size(); ++i)
        {
          const unsigned coordB(theOffsets[row%2][0][i]+coord);
          if(std::find(pairCoords.begin(), pairCoords.end(), coordB) !=
             pairCoords.end())
            {
              return false;
            }
        }
      for(unsigned i(0); i != theOffsets[row%2][0].size(); ++i)
        {
          const unsigned coordB(theOffsets[row%2][0][i]+coord);
          Voxel& aVoxel(theLattice[coordB+lipStartCoord]);
          for(unsigned j(0); j != aVoxel.diffuseSize; ++j)
            {
              const unsigned coordC(aVoxel.adjoiningCoords[j]-lipStartCoord);
              if(std::find(pairCoords.begin(), pairCoords.end(), coordC) !=
                 pairCoords.end())
                {
                  return true;
                }
            }
        }
      return false;
    }
  void setIntersectOffsets(const std::vector<std::vector<int> > anAdjoinOffsets,
                           const std::vector<int> aRowOffsets,
                           Species* aLipid, const Point& aLipidStart, 
                           const unsigned aVacantRows,
                           const unsigned aVacantCols,
                           const double nLipidRadius,
                           const double aSubunitAngle, Point& aSurfaceNormal,
                           Point& widthVector, Point& lengthVector)
    {
      theWidthVector = widthVector;
      theLengthVector = lengthVector;
      setAdjoinOffsets(anAdjoinOffsets, aRowOffsets);
      double nDist((aLipid->getMoleculeRadius()+theMoleculeRadius)/
                   (2*theVoxelRadius));
      if(aLipid->getMoleculeRadius() >= theMoleculeRadius)
        {
          nDist = 1.5*theMoleculeRadius/(2*theVoxelRadius);
        }
      theRegLatticeCoord = lipRows/2*lipCols+lipCols/2;
      theOffsets.resize(2);
      unsigned rowA(aVacantRows/2);
      unsigned rowB(rowA+1);
      if(rowB >= aVacantRows && rowA > 0)
        {
          rowB = rowA-1;
        }
      else if(!rowA)
        {
          rowB = 0;
        }
      unsigned aCol(aVacantCols/2);
      unsigned coordA(rowA*aVacantCols+aCol);
      unsigned coordB(rowB*aVacantCols+aCol);
      theTarOffsets.resize(2);
      theSrcOffsets.resize(2);
      theRotOffsets.resize(2);
      setOffsets(rowA, coordA, nDist, aLipidStart, nLipidRadius, aSubunitAngle,
                 aSurfaceNormal, widthVector, lengthVector, theOffsets[rowA%2]);
      setOffsets(rowB, coordB, nDist, aLipidStart, nLipidRadius, aSubunitAngle,
                 aSurfaceNormal, widthVector, lengthVector, theOffsets[rowB%2]);
      theTotalLipidSites = theOffsets[0][0].size();
    }
  void setOffsets(const unsigned row, const unsigned coordA, const double nDist,
                  const Point& aLipidStart, const double nLipidRadius,
                  const double aSubunitAngle, Point& aSurfaceNormal,  
                  Point& widthVector, Point& lengthVector,
                  std::vector<std::vector<int> >& refOffsets)
    {
      theRotateSize = 16;
      refOffsets.resize(theRotateSize);
      theTarOffsets[row%2].resize(theRotateSize);
      theSrcOffsets[row%2].resize(theRotateSize);
      theRotOffsets[row%2].resize(theRotateSize);
      double angle(0);
      double incAngle(M_PI*2/theRotateSize);
      for(unsigned i(0); i != theRotateSize; ++i)
        {
          setCoordOffsets(coordA, coordA, nDist, aLipidStart, nLipidRadius,
                          aSubunitAngle, aSurfaceNormal, widthVector, lengthVector,
                          angle, refOffsets[i]);
          setWalkOffsets(row, coordA, nDist, aLipidStart, nLipidRadius,
                         aSubunitAngle, aSurfaceNormal, widthVector, lengthVector,
                         angle, refOffsets[i], theTarOffsets[row%2][i],
                         theSrcOffsets[row%2][i]);
          theRotOffsets[row%2][i].resize(2);
          if(i)
            {
              setDiffOffsets(refOffsets[i-1], refOffsets[i],
                             theRotOffsets[row%2][i-1][1]);
              setDiffOffsets(refOffsets[i], refOffsets[i-1],
                             theRotOffsets[row%2][i][0]);
            }
          angle += incAngle;
        }
      setDiffOffsets(refOffsets[theRotateSize-1], refOffsets[0],
                     theRotOffsets[row%2][theRotateSize-1][1]);
      setDiffOffsets(refOffsets[0], refOffsets[theRotateSize-1],
                     theRotOffsets[row%2][0][0]);
    }
  void setWalkOffsets(const unsigned row, const unsigned coordA,
                      const double nDist, const Point& aLipidStart,
                      const double nLipidRadius, const double aSubunitAngle,
                      Point& aSurfaceNormal, Point& widthVector,
                      Point& lengthVector, const double angle,
                      std::vector<int>& prevOffsets,
                      std::vector<std::vector<int> >& aTarOffsets,
                      std::vector<std::vector<int> >& aSrcOffsets)
    {
      aTarOffsets.resize(theDiffuseSize);
      aSrcOffsets.resize(theDiffuseSize);
      for(unsigned i(0); i != theDiffuseSize; ++i)
        {
          const unsigned coordB(theAdjoinOffsets[row%2][i]+coordA);
          std::vector<int> nextOffsets;
          setCoordOffsets(coordB, coordA, nDist, aLipidStart, nLipidRadius,
                          aSubunitAngle, aSurfaceNormal, widthVector, lengthVector,
                          angle, nextOffsets);
          setDiffOffsets(prevOffsets, nextOffsets, aTarOffsets[i]);
          setDiffOffsets(nextOffsets, prevOffsets, aSrcOffsets[i]);
        }
    }
  void setDiffOffsets(std::vector<int>& srcOffsets,
                      std::vector<int>& tarOffsets,
                      std::vector<int>& aWalkOffsets)
    {
      for(unsigned i(0); i != tarOffsets.size(); ++i)
        {
          if(std::find(srcOffsets.begin(), srcOffsets.end(), tarOffsets[i]) ==
             srcOffsets.end())
            {
              aWalkOffsets.push_back(tarOffsets[i]);
            }
        }
    }
  void setCoordOffsets(const unsigned coordA, const unsigned coordB,
                       const double nDist, const Point& aLipidStart,
                       const double nLipidRadius, const double aSubunitAngle,
                       Point& aSurfaceNormal, Point& widthVector,
                       Point& lengthVector, const double aRotateAngle,
                       std::vector<int>& anIntersectOffsets)
    {
      std::vector<unsigned> anIntersectLipids;
      getIntersectLipidsRegular(coordA+vacStartCoord, nDist, aLipidStart,
                                nLipidRadius, aSubunitAngle, aSurfaceNormal,
                                widthVector, lengthVector, aRotateAngle,
                                anIntersectLipids);
      anIntersectOffsets.resize(0);
      for(unsigned i(0); i != anIntersectLipids.size(); ++i)
        {
          anIntersectOffsets.push_back(long(anIntersectLipids[i])-long(coordB));
        }
    }
  void getIntersectLipidsRegular(const unsigned coordA, const double nDist,
                                 const Point& aLipidStart,
                                 const double nLipidRadius,
                                 const double aSubunitAngle,
                                 Point& aSurfaceNormal,
                                 Point& widthVector,
                                 Point& lengthVector,
                                 const double aRotateAngle,
                                 std::vector<unsigned>& anIntersectLipids)
    {
      Point& pointA(*theLattice[coordA].point);
      //The rectangle shape of the rotating subunit on the 2D surface
      //(assuming the molecule is in rectangular shape)
      Point maxRect(disp(pointA, widthVector, nDist*sin(aSubunitAngle)));
      Point minRect(disp(pointA, widthVector, -nDist*sin(aSubunitAngle)));
      maxRect = disp(maxRect, lengthVector, nDist*cos(aSubunitAngle));
      minRect = disp(minRect, lengthVector, -nDist*cos(aSubunitAngle));
      /*
      double maxRectY(pointA.y+nDist*sin(aSubunitAngle));
      double minRectY(pointA.y-nDist*sin(aSubunitAngle));
      double maxRectZ(pointA.z+nDist*cos(aSubunitAngle));
      double minRectZ(pointA.z-nDist*cos(aSubunitAngle));
      */
      double minWidth(point2lineDist(pointA, lengthVector, aLipidStart));
      double minLength(point2lineDist(pointA, widthVector, aLipidStart));
      double maxWidth(minWidth+nDist*2);
      double maxLength(minLength+nDist*2);
      minWidth = std::max(0.0, minWidth-nDist*2);
      minLength = std::max(0.0, minLength-nDist*2);
      unsigned rowStart((unsigned)std::max(minWidth/(nLipidRadius*sqrt(3))-1,
                                           0.0));
      unsigned rowEnd((unsigned)std::min(maxWidth/(nLipidRadius*sqrt(3))+1,
                                         double(lipRows)));
      unsigned colStart((unsigned)std::max(minLength/(nLipidRadius*2)-1, 0.0));
      unsigned colEnd((unsigned)std::min(maxLength/(nLipidRadius*2)+1,
                                         double(lipCols)));
      for(unsigned i(rowStart); i != rowEnd; ++i)
        {
          for(unsigned j(colStart); j != colEnd; ++j)
            {
              unsigned coord(i*lipCols+j);
              Point pointB(*theLattice[coord+lipStartCoord].point);
              rotatePointAlongVector(pointB, pointA, aSurfaceNormal,
                                     aRotateAngle);
              if(aSubunitAngle)
                {
                  if(pointB.y >= minRect.y && pointB.y <= maxRect.y &&
                     pointB.x >= minRect.x && pointB.x <= maxRect.x &&
                     pointB.z >= minRect.z && pointB.z <= maxRect.z)
                    {
                      anIntersectLipids.push_back(coord);
                    }
                }
              else if(distance(pointA, pointB) < nDist)
                {
                  anIntersectLipids.push_back(coord);
                }
            }
        }
    }
  void setMultiMultiReactant(unsigned anID)
    {
      isMultiMultiReactive = true;
      isMultiMultiReactantID[anID] = true;
      for(unsigned i(0); i != isMultiscaleBoundID.size(); ++i)
        {
          if(isMultiscaleBoundID[i])
            {
              theSpecies[anID]->setMultiMultiReactant(i);
            }
        }
    }
  void setMultiscaleBindIDs(unsigned subID, unsigned prodID)
    {
      isMultiscaleBoundID[prodID] = true;
      isMultiscaleBinderID[subID] = true;
      if(isMultiMultiReactive)
        {
          for(unsigned i(0); i != isMultiMultiReactantID.size(); ++i)
            {
              if(isMultiMultiReactantID[i])
                {
                  theSpecies[i]->setMultiMultiReactant(prodID);
                }
            }
        }
    }
  void setMultiscaleUnbindIDs(unsigned subID, unsigned prodID)
    {
      isMultiscaleBoundID[subID] = true;
      theMultiscaleUnbindIDs[subID] = prodID;
    }
  //Get the fraction of number of nanoscopic molecules (anID) within the
  //multiscale molecule (index):
  double getMultiscaleBoundFraction(unsigned index, unsigned anID)
    {
      double fraction(0);
      if(isMultiscale)
        {
          const unsigned coordA(getCoord(index)-vacStartCoord);
          //Need to optimize this, which is a big bottleneck:
          if(isRegularLattice)
            {
              const int rowA(coordA/lipCols);
              const std::vector<int>& anOffsetsA(theOffsets[rowA%2][
                                                 theTags[index].rotIndex]);
              unsigned size(0);
              for(unsigned i(0); i != anOffsetsA.size(); ++i)
                {
                  const int offsetRow((anOffsetsA[i]+theRegLatticeCoord)/
                                      lipCols-theRegLatticeCoord/lipCols);
                  int coord(coordA+anOffsetsA[i]);
                  if(isInLattice(coord, offsetRow+rowA))
                    {
                      ++size;
                      if(getID(theLattice[coord+lipStartCoord]) == anID)
                        {
                          fraction += 1;
                        }
                    }
                }
              fraction /= size;
            }
          else
            {
              for(unsigned i(0); i != theIntersectLipids[coordA].size(); ++i)
                {
                  unsigned aCoord(theIntersectLipids[coordA][i]+lipStartCoord);
                  if(getID(theLattice[aCoord]) == anID)
                    {
                      fraction += 1;
                    }
                }
              fraction /= theIntersectLipids[coordA].size();
            }
        }
      return fraction;
    }
  unsigned getPopulatableSize()
    {
      if(isMultiscale)
        {
          if(!getIsPopulated())
            {
              cout << "The multiscale species:" << 
                getVariable()->getFullID().asString() << " has not yet " <<
                "been populated, but it being populated on." << std::endl;
            }
          thePopulatableCoords.resize(0);
          for(unsigned i(0); i != theMoleculeSize; ++i)
            {
              const unsigned coordA(getCoord(i)-vacStartCoord);
              if(isRegularLattice)
                {
                  const int rowA(coordA/lipCols);
                  const std::vector<int>& anOffsetsA(theOffsets[
                                               rowA%2][theTags[i].rotIndex]);
                  for(unsigned j(0); j != anOffsetsA.size(); ++j)
                    {
                      const int offsetRow((anOffsetsA[j]+theRegLatticeCoord)/
                                          lipCols-theRegLatticeCoord/lipCols);
                      int coord(coordA+anOffsetsA[j]);
                      if(isInLattice(coord, offsetRow+rowA))
                        {
                          coord += lipStartCoord;
                          if(getID(theLattice[coord]) == theID)
                            {
                              thePopulatableCoords.push_back(coord);
                            }
                        }
                    }
                }
              else
                {
                  for(unsigned j(0); j != theIntersectLipids[coordA].size();
                      ++j)
                    {
                      unsigned aCoord(theIntersectLipids[coordA][j]+
                                      lipStartCoord);
                      if(getID(theLattice[aCoord]) == theID)
                        {
                          thePopulatableCoords.push_back(aCoord);
                        }
                    }
                }
            }
          return thePopulatableCoords.size();
        }
      //Required by populate dense because some comp vacant voxels would
      //have become interface species voxels and no longer populatable:
      else if(isVacant)
        {
          unsigned aSize(0);
          for(unsigned i(0); i != theMoleculeSize; ++i)
            {
              if(getID(theMolecules[i]) == theID)
                {
                  ++aSize;
                }
            }
          return aSize;
        }
      return theMoleculeSize;
    }
  Voxel* getRandomValidMolecule()
    {
      if(isMultiscale)
        {
          unsigned index(0);
          do
            {
              index = theRng.Integer(thePopulatableCoords.size());
            }
          while(getID(theLattice[thePopulatableCoords[index]]) != theID);
          return &theLattice[thePopulatableCoords[index]];
        }
      else
        {
          const unsigned anIndex(getRandomValidIndex());
          if(anIndex == theMoleculeSize)
            {
              return NULL;
            }
          return theMolecules[anIndex];
        }
    }
  unsigned getPopulatableCoord(unsigned index)
    {
      if(isMultiscale)
        {
          return thePopulatableCoords[index];
        }
      return getCoord(index);
    }
  void setMultiscaleVacantSpecies(Species* aSpecies)
    {
      theMultiscaleVacantSpecies = aSpecies;
    }
  const unsigned getVacantIdx()
    {
      return theVacantIdx;
    }
  //Can aVoxel be populated by this species:
  bool isPopulatable(Voxel* aVoxel)
    {
      if(isMultiscale)
        {
          if(isMultiIntersectCoord(aVoxel->coord))
            {
              return false;
            }
        }
      else if(getID(aVoxel) != theVacantID)
        {
          return false;
        }
      return true;
    }
  //Can aVoxel of this species replaced by aSpecies:
  //TODO: simplify this:
  bool isReplaceable(Voxel* aVoxel, Species* aSpecies)
    {
      if(getComp() != aSpecies->getComp() &&
         theID != aSpecies->getVacantID())
        {
          return false;
        }
      if(aSpecies->getIsMultiscale())
        {
          if(aSpecies->isMultiIntersectCoord(aVoxel->coord))
            {
              return false;
            }
        }
      if(getVacantID() != aSpecies->getID() && 
         aSpecies->getVacantID() != theID && 
         getVacantID() != aSpecies->getVacantID())
        {
          return false;
        }
      return true;
    }
  bool isReplaceable(Species* aSpecies)
    {
      if(getComp() != aSpecies->getComp() &&
         theID != aSpecies->getVacantID())
        {
          return false;
        }
      if(getVacantID() != aSpecies->getID() && 
         aSpecies->getVacantID() != theID && 
         getVacantID() != aSpecies->getVacantID())
        {
          return false;
        }
      return true;
    }
  unsigned getMultiscaleStructureSize()
    {
      theMultiscaleStructureCoords.resize(0);
      //TODO: should be able to optimize this if you don't
      //want to use this for bug checking:
      for(unsigned i(0); i != theMultiscaleVacantSpecies->size(); ++i)
        {
          unsigned coord(theMultiscaleVacantSpecies->getCoord(i));
          if(getID(theLattice[coord]) == theID)
            {
              theMultiscaleStructureCoords.push_back(coord);
            }
        }
      return theMultiscaleStructureCoords.size();
    }
  Point& getMultiscaleStructurePoint(unsigned index)
    {
      return *theLattice[theMultiscaleStructureCoords[index]].point;
    }
private:
  bool isCentered;
  bool isCompVacant;
  bool isDeoligomerize;
  bool isDiffusing;
  bool isDiffusiveVacant;
  bool isFixedAdjoins;
  bool isGaussianPopulation;
  bool isInContact;
  bool isInterface;
  bool isMultiMultiReactive;
  bool isMultiscale;
  bool isMultiscaleComp;
  bool isOffLattice;
  bool isOnMultiscale;
  bool isOrigins;
  bool isPeriodic;
  bool isPolymer;
  bool isReactiveVacant;
  bool isRegularLattice;
  bool isSubunitInitialized;
  bool isTag;
  bool isTagged;
  bool isVacant;
  const unsigned short theID;
  int lipCols;
  int lipRows;
  unsigned lipStartCoord;
  unsigned theAdjoiningCoordSize;
  unsigned theDiffuseSize;
  unsigned theDimension;
  unsigned theInitCoordSize;
  unsigned theMoleculeSize;
  unsigned theNullCoord;
  unsigned theNullID;
  unsigned theRegLatticeCoord;
  unsigned theStride;
  unsigned theRotateSize;
  unsigned long theSpeciesCollisionCnt;
  unsigned theVacantIdx;
  unsigned vacCols;
  unsigned vacRows;
  unsigned vacStartCoord;
  unsigned theVacantID;
  unsigned theTotalLipidSites;
  int thePolymerDirectionality;
  double D;
  double nDiffuseRadius;
  double theDiffuseRadius;
  double theDiffusionInterval;
  double theMoleculeRadius;
  double theVoxelRadius;
  double theWalkProbability;
  double theWalkPropensity;
  Point theLengthVector;
  Point theWidthVector;
  RandomLib::Random& theRng;
  Species* theTrailSpecies;
  Species* theVacantSpecies;
  Species* theMultiscaleVacantSpecies;
  Species* theDeoligomerizedProduct;
  Comp* theComp;
  MoleculePopulateProcess* thePopulateProcess;
  SpatiocyteStepper* theStepper;
  Variable* theVariable;
  Tag theNullTag;
  SpatiocyteDebug cout;
  std::vector<int> theRowOffsets;
  std::vector<std::vector<int> > theAdjoinOffsets;
  std::vector<std::vector<std::vector<int> > > theProductPairOffsets;
  std::vector<std::vector<std::vector<int> > > theOffsets;
  std::vector<std::vector<std::vector<std::vector<int> > > > theTarOffsets;
  std::vector<std::vector<std::vector<std::vector<int> > > > theSrcOffsets;
  std::vector<std::vector<std::vector<std::vector<int> > > > theRotOffsets;
  std::vector<bool> isFinalizeReactions;
  std::vector<bool> isMultiscaleBinderID;
  std::vector<bool> isMultiMultiReactantID;
  std::vector<bool> isMultiscaleBoundID;
  std::vector<bool> isProductPair;
  std::vector<unsigned> collisionCnts;
  std::vector<unsigned> theCoords;
  std::vector<unsigned> theMultiscaleUnbindIDs;
  std::vector<unsigned> thePopulatableCoords;
  std::vector<unsigned> theMultiscaleStructureCoords;
  std::vector<Tag> theTags;
  std::vector<double> theBendAngles;
  std::vector<double> theReactionProbabilities;
  std::vector<Voxel*> theMolecules;
  std::vector<Voxel*>* theCompVoxels;
  std::vector<Species*> theDiffusionInfluencedReactantPairs;
  std::vector<Species*> theTaggedSpeciesList;
  std::vector<Species*> theTagSpeciesList;
  std::vector<DiffusionInfluencedReactionProcess*> 
    theDiffusionInfluencedReactions;
  std::vector<SpatiocyteProcess*> theInterruptedProcessesAddMolecule;
  std::vector<SpatiocyteProcess*> theInterruptedProcessesRemoveMolecule;
  std::vector<SpatiocyteProcess*> theInterruptedProcessesEndDiffusion;
  std::vector<Voxel>& theLattice;
  std::vector<std::vector<unsigned> > theIntersectLipids;
  std::vector<unsigned> theBoundCnts;
  std::vector<Species*>& theSpecies;
  std::vector<unsigned> theIdxList;
  std::vector<unsigned> theBindList;
  std::vector<unsigned> theVacantList;
  std::vector<unsigned> theUnbindList;
  std::vector<std::vector<double> > theInterfaceConsts;
};

}

#endif /* __SpatiocyteSpecies_hpp */
