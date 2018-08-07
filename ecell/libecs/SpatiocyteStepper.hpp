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


#ifndef __SpatiocyteStepper_hpp
#define __SpatiocyteStepper_hpp

#include <libecs/Stepper.hpp>
#include <RandomLib/Random.hpp>
#include <random_simd/ranvec1.h>
#include <libecs/SpatiocyteCommon.hpp>
#include <libecs/SpatiocyteDebug.hpp>

namespace libecs
{

LIBECS_DM_CLASS(SpatiocyteStepper, Stepper)
{ 
public: 
  LIBECS_DM_OBJECT(SpatiocyteStepper, Stepper)
    {
      INHERIT_PROPERTIES(Stepper);
      PROPERTYSLOT_SET_GET(Integer, Seed);
      PROPERTYSLOT_SET_GET(Integer, LatticeType);
      PROPERTYSLOT_SET_GET(Integer, SearchVacant);
      PROPERTYSLOT_SET_GET(Integer, DebugLevel);
      PROPERTYSLOT_SET_GET(Integer, RemoveSurfaceBias);
      PROPERTYSLOT_SET_GET(Real, VoxelRadius);
    }
  SIMPLE_SET_GET_METHOD(Integer, Seed); 
  SIMPLE_SET_GET_METHOD(Integer, LatticeType); 
  SIMPLE_SET_GET_METHOD(Integer, SearchVacant); 
  SIMPLE_SET_GET_METHOD(Integer, DebugLevel); 
  SIMPLE_SET_GET_METHOD(Integer, RemoveSurfaceBias); 
  SIMPLE_SET_GET_METHOD(Real, VoxelRadius); 
  SpatiocyteStepper():
    isInitialized(false),
    isPeriodicEdge(false),
    SearchVacant(false),
    RemoveSurfaceBias(false),
    Seed(1),
    DebugLevel(1),
    LatticeType(HCP_LATTICE),
    theMoleculeID(0),
    nVoxelRadius(0.5),
    VoxelRadius(10e-9),
    theRanSimd(1),
    cout(std::cout) {}
  virtual ~SpatiocyteStepper() {}
  virtual void initialize();
  /*
  bool isDependentOn(const Stepper*)
    {
      return false;
    }
    */
  // need to check interrupt when we suddenly stop the simulation, do we
  // need to update the priority queue?
  virtual void interrupt(double);
  virtual void interruptAllProcesses(const double);
  virtual void step();
  void createSpecies();
  Species* addSpecies(Variable*);
  Species* createSpecies(System*, String);
  static Variable* createVariable(const String&, const System*, Model*);
  static Process* createProcess(const String&, const String&, const System*,
                                Model*);
  Species* getSpecies(Variable*);
  std::vector<Species*> getSpecies();
  Point coord2point(unsigned);
  bool setSurfaceDiffuseSize(const unsigned, const unsigned, const unsigned,
                             const unsigned, const bool, const unsigned);
  Species* id2species(unsigned short);
  Comp* id2Comp(unsigned short);
  void coord2global(unsigned, unsigned&, unsigned&, unsigned&);
  void point2global(Point, unsigned&, unsigned&, unsigned&);
  Comp* system2Comp(System*);
  bool isBoundaryCoord(unsigned, unsigned);
  unsigned getPeriodicCoord(unsigned, unsigned, Origin*);
  unsigned global2coord(unsigned, unsigned, unsigned);
  Point getPeriodicPoint(unsigned, unsigned, Origin*);
  void checkLattice();
  void checkSpecies();
  void setPeriodicEdge();
  void reset(int);
  unsigned getNewMoleculeID();
  unsigned getRowSize();
  unsigned getLayerSize();
  unsigned getColSize();
  unsigned getLatticeSize();
  unsigned short getNullID();
  Point getCenterPoint();
  double getNormalizedVoxelRadius();
  unsigned point2coord(Point&);
  std::vector<Comp*> const& getComps() const;
  Species* variable2species(const Variable*) const;
  void rotate(Comp*, Point&, int sign=1);
  void rotateX(double, Point*, int sign=1);
  void rotateY(double, Point*, int sign=1);
  void rotateZ(double, Point*, int sign=1);
  bool isPeriodicEdgeCoord(unsigned, Comp*);
  bool isRemovableEdgeCoord(unsigned, Comp*);
  double getRowLength();
  double getColLength();
  double getLayerLength();
  double getMinLatticeSpace();
  void updateSpecies();
  void finalizeSpecies();
  void resetVariables();
  void enlistModelVariables();
  unsigned getStartCoord();
  unsigned getID(const Voxel&) const;
  unsigned getID(const Voxel*) const;
  virtual GET_METHOD(Real, TimeScale) { return 0; }
  Voxel* getVoxel(const unsigned int& aCoord) { return &theLattice[aCoord]; }
  void addInterruptedProcess(SpatiocyteProcess*);
  static Variable* getVariable(System*, String const&);
  SpatiocyteProcess* getCurrentProcess();
  Comp* getRootComp();
private:
  void addSurfaceAdjoins(const unsigned, const Comp*);
  void interruptProcesses(const double);
  void setCompsCenterPoint();
  void setCuboidCompsCenterPoint();
  void setIntersectingCompartmentList();
  void setIntersectingParent();
  void setIntersectingPeers();
  void printProcessParameters();
  void printSpeciesIDs();
  void checkSurfaceComp();
  void shuffleAdjoiningCoords();
  void setLatticeProperties();
  void checkModel();
  void resizeProcessLattice();
  void initPriorityQueue();
  void initializeFirst();
  void initializeSecond();
  void initializeCompartments();
  void initializeCompartmentsOnce();
  void initializeThird();
  void initializeBeforePopulate();
  void initializeFourth();
  void initializeFifth();
  void initializeLastOnce();
  void storeSimulationParameters();
  void setSystemSize(System*, double);
  void printSimulationParameters();
  void setCompProperties();
  void setCompAttributes(System*, Comp*);
  void initSpecies();
  void readjustSurfaceBoundarySizes();
  void constructLattice();
  void compartmentalizeLattice();
  void concatenatePeriodicSurfaces();
  void registerComps();
  void setCompsProperties();
  void setCompVoxelProperties();
  void populateComps();
  void broadcastLatticeProperties();
  void clearComps();
  void clearComp(Comp*);
  void populateComp(Comp*);
  void populateSpeciesDense(std::vector<Species*>&, Species*, unsigned,
                            unsigned);
  void populateSpeciesSparse(std::vector<Species*>&, Species*, unsigned);
  void registerCompSpecies(Comp*);
  void setCompProperties(Comp*);
  void setDiffusiveComp(Comp*);
  void setCompCenterPoint(Comp*);
  void setLineVoxelProperties(Comp*);
  void setLineCompProperties(Comp*);
  void setSurfaceVoxelProperties(Comp*);
  void setSurfaceCompProperties(Comp*);
  void setVolumeCompProperties(Comp*);
  void concatenateVoxel(Voxel&, unsigned, unsigned, unsigned);
  void concatenateLayers(Voxel&, unsigned, unsigned, unsigned, unsigned);
  void concatenateRows(Voxel&, unsigned, unsigned, unsigned, unsigned);
  void concatenateCols(Voxel&, unsigned, unsigned, unsigned, unsigned);
  void replaceVoxel(unsigned, unsigned);
  void replaceUniVoxel(unsigned, unsigned);
  void setMinMaxDimensions(unsigned, Comp*);
  void setMinMaxDimensionsSurfaceSub(Comp*);
  bool isInsideCoord(unsigned, Comp*, double);
  bool isSurfaceVoxel(Voxel&, unsigned, Comp*);
  bool isLineVoxel(Voxel&, unsigned, Comp*);
  bool isEnclosedSurfaceVoxel(Voxel&, unsigned, Comp*);
  bool isEnclosedRootSurfaceVoxel(Voxel&, unsigned, Comp*, Comp*);
  bool isPeerCoord(unsigned, Comp*);
  bool isLowerPeerCoord(unsigned, Comp*);
  bool isRootSurfaceVoxel(Voxel&, unsigned, Comp*);
  bool isVoxelTouchingComp(Voxel&, Comp*);
  bool isParentSurfaceVoxel(Voxel&, unsigned, Comp*);
  bool compartmentalizeVoxel(unsigned, Comp*);
  double getCuboidSpecArea(Comp*);
  unsigned coord2row(unsigned);
  unsigned coord2col(unsigned);
  unsigned coord2layer(unsigned);
  Comp* registerComp(System*, std::vector<Comp*>*);
  void rotateCompartment(Comp*);
private:
  bool isInitialized;
  bool isPeriodicEdge;
  bool SearchVacant;
  bool RemoveSurfaceBias;
  unsigned short theNullID;
  unsigned Seed;
  unsigned DebugLevel;
  unsigned LatticeType; 
  unsigned theAdjoiningCoordSize;
  unsigned theBioSpeciesSize;
  unsigned theCellShape;
  unsigned theColSize;
  unsigned theLayerSize;
  unsigned theNullCoord;
  unsigned theRowSize;
  unsigned theStride;
  unsigned theMoleculeID;
  double nVoxelRadius;
  double theHCPl;
  double theHCPx;
  double theHCPy;
  double VoxelRadius; //r_v
  Point theCenterPoint;
  ProcessPriorityQueue thePriorityQueue; 
  SpatiocyteDebug cout;
  std::vector<unsigned> initVariableValues;
  std::vector<Species*>::iterator variable2ispecies(Variable*);
  std::vector<Species*> theSpecies;
  std::vector<Variable*> theVariables;
  std::vector<Comp*> theComps;
  std::vector<Voxel> theLattice;
  std::vector<Process*> theExternInterruptedProcesses;
  RandomLib::Random theRan;
  Ranvec1 theRanSimd;
  std::vector<SpatiocyteProcess*> theInterruptedProcesses;
  std::vector<SpatiocyteProcess*> theSpatiocyteProcesses;
};

}

#endif /* __SpatiocyteStepper_hpp */

