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


#ifndef __CompartmentProcess_hpp
#define __CompartmentProcess_hpp

#include <sstream>
#include <libecs/SpatiocyteProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>

/*Options---------------
 RegularLattice = 1, only for flat surface with voxels directly accessible by
                     coordinates. For non-multiscale diffusion (one molecule,
                     one voxel), it uses walkRegular() method, which is slower 
                     than RegularLattice=0 (uses the standard walk() method).
 Periodic = 1, set periodic boundary condition at the edges of the 1D or 2D
               compartment.
*/

LIBECS_DM_CLASS(CompartmentProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(CompartmentProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Integer, Autofit);
      PROPERTYSLOT_SET_GET(Integer, BindingDirection);
      PROPERTYSLOT_SET_GET(Integer, Filaments);
      PROPERTYSLOT_SET_GET(Integer, Periodic);
      PROPERTYSLOT_SET_GET(Integer, PlaneXY);
      PROPERTYSLOT_SET_GET(Integer, PlaneXZ);
      PROPERTYSLOT_SET_GET(Integer, PlaneYZ);
      PROPERTYSLOT_SET_GET(Integer, RegularLattice);
      PROPERTYSLOT_SET_GET(Integer, Subunits);
      PROPERTYSLOT_SET_GET(Integer, DissociationDirection);
      PROPERTYSLOT_SET_GET(Integer, Verbose);
      PROPERTYSLOT_SET_GET(Real, DiffuseRadius); //off-lattice voxel radius
      PROPERTYSLOT_SET_GET(Real, Length);
      PROPERTYSLOT_SET_GET(Real, LipidRadius); //radius of lipid voxels
      PROPERTYSLOT_SET_GET(Real, OriginX);
      PROPERTYSLOT_SET_GET(Real, OriginY);
      PROPERTYSLOT_SET_GET(Real, OriginZ);
      PROPERTYSLOT_SET_GET(Real, RotateX);
      PROPERTYSLOT_SET_GET(Real, RotateY);
      PROPERTYSLOT_SET_GET(Real, RotateZ);
      PROPERTYSLOT_SET_GET(Real, SubunitAngle);
      //SubunitRadius is the actual radius of a molecule.
      //For a normal molecule the SubunitRadius = DiffuseRadius.
      //For a multiscale molecule, the SubunitRadius can be larger
      //than the DiffuseRadius:
      PROPERTYSLOT_SET_GET(Real, SubunitRadius);
      PROPERTYSLOT_SET_GET(Real, Width);
    }
  CompartmentProcess():
    isCompartmentalized(false),
    Autofit(1),
    BindingDirection(2), //bidirectional
    Filaments(0),
    LipidCols(0),
    LipidRows(0),
    Periodic(0),
    PlaneXY(1),
    PlaneXZ(0),
    PlaneYZ(0),
    RegularLattice(1),
    Subunits(0),
    DissociationDirection(2), //bidirectional
    Verbose(0),
    theDiffuseSize(6),
    theDimension(1),
    DiffuseRadius(0),
    Length(0),
    LipidRadius(0),
    nVoxelRadius(0.5),
    OriginX(0),
    OriginY(0),
    OriginZ(0),
    RotateX(0),
    RotateY(0),
    RotateZ(0),
    SubunitAngle(0),
    SubunitRadius(0),
    Width(0),
    theLipidSpecies(NULL),
    theVacantSpecies(NULL) {}
  virtual ~CompartmentProcess() {}
  SIMPLE_SET_GET_METHOD(Integer, Autofit);
  SIMPLE_SET_GET_METHOD(Integer, BindingDirection);
  SIMPLE_SET_GET_METHOD(Integer, Filaments);
  SIMPLE_SET_GET_METHOD(Integer, Periodic);
  SIMPLE_SET_GET_METHOD(Integer, PlaneXY);
  SIMPLE_SET_GET_METHOD(Integer, PlaneXZ);
  SIMPLE_SET_GET_METHOD(Integer, PlaneYZ);
  SIMPLE_SET_GET_METHOD(Integer, RegularLattice);
  SIMPLE_SET_GET_METHOD(Integer, Subunits);
  SIMPLE_SET_GET_METHOD(Integer, DissociationDirection);
  SIMPLE_SET_GET_METHOD(Integer, Verbose);
  SIMPLE_SET_GET_METHOD(Real, DiffuseRadius);
  SIMPLE_SET_GET_METHOD(Real, Length);
  SIMPLE_SET_GET_METHOD(Real, LipidRadius);
  SIMPLE_SET_GET_METHOD(Real, OriginX);
  SIMPLE_SET_GET_METHOD(Real, OriginY);
  SIMPLE_SET_GET_METHOD(Real, OriginZ);
  SIMPLE_SET_GET_METHOD(Real, RotateX);
  SIMPLE_SET_GET_METHOD(Real, RotateY);
  SIMPLE_SET_GET_METHOD(Real, RotateZ);
  SIMPLE_SET_GET_METHOD(Real, SubunitAngle);
  SIMPLE_SET_GET_METHOD(Real, SubunitRadius);
  SIMPLE_SET_GET_METHOD(Real, Width);
  virtual void prepreinitialize();
  virtual void initialize();
  virtual void initializeFirst();
  virtual void initializeCompartmentOnce();
  virtual unsigned getLatticeResizeCoord(unsigned);
  virtual void initializeCompartment();
  virtual void printParameters();
  virtual void updateResizedLattice();
  virtual void setCompartmentDimension();
  virtual void initializeVectors();
  virtual void connectFilaments(unsigned, unsigned, unsigned);
  virtual void setSubunitStart();
  virtual void elongateFilaments(Species*, unsigned, unsigned, unsigned,
                                 double);
  virtual void initializeFilaments(Point&, unsigned, unsigned, double, Species*,
                                   unsigned);
  virtual void addSurfaceIntersectInterfaceVoxel(Voxel&, Point&);
  virtual bool isInside(Point&);
  virtual bool isOnAboveSurface(Point&);
  virtual double getDisplacementToSurface(Point&);
  void connectSubunit(unsigned, unsigned, unsigned, unsigned);
  void addInterfaceVoxel(Voxel&);
  void setNearestSubunit(const unsigned, const unsigned);
  void addSortedSubunitInterface(const unsigned, const unsigned, const double,
                                 const unsigned);
  unsigned getNearestInterface(const unsigned, double&);
  void setSubunitInterfaces();
  void addSubunitInterfaces();
  void setNearestSubunitForOrphanInterfaces();
  void setNearestInterfaceForOrphanSubunits();
  void setVacantCompSpeciesProperties();
  void setLipidCompSpeciesProperties();
  void setDiffuseSize(unsigned, unsigned);
  void interfaceSubunits();
  void addFirstInterface();
  void enlistOrphanSubunitInterfaceVoxels();
  virtual void extendInterfacesOverSurface();
  void connectSubunitInterfaceAdjoins();
  void rotateAsParent(Point&);
  void rotate(Point&);
  void addAdjoin(Voxel&, unsigned);
  void setSpeciesIntersectLipids();
  void setSpeciesIntersectVacants();
  void getStartVoxelPoint(Point&, Point&, Point&);
  void setAdjoinOffsets(const unsigned);
  int getCoefficient(Species*);
  void allocateGrid();
  void setGrid(Species*, std::vector<std::vector<unsigned> >&, unsigned);
  bool setSubunitInterfaceVoxel(const unsigned, const double, 
                                const bool isSingle=false);
  Species* coefficient2species(int);
  Voxel* getNearestVoxelToSubunit(const unsigned, double&, const bool);
  Voxel* getNearestVoxelToPoint(Point&, double&, const bool);
  Voxel* getNearestVoxelToSurface(const unsigned, double&, const bool);
  Voxel* addCompVoxel(unsigned, unsigned, Point&, Species*, unsigned, unsigned);
  virtual void removeAdjoinsFromNonBindingSide(Voxel&);
  void setSubunitBindFractions();
  virtual void setCompSubunitBindFractions();
  void initializeConstants();
  void addSubunitBinder(unsigned, unsigned);
protected:
  void setInterfaceConsts();
  bool isThisCompInterface(const unsigned);
  void removeInterfaceCompVoxels();
  bool isDissociationSide(const unsigned);
  bool isBindingSide(const unsigned);
  bool isCompartmentalized;
  unsigned Autofit;
  unsigned BindingDirection;
  unsigned endCoord;
  unsigned Filaments;
  unsigned LipidCols;
  unsigned LipidRows;
  unsigned lipStartCoord;
  unsigned Periodic;
  unsigned PlaneXY;
  unsigned PlaneXZ;
  unsigned PlaneYZ;
  unsigned subStartCoord;
  unsigned RegularLattice;
  unsigned Subunits;
  unsigned DissociationDirection;
  unsigned intSize;
  unsigned intStartIndex;
  unsigned theDiffuseSize;
  unsigned theDimension;
  unsigned vacStartIndex;
  unsigned Verbose;
  double DiffuseRadius;
  double nGridSize;
  double Height;
  double nHeight;
  double Length;
  double lengthDisplace;
  double lengthDisplaceOpp;
  double LipidRadius;
  double nDiffuseRadius;
  double nLength;
  double nLipidRadius;
  double nSubunitRadius;
  double nVoxelRadius;
  double nMaxRadius;
  double nMinRadius;
  double nWidth;
  double OriginX;
  double OriginY;
  double OriginZ;
  double RotateX;
  double RotateY;
  double RotateZ;
  double SubunitAngle;
  double SubunitRadius;
  double surfaceDisplace;
  double VoxelRadius;
  double Width;
  double widthDisplace;
  double widthDisplaceOpp;
  IntPoint parentGrid;
  Point parentOrigin;
  Point parentVectorX;
  Point parentVectorY;
  Point parentVectorZ;
  Point heightEnd;
  Point heightVector;
  Point lengthEnd;
  Point lengthStart;
  Point lengthVector;
  Point lipidStart;
  Point Origin;
  Point subunitStart;
  Point subunitEnd;
  Point surfaceNormal;
  Point widthEnd;
  Point widthVector;
  Comp* theComp;
  Species* theLipidSpecies;
  Species* theInterfaceSpecies;
  Species* theVacantSpecies;
  Variable* theInterfaceVariable;
  Variable* theLipidVariable;
  Variable* theVacantVariable;
  std::vector<double> subunitBindFractions;
  std::vector<std::vector<unsigned> > theVacGrid;
  std::vector<std::vector<unsigned> > theLipGrid;
  std::vector<Point> thePoints;
  std::vector<Species*> theLipidCompSpecies;
  std::vector<Species*> theVacantCompSpecies;
  std::vector<std::vector<unsigned> > subunitInterfaces;
  std::vector<std::vector<unsigned> > subunitBinders;
  //Contains subunits connected by the interface:
  std::vector<std::vector<unsigned> > interfaceSubs;
  //Contains comp vacants that can bind to the subunits connected by the
  //interface:
  std::vector<std::vector<unsigned> > interfaceBinders;
  std::vector<std::vector<double> > subunitInterfaceDists;
  std::vector<std::vector<int> > theAdjoinOffsets;
  std::vector<int> theRowOffsets;
};

#endif /* __CompartmentProcess_hpp */



