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


#ifndef __FilamentProcess_hpp
#define __FilamentProcess_hpp

#include <sstream>
#include <libecs/CompartmentProcess.hpp>

namespace libecs
{

LIBECS_DM_CLASS(FilamentProcess, CompartmentProcess)
{ 
public:
  LIBECS_DM_OBJECT(FilamentProcess, Process)
    {
      INHERIT_PROPERTIES(CompartmentProcess);
      PROPERTYSLOT_SET_GET(Integer, LineX);
      PROPERTYSLOT_SET_GET(Integer, LineY);
      PROPERTYSLOT_SET_GET(Integer, LineZ);
    }
  FilamentProcess():
    LineX(1),
    LineY(0),
    LineZ(0),
    Radius(0),
    theMinusSpecies(NULL),
    thePlusSpecies(NULL)
  {
    Filaments = 1;
    //Both FilamentProcess and MicrotubuleProcess are not RegularLattice
    //because each subunit may have different diffuseSize
    RegularLattice = 0;
    Autofit = 0;
  }
  virtual ~FilamentProcess() {}
  SIMPLE_SET_GET_METHOD(Integer, LineX);
  SIMPLE_SET_GET_METHOD(Integer, LineY);
  SIMPLE_SET_GET_METHOD(Integer, LineZ);
  virtual void initialize();
  virtual void initializeFirst();
  virtual unsigned getLatticeResizeCoord(unsigned);
  virtual void setCompartmentDimension();
  virtual void initializeVectors();
  virtual void initializeFilaments(Point&, unsigned, unsigned, double, Species*,
                                   unsigned);
  virtual void initializeCompartment();
  virtual void setSubunitStart();
  virtual void connectFilaments(unsigned, unsigned, unsigned);
  virtual void elongateFilaments(Species*, unsigned, unsigned, unsigned,
                                 double);
  virtual bool isInside(Point&);
  virtual bool isOnAboveSurface(Point&);
  virtual void extendInterfacesOverSurface();
  virtual double getDisplacementToSurface(Point&);
  //virtual void addLineIntersectInterfaceVoxel(Voxel&, Point&, const bool, const bool);
  void connectTrailSubunits(unsigned, unsigned, unsigned);
  void setTrailSize(unsigned, unsigned);
  void extendInterfacesOverFilamentSurface(const bool);
  void populateMinusPlusSpecies();
  virtual void setCompSubunitBindFractions();
protected:
  unsigned getAdjoiningInterfaceCnt(Voxel&);
  bool getFilamentAdjoin(Voxel*, const bool, const unsigned, const unsigned,
               const double, const double, Point&, double&, double&, Voxel**);
  bool isAdjoin(Voxel*, Voxel*);
  double getMinDistanceFromLineEnd(Point&);
  unsigned LineX;
  unsigned LineY;
  unsigned LineZ;
  double heightDisplace;
  double nRadius;
  double Radius;
  Point Minus; //Minus end
  Point Plus; //Plus end
  Species* theMinusSpecies;
  Species* thePlusSpecies;
  std::vector<Voxel*> theMinusVoxels;
  std::vector<Voxel*> thePlusVoxels;
};

}

#endif /* __FilamentProcess_hpp */




