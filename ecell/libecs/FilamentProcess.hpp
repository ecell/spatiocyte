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
    }
  FilamentProcess():
    theMinusSpecies(NULL),
    thePlusSpecies(NULL)
  {
    Filaments = 1;
    /*
    SurfaceDirection = 0;
    Autofit = 0;
    Subunits = 1;
    RegularLattice = 0;
    */
  }
  virtual ~FilamentProcess() {}
  virtual void prepreinitialize();
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
  virtual void addSurfaceIntersectInterfaceVoxel(Voxel&, Point&);
  virtual bool isInside(Point&);
  virtual bool isOnAboveSurface(Point&);
  virtual double getDistanceToSurface(Point&);
  void connectTrailSubunits(unsigned, unsigned, unsigned);
  void setTrailSize(unsigned, unsigned);
protected:
  double nRadius;
  double Radius;
  Point Minus; //Minus end
  Point Plus; //Plus end
  Species* theMinusSpecies;
  Species* thePlusSpecies;
  std::vector<Species*> theBindingSpecies;
};

}

#endif /* __FilamentProcess_hpp */




