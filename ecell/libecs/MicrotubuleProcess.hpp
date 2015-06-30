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


#ifndef __MicrotubuleProcess_hpp
#define __MicrotubuleProcess_hpp

#include <sstream>
#include <libecs/FilamentProcess.hpp>

namespace libecs
{

LIBECS_DM_CLASS(MicrotubuleProcess, FilamentProcess)
{ 
public:
  LIBECS_DM_OBJECT(MicrotubuleProcess, Process)
    {
      INHERIT_PROPERTIES(FilamentProcess);
      PROPERTYSLOT_SET_GET(Real, MonomerPitch);
      //Radius is the radius of the microtubule:
      PROPERTYSLOT_SET_GET(Real, Radius);
    }
  MicrotubuleProcess():
    MonomerPitch(4e-9)
  {
    //Both FilamentProcess and MicrotubuleProcess are not RegularLattice
    //because each subunit may have different diffuseSize
    RegularLattice = 0;
    Filaments = 13;
    DiffuseRadius = 8e-9/2;
    SubunitRadius = DiffuseRadius;
    //Only allow kinesins to bind from the outside of the MT surface, not
    //from the cylindrical tube inside MT
    BindingDirection = 0; 
    DissociationDirection = 0; 
  }
  virtual ~MicrotubuleProcess() {}
  SIMPLE_SET_GET_METHOD(Real, MonomerPitch);
  SIMPLE_SET_GET_METHOD(Real, Radius);
  virtual void initialize();
  virtual void initializeVectors();
  virtual void initializeCompartment();
  virtual void initializeFilaments(Point&, unsigned, unsigned, double, Species*,
                                   unsigned);
  virtual void extendInterfacesOverSurface();
  virtual bool isInside(Point&);
  virtual double getDisplacementToSurface(Point&);
protected:
  unsigned getNearestFilament(Point&);
  double MonomerPitch;
  double nMonomerPitch;
  std::vector<Point> theMinusPoints;
  std::vector<Species*> theKinesinSpecies;
};

}

#endif /* __MicrotubuleProcess_hpp */




