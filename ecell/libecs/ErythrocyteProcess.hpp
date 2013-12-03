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


#ifndef __ErythrocyteProcess_hpp
#define __ErythrocyteProcess_hpp

#include <sstream>
#include <libecs/CompartmentProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>  

namespace libecs
{

LIBECS_DM_CLASS(ErythrocyteProcess, CompartmentProcess)  
{
public:
  LIBECS_DM_OBJECT(ErythrocyteProcess, Process)
    {
      INHERIT_PROPERTIES(CompartmentProcess);
      PROPERTYSLOT_SET_GET(Real, EdgeLength);
    }
  ErythrocyteProcess():
    EdgeLength(60e-9),
    theEdgeSpecies(NULL),
    theVertexSpecies(NULL) {}
  virtual ~ErythrocyteProcess() {}
  SIMPLE_SET_GET_METHOD(Real, EdgeLength);
  virtual void prepreinitialize();
  virtual void initialize();
  virtual void initializeFirst();
  virtual unsigned getLatticeResizeCoord(unsigned);
  virtual void initializeVectors();
  virtual void initializeFilaments();
  virtual void initializeThird();
  virtual void normalize(Point&);
  virtual bool isInsidePlane(Point&, Point&, Point&);
  virtual bool isOnPlane(Point&, Point&, Point&, unsigned int);
  virtual bool isOnLowerPlanes(Point&, unsigned int, Point&);
  virtual bool isOnUpperPlanes(Point&, unsigned int, Point&);
  virtual unsigned int getIntersectCount(Point&, unsigned int);
  virtual void rotatePointAlongVector(Point&, Point&, Point&, double);
  virtual void printParameters();
protected:
  double EdgeLength;
  double VoxelDiameter;
  double TriangleAltitude;
  Point Y; //Direction vector along the rotated positive y-axis
  Point X; //Direction vector along the rotated positive x-axis
  Point R; //Direction vector along rotated north east
  Point L; //Direction vector along rotated north west
  Point C; //Center point
  Species* theEdgeSpecies;
  Species* theVertexSpecies;
  std::vector<unsigned int> filamentCoords;
  std::vector<Species*> theEdgeCompSpecies;
  std::vector<Species*> theVertexCompSpecies;
};

}

#endif /* __ErythrocyteProcess_hpp */
