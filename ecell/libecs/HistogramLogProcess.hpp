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


#ifndef __HistogramLogProcess_hpp
#define __HistogramLogProcess_hpp

#include <fstream>
#include <math.h>
#include <libecs/IteratingLogProcess.hpp>
#include <libecs/SpatiocyteSpecies.hpp>

namespace libecs
{

LIBECS_DM_CLASS(HistogramLogProcess, IteratingLogProcess)
{ 
public:
  LIBECS_DM_OBJECT(HistogramLogProcess, Process)
    {
      INHERIT_PROPERTIES(IteratingLogProcess);
      PROPERTYSLOT_SET_GET(Integer, Bins);
      PROPERTYSLOT_SET_GET(Integer, StartBin);
      PROPERTYSLOT_SET_GET(Integer, Density);
      PROPERTYSLOT_SET_GET(Real, Radius);
      PROPERTYSLOT_SET_GET(Real, InnerRadius);
      PROPERTYSLOT_SET_GET(Real, Length);
      PROPERTYSLOT_SET_GET(Real, OriginX);
      PROPERTYSLOT_SET_GET(Real, OriginY);
      PROPERTYSLOT_SET_GET(Real, OriginZ);
      PROPERTYSLOT_SET_GET(Real, OuterRadius);
      PROPERTYSLOT_SET_GET(Real, RadialHeight);
      PROPERTYSLOT_SET_GET(Real, RotateX);
      PROPERTYSLOT_SET_GET(Real, RotateY);
      PROPERTYSLOT_SET_GET(Real, RotateZ);
    }
  SIMPLE_SET_GET_METHOD(Integer, Bins);
  SIMPLE_SET_GET_METHOD(Integer, StartBin);
  SIMPLE_SET_GET_METHOD(Integer, Density);
  SIMPLE_SET_GET_METHOD(Real, InnerRadius);
  SIMPLE_SET_GET_METHOD(Real, Radius);
  SIMPLE_SET_GET_METHOD(Real, Length);
  SIMPLE_SET_GET_METHOD(Real, OriginX);
  SIMPLE_SET_GET_METHOD(Real, OriginY);
  SIMPLE_SET_GET_METHOD(Real, OriginZ);
  SIMPLE_SET_GET_METHOD(Real, OuterRadius);
  SIMPLE_SET_GET_METHOD(Real, RadialHeight);
  SIMPLE_SET_GET_METHOD(Real, RotateX);
  SIMPLE_SET_GET_METHOD(Real, RotateY);
  SIMPLE_SET_GET_METHOD(Real, RotateZ);
  HistogramLogProcess():
    Density(1),
    Bins(1),
    StartBin(0),
    InnerRadius(0),
    OriginX(0),
    OriginY(0),
    OriginZ(0),
    OuterRadius(libecs::INF),
    RadialHeight(0),
    Radius(12.5e-9),
    RotateX(0),
    RotateY(0),
    RotateZ(0),
    theCenterSpecies(NULL),
    theFirstBinSpecies(NULL),
    theLastBinSpecies(NULL),
    theMarkerSpecies(NULL)
  {
    FileName = "HistogramLog.csv";
    LogStart = 1e-8;
  }
  virtual ~HistogramLogProcess() {}
  virtual void initialize();
  virtual void initializeLastOnce();
  virtual void saveBackup();
  virtual void logValues();
  virtual void logCollision();
  virtual void logDensity();
  virtual void initLogValues();
  virtual void saveFileHeader(std::ofstream&);
  virtual void saveFileData(std::ofstream&, const unsigned);
  virtual void saveATimePoint(std::ofstream&, const double, const unsigned,
                              const unsigned);
  void initializeVectors();
  void setVacantSizes();
  bool isInside(unsigned int&, Point);
  bool isInsideLength(unsigned int&, Point);
  bool isInsideRadial(unsigned int&, Point);
protected:
  unsigned Density;
  unsigned Bins;
  unsigned StartBin;
  double binInterval;
  double InnerRadius;
  double Length;
  double nLength;
  double nHeight;
  double Width;
  double nWidth;
  double nRadius;
  double OriginX;
  double OriginY;
  double OriginZ;
  double OuterRadius;
  double RadialHeight;
  double Radius;
  double RotateX;
  double RotateY;
  double RotateZ;
  double VoxelDiameter;
  Point lengthVector; 
  Point heightVector; 
  Point widthVector; 
  Point MinusLength; 
  Point MinusHeight; 
  Point MinusWidth; 
  Point SuperOrigin;
  Point CompOrigin;
  Species* theCenterSpecies;
  Species* theFirstBinSpecies;
  Species* theLastBinSpecies;
  Species* theMarkerSpecies;
  std::vector<std::vector<std::vector<double> > > theLogValues;
  std::vector<std::vector<double> > theVacantSizes;
};

}

#endif /* __HistogramLogProcess_hpp */
