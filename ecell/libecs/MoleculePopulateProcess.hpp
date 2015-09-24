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


#ifndef __MoleculePopulateProcess_hpp
#define __MoleculePopulateProcess_hpp

#include <libecs/SpatiocyteProcess.hpp>

namespace libecs
{

LIBECS_DM_CLASS(MoleculePopulateProcess, SpatiocyteProcess)
{ 
public:
  LIBECS_DM_OBJECT(MoleculePopulateProcess, Process)
    {
      INHERIT_PROPERTIES(Process);
      PROPERTYSLOT_SET_GET(Real, OriginX);
      PROPERTYSLOT_SET_GET(Real, OriginY);
      PROPERTYSLOT_SET_GET(Real, OriginZ);
      PROPERTYSLOT_SET_GET(Real, GaussianSigma);
      PROPERTYSLOT_SET_GET(Real, ResetTime);
      PROPERTYSLOT_SET_GET(Real, StartTime);
      PROPERTYSLOT_SET_GET(Real, UniformLengthX);
      PROPERTYSLOT_SET_GET(Real, UniformLengthY);
      PROPERTYSLOT_SET_GET(Real, UniformLengthZ);
      PROPERTYSLOT_SET_GET(Real, UniformLength);
      PROPERTYSLOT_SET_GET(Real, UniformWidth);
      PROPERTYSLOT_SET_GET(Real, UniformHeight);
      PROPERTYSLOT_SET_GET(Real, UniformRadiusWidth); //Inner radius of 
                                                      //circular population
      PROPERTYSLOT_SET_GET(Real, UniformRadiusXY); //For circular population
      PROPERTYSLOT_SET_GET(Real, UniformRadiusXZ); //For circular population
      PROPERTYSLOT_SET_GET(Real, UniformRadiusYZ); //For circular population
      //Populate along the length of the compartment in bins given as
      //fraction of total available vacant voxels in the bin 
      PROPERTYSLOT_SET_GET(Polymorph, LengthBinFractions);
    }
  MoleculePopulateProcess():
    GaussianSigma(0),
    OriginX(0),
    OriginY(0),
    OriginZ(0),
    ResetTime(libecs::INF),
    StartTime(0),
    UniformLength(1),
    UniformWidth(1),
    UniformHeight(1),
    UniformRadiusWidth(-1),
    UniformRadiusXY(0),
    UniformRadiusXZ(0),
    UniformRadiusYZ(0),
    UniformLengthX(1),
    UniformLengthY(1),
    UniformLengthZ(1) {}
  virtual ~MoleculePopulateProcess() {}
  SIMPLE_SET_GET_METHOD(Real, OriginX);
  SIMPLE_SET_GET_METHOD(Real, OriginY);
  SIMPLE_SET_GET_METHOD(Real, OriginZ);
  SIMPLE_SET_GET_METHOD(Real, GaussianSigma);
  SIMPLE_SET_GET_METHOD(Real, ResetTime);
  SIMPLE_SET_GET_METHOD(Real, StartTime);
  SIMPLE_SET_GET_METHOD(Real, UniformLength);
  SIMPLE_SET_GET_METHOD(Real, UniformWidth);
  SIMPLE_SET_GET_METHOD(Real, UniformHeight);
  SIMPLE_SET_GET_METHOD(Real, UniformLengthX);
  SIMPLE_SET_GET_METHOD(Real, UniformLengthY);
  SIMPLE_SET_GET_METHOD(Real, UniformLengthZ);
  SIMPLE_SET_GET_METHOD(Real, UniformRadiusWidth);
  SIMPLE_SET_GET_METHOD(Real, UniformRadiusXY);
  SIMPLE_SET_GET_METHOD(Real, UniformRadiusXZ);
  SIMPLE_SET_GET_METHOD(Real, UniformRadiusYZ);
  SIMPLE_GET_METHOD(Polymorph, LengthBinFractions);
  void setLengthBinFractions(const Polymorph& aValue)
    {
      LengthBinFractions = aValue;
      PolymorphVector aValueVector(aValue.as<PolymorphVector>());
      for(unsigned i(0); i != aValueVector.size(); ++i)
        {
          theLengthBinFractions.push_back(aValueVector[i].as<double>());
        }
    }
  virtual void initialize();
  virtual void initializeSecond();
  virtual void initializeFourth();
  virtual void populateGaussian(Species*);
  virtual void populateUniformDense(Species*, unsigned int*, unsigned int*);
  virtual void populateUniformSparse(Species*);
  virtual void populateBinFractions(Species*);
  virtual void populateUniformRanged(Species*);
  virtual void populateUniformOnDiffusiveVacant(Species*);
  virtual void populateUniformOnMultiscale(Species*);
  virtual void fire();
  virtual void initializeFifth()
    {
      theInterval = ResetTime;
      theTime = StartTime+theInterval; 
      thePriorityQueue->move(theQueueID);
    }
  void checkProcess();
  virtual int getPriority()
    {
      return SpatiocyteProcess::getPriority();
    }
protected:
  double GaussianSigma;
  double OriginX;
  double OriginY;
  double OriginZ;
  double ResetTime;
  double StartTime;
  double UniformLength;
  double UniformWidth;
  double UniformHeight;
  double UniformLengthX;
  double UniformLengthY;
  double UniformLengthZ;
  double UniformRadiusWidth;
  double UniformRadiusXY;
  double UniformRadiusXZ;
  double UniformRadiusYZ;
  std::vector<double> theLengthBinFractions;
  Polymorph LengthBinFractions;
};

}

#endif /* __MoleculePopulateProcess_hpp */
