//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-CELL Simulation Environment package
//
//                Copyright (C) 2002 Keio University
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-CELL is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// E-CELL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with E-CELL -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Kouichi Takahashi <shafi@e-cell.org> at
// E-CELL Project, Lab. for Bioinformatics, Keio University.
//

#ifndef __FIXEDSODE1_HPP
#define __FIXEDSODE1_HPP

#include <gsl/gsl_linalg.h>

#include "libecs/DifferentialStepper.hpp"

USE_LIBECS;

DECLARE_VECTOR( Int, IntVector );

LIBECS_DM_CLASS( FixedDAE1Stepper, DifferentialStepper )
{

public:

  LIBECS_DM_OBJECT( FixedDAE1Stepper, Stepper )
    {
      INHERIT_PROPERTIES( DifferentialStepper );

      PROPERTYSLOT_SET_GET( Real, PerturbationRate );
      PROPERTYSLOT_SET_GET( Real, Tolerance );
    }

  FixedDAE1Stepper( void );
  
  virtual ~FixedDAE1Stepper( void );

  SET_METHOD( Real, PerturbationRate )
  {
    thePerturbationRate = value;
  }

  GET_METHOD( Real, PerturbationRate )
  {
    return thePerturbationRate;
  }

  SET_METHOD( Real, Tolerance )
  {
    theTolerance = value;
  }

  GET_METHOD( Real, Tolerance )
  {
    return theTolerance;
  }

  virtual void initialize();

  virtual void step();

  void calculateVelocityVector();
  void calculateJacobian();

  void checkDependency();
  void resetVelocity();

  const Real solve();

  /**
     Only for debugging
  */

  void printJacobianMatrix()
  {
    for ( UnsignedInt i( 0 ); i < theSystemSize; i++ )
      for ( UnsignedInt j( 0 ); j < theSystemSize; j++ )
    	std::cout << "m(" << i << "," << j << ") = " 
    		  << gsl_matrix_get( theJacobianMatrix, i, j ) 
		  << ":" << theVariableVector[ j ]->getID() << std::endl;
  }

protected:

  UnsignedInt     theSystemSize;
  Real            thePerturbationRate;
  Real            theTolerance;

  RealVector                   theEachVelocityBuffer;
  // std::vector<ProcessVector>
  std::vector<IntVector>       theDependentProcessVector;
  std::vector<IntVector>       theDependentVariableVector;

  gsl_matrix*         theJacobianMatrix;
  gsl_vector*         theVelocityVector;
  gsl_vector*         theSolutionVector;
  gsl_permutation*    thePermutation;

  IntVector       theContinuousVariableVector;
  RealVector      theDiscreteActivityBuffer;
};

#endif /* __FIXEDSODE1_HPP */