//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 1996-2002 Keio University
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
// written by Kouichi Takahashi <shafi@e-cell.org>,
// E-Cell Project, Institute for Advanced Biosciences, Keio University.
//

#ifndef __DIFFERENTIALSTEPPER_HPP
#define __DIFFERENTIALSTEPPER_HPP

#include "libecs.hpp"
#include "Stepper.hpp"


namespace libecs
{

  /** @addtogroup stepper
   *@{
   */

  /** @file */

  /**
     DIFFERENTIAL EQUATION SOLVER


  */

  DECLARE_CLASS( DifferentialStepper );

  LIBECS_DM_CLASS( DifferentialStepper, Stepper )
  {

  public:

    LIBECS_DM_OBJECT_ABSTRACT( DifferentialStepper )
      {
	INHERIT_PROPERTIES( Stepper );

	// FIXME: load/save ??
	PROPERTYSLOT( Real, StepInterval,
		      &DifferentialStepper::initializeStepInterval,
		      &DifferentialStepper::getStepInterval );
	
	PROPERTYSLOT_GET_NO_LOAD_SAVE( Real, NextStepInterval );
      }

    class Interpolant
      :
      public libecs::Interpolant
    {
    public:
      Interpolant( DifferentialStepperRef aStepper, 
		     VariablePtr const aVariablePtr )
	:
	libecs::Interpolant( aVariablePtr ),
	theStepper( aStepper ),
	theIndex( theStepper.getVariableIndex( aVariablePtr ) )
      {
	; // do nothing
      }

      virtual const Real getDifference( RealParam aTime, 
					RealParam anInterval )
      {
	// First order interpolation.  This should be overridden in
	// higher order DifferentialSteppers.
	return theStepper.getVelocityBuffer()[ theIndex ] * anInterval;
      }
      
    protected:

      DifferentialStepperRef    theStepper;
      VariableVector::size_type theIndex;

    };


  public:

    DifferentialStepper();
    virtual ~DifferentialStepper();

    /**
       Override setStepInterval() for theTolerableStepInterval.

    virtual SET_METHOD( Real, StepInterval )
    {
      Stepper::setStepInterval( value );
    }

    */

    SET_METHOD( Real, NextStepInterval )
    {
      theNextStepInterval = value;
    }

    GET_METHOD( Real, NextStepInterval )
    {
      return theNextStepInterval;
    }

    void initializeStepInterval( RealParam aStepInterval )
    {
      setStepInterval( aStepInterval );
      setOriginalStepInterval( aStepInterval );
      setTolerableStepInterval( aStepInterval );
      setNextStepInterval( aStepInterval );
    }

    void resetAll();
    void interIntegrate();
 
    virtual void initialize();

    virtual void reset();

    virtual void interrupt( StepperPtr const aCaller );


    RealVectorCref getVelocityBuffer() const
    {
      return theVelocityBuffer;
    }

    virtual InterpolantPtr createInterpolant( VariablePtr aVariable )
    {
      return new DifferentialStepper::Interpolant( *this, aVariable );
    }

  protected:

    const bool isExternalErrorTolerable() const;

  protected:

    RealVector theVelocityBuffer;
    bool                   theStateFlag;

  private:

    Real theNextStepInterval;

  };


  /**
     ADAPTIVE STEPSIZE DIFFERENTIAL EQUATION SOLVER


  */

  DECLARE_CLASS( AdaptiveDifferentialStepper );

  LIBECS_DM_CLASS( AdaptiveDifferentialStepper, DifferentialStepper )
  {

  public:

    LIBECS_DM_OBJECT_ABSTRACT( AdaptiveDifferentialStepper )
      {
	INHERIT_PROPERTIES( DifferentialStepper );

	PROPERTYSLOT_SET_GET( Real, Tolerance );
	PROPERTYSLOT_SET_GET( Real, AbsoluteToleranceFactor );
	PROPERTYSLOT_SET_GET( Real, StateToleranceFactor );
	PROPERTYSLOT_SET_GET( Real, DerivativeToleranceFactor );

	PROPERTYSLOT( Integer, IsEpsilonChecked,
		      &AdaptiveDifferentialStepper::setEpsilonChecked,
		      &AdaptiveDifferentialStepper::isEpsilonChecked );
	PROPERTYSLOT_SET_GET( Real, AbsoluteEpsilon );
	PROPERTYSLOT_SET_GET( Real, RelativeEpsilon );

	PROPERTYSLOT_GET_NO_LOAD_SAVE( Real, MaxErrorRatio );
	PROPERTYSLOT_GET_NO_LOAD_SAVE( Integer,  Order );
      }


    class Interpolant
      :
      public libecs::Interpolant
    {
    public:

      Interpolant( AdaptiveDifferentialStepperRef aStepper, 
		     VariablePtr const aVariablePtr )
	:
	libecs::Interpolant( aVariablePtr ),
	theStepper( aStepper ),
	theIndex( theStepper.getVariableIndex( aVariablePtr ) )
      {
	; // do nothing
      }

      virtual const Real getDifference( RealParam aTime, 
					RealParam anInterval )
      {
	const Real aTolerableStepInterval
	  ( theStepper.getTolerableStepInterval() );

	const Real aTimeInterval( aTime - theStepper.getCurrentTime() );

	const Real theta( FMA( aTimeInterval, 2.0, - anInterval )
			  / aTolerableStepInterval );

	const Real k1( theStepper.getK1()[ theIndex ] );
	const Real k2( theStepper.getVelocityBuffer()[ theIndex ] );

	return ( FMA( theta, ( k2 - k1 ), k1 ) * anInterval );
      }

    protected:

      AdaptiveDifferentialStepperRef theStepper;
      VariableVector::size_type      theIndex;
    };

  public:

    AdaptiveDifferentialStepper();
    virtual ~AdaptiveDifferentialStepper();

    /**
       Adaptive stepsize control.

       These methods are for handling the standerd error control objects.
    */

    SET_METHOD( Real, Tolerance )
    {
      theTolerance = value;
    }

    GET_METHOD( Real, Tolerance )
    {
      return theTolerance;
    }

    SET_METHOD( Real, AbsoluteToleranceFactor )
    {
      theAbsoluteToleranceFactor = value;
    }

    GET_METHOD( Real, AbsoluteToleranceFactor )
    {
      return theAbsoluteToleranceFactor;
    }

    SET_METHOD( Real, StateToleranceFactor )
    {
      theStateToleranceFactor = value;
    }

    GET_METHOD( Real, StateToleranceFactor )
    {
      return theStateToleranceFactor;
    }

    SET_METHOD( Real, DerivativeToleranceFactor )
    {
      theDerivativeToleranceFactor = value;
    }

    GET_METHOD( Real, DerivativeToleranceFactor )
    {
      return theDerivativeToleranceFactor;
    }

    SET_METHOD( Real, MaxErrorRatio )
    {
      theMaxErrorRatio = value;
    }

    GET_METHOD( Real, MaxErrorRatio )
    {
      return theMaxErrorRatio;
    }

    /**
       check difference in one step
    */

    SET_METHOD( Integer, EpsilonChecked )
    {
      if ( value > 0 ) {
	theEpsilonChecked = true;
      }
      else {
	theEpsilonChecked = false;
      }
    }

    const Integer isEpsilonChecked() const
    {
      return theEpsilonChecked;
    }

    SET_METHOD( Real, AbsoluteEpsilon )
    {
      theAbsoluteEpsilon = value;
    }

    GET_METHOD( Real, AbsoluteEpsilon )
    {
      return theAbsoluteEpsilon;
    }

    SET_METHOD( Real, RelativeEpsilon )
    {
      theRelativeEpsilon = value;
    }

    GET_METHOD( Real, RelativeEpsilon )
    {
      return theRelativeEpsilon;
    }

    virtual GET_METHOD( Integer, Order )
    { 
      return 1; 
    }

    virtual void initialize();
    virtual void step();
    virtual bool calculate() = 0;

    RealVectorCref getK1() const
    {
      return theK1;
    }

    virtual InterpolantPtr createInterpolant( VariablePtr aVariable )
    {
      return new AdaptiveDifferentialStepper::Interpolant( *this, aVariable );
    }

  protected:

    Real safety;
    RealVector theK1;

  private:

    Real theTolerance;
    Real theAbsoluteToleranceFactor;
    Real theStateToleranceFactor;
    Real theDerivativeToleranceFactor;

    bool    theEpsilonChecked;
    Real    theAbsoluteEpsilon;
    Real    theRelativeEpsilon;

    Real theMaxErrorRatio;
  };


} // namespace libecs

#endif /* __DIFFERENTIALSTEPPER_HPP */



/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/
