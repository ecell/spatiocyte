//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-CELL Simulation Environment package
//
//                Copyright (C) 1996-2002 Keio University
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

#include "FullID.hpp"

#include "DiscreteEventStepper.hpp"


namespace libecs
{

  LIBECS_DM_INIT_STATIC( DiscreteEventStepper, Stepper );


  //////////////////// DiscreteEventStepper

  DiscreteEventStepper::DiscreteEventStepper()
    :
    theTimeScale( 0.0 ),
    theTolerance( 0.0 )
  {
    ; // do nothing
  }


  GET_METHOD_DEF( String, LastProcess, DiscreteEventStepper )
  {
    return theLastProcess->getFullID().getString();
  }

  void DiscreteEventStepper::initialize()
  {
    Stepper::initialize();

    // dynamic_cast theProcessVector of this Stepper to DiscreteEventProcess,
    // and store it in theDiscreteEventProcessVector.
    theDiscreteEventProcessVector.clear();
    try
      {
	std::transform( theProcessVector.begin(), theProcessVector.end(),
			std::back_inserter( theDiscreteEventProcessVector ),
			DynamicCaster<DiscreteEventProcessPtr,ProcessPtr>() );
      }
    catch( const libecs::TypeError& )
      {
	THROW_EXCEPTION( InitializationFailed,
			 String( getClassName() ) + 
			 ": Only DiscreteEventProcesses are allowed to exist "
			 "in this Stepper." );
      }

    // optimization: sort by memory address
    std::sort( theDiscreteEventProcessVector.begin(), 
	       theDiscreteEventProcessVector.end() );


    // (1) check Process dependency
    // (2) update step interval of each Process
    // (3) construct the priority queue (scheduler)
    thePriorityQueue.clear();
    const Real aCurrentTime( getCurrentTime() );
    for( DiscreteEventProcessVector::const_iterator i( theDiscreteEventProcessVector.begin() );
	 i != theDiscreteEventProcessVector.end(); ++i )
      {      
	DiscreteEventProcessPtr anDiscreteEventProcessPtr( *i );
	
	// check Process dependencies
	anDiscreteEventProcessPtr->clearDependentProcessVector();
	// here assume aCoefficient != 0
	for( DiscreteEventProcessVector::const_iterator 
	       j( theDiscreteEventProcessVector.begin() );
	     j != theDiscreteEventProcessVector.end(); ++j )
	  {
	    DiscreteEventProcessPtr const anDiscreteEventProcess2Ptr( *j );
	  
	    if( anDiscreteEventProcessPtr->
		checkProcessDependency( anDiscreteEventProcess2Ptr ) )
	      {
		anDiscreteEventProcessPtr->
		  addDependentProcess( anDiscreteEventProcess2Ptr );
	      }
	  }

	// warning: implementation dependent
	// here we assume size() is the index of the newly pushed element
	const Int anIndex( thePriorityQueue.size() );

	anDiscreteEventProcessPtr->setIndex( anIndex );
	anDiscreteEventProcessPtr->updateStepInterval();
	thePriorityQueue.push( StepperEvent( anDiscreteEventProcessPtr->getStepInterval()
					+ aCurrentTime,
					anDiscreteEventProcessPtr ) );
      }

    // here all the DiscreteEventProcesses are updated, then set new
    // step interval and reschedule this stepper.
    // That means, this Stepper doesn't necessary steps immediately
    // after initialize()
    StepperEventCref aTopEvent( thePriorityQueue.top() );
    const Real aNewTime( aTopEvent.getTime() );

    setStepInterval( aNewTime - aCurrentTime );
    getModel()->reschedule( this );

  }
  

  // this doesn't necessary occur at the first step of the simulation,
  // and imediately after initialize(), because initialize() recalculates
  // all propensities and reschedules this stepper.
  void DiscreteEventStepper::step()
  {
    StepperEventCref anEvent( thePriorityQueue.top() );

    DiscreteEventProcessPtr const aMuProcess( anEvent.getProcess() );
    aMuProcess->process();
    theLastProcess = aMuProcess;

    const Real aCurrentTime( getCurrentTime() );

    // Update relevant processes
    DiscreteEventProcessVectorCref 
      aDependentProcessVector( aMuProcess->getDependentProcessVector() );
    for ( DiscreteEventProcessVectorConstIterator 
	    i( aDependentProcessVector.begin() );
	  i != aDependentProcessVector.end(); ++i ) 
      {
	DiscreteEventProcessPtr const anAffectedProcess( *i );
	anAffectedProcess->updateStepInterval();
	const Real aStepInterval( anAffectedProcess->getStepInterval() );
	// aTime is time in the priority queue

	Int anIndex( anAffectedProcess->getIndex() );
	thePriorityQueue.changeOneKey( anIndex,
				       StepperEvent( aStepInterval + aCurrentTime,
						anAffectedProcess ) );
      }

    StepperEventCref aTopEvent( thePriorityQueue.top() );
    const Real aNextStepInterval( aTopEvent.getTime() - aCurrentTime );

    DiscreteEventProcessPtr const aNewTopProcess( aTopEvent.getProcess() );

    // Calculate new timescale.
    // Some benchmark showed this doesn't affect performance.
    theTimeScale = theTolerance * aNewTopProcess->getTimeScale();

    setStepInterval( aNextStepInterval );
  }


  void DiscreteEventStepper::interrupt( StepperPtr const aCaller )
  {
    // update step intervals of DiscreteEventProcesses
    const Real aCurrentTime( aCaller->getCurrentTime() );
    for( DiscreteEventProcessVector::const_iterator i( theDiscreteEventProcessVector.begin() );
	 i != theDiscreteEventProcessVector.end(); ++i )
      {      
	DiscreteEventProcessPtr const anDiscreteEventProcessPtr( *i );
	
	anDiscreteEventProcessPtr->updateStepInterval();
	const Real aStepInterval( anDiscreteEventProcessPtr->getStepInterval() );

	thePriorityQueue.changeOneKey( anDiscreteEventProcessPtr->getIndex(),
				       StepperEvent( aStepInterval + aCurrentTime,
						     anDiscreteEventProcessPtr ) );
      }

    StepperEventCref aTopEvent( thePriorityQueue.top() );
    const Real aNewTime( aTopEvent.getTime() );

    // reschedule this Stepper to aNewStepInterval past current time of
    // the interruption caller.
    setStepInterval( aNewTime - getCurrentTime() );
    getModel()->reschedule( this );
  }


} // namespace libecs


/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/

