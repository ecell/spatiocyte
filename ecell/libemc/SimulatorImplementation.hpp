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


#ifndef __SIMULATORIMPLEMENTATION_HPP
#define __SIMULATORIMPLEMENTATION_HPP

#include "libecs/libecs.hpp"

#include "libemc.hpp"

namespace libemc
{

  /** @addtogroup libemc_module 
   * @{ 
   */ 

  /** @file */

  class EventHandler
    :
    public std::unary_function<void,void> 
  {
  public:
    EventHandler() {}
    virtual ~EventHandler() {}

    virtual void operator()( void ) const = 0;
  };

  class PendingEventChecker
    :
    public std::unary_function<bool,void>
  {
  public:
    PendingEventChecker() {}
    virtual ~PendingEventChecker() {}

    virtual bool operator()( void ) const = 0;
  };

  class DefaultPendingEventChecker
    :
    public PendingEventChecker
  {
  public:
    DefaultPendingEventChecker() {}
    //    virtual ~DefaultPendingEventChecker() {}

    virtual bool operator()( void ) const
    {
      return false;
    }
  };



  /**
     Pure virtual base class (interface definition) of simulator
     implementation.
  */

  class SimulatorImplementation
  {

  public:

    SimulatorImplementation() {}
    virtual ~SimulatorImplementation() {}


    virtual void createStepper( libecs::StringCref         aClassname,
				libecs::StringCref         anId ) = 0;

    virtual const libecs::Polymorph getStepperList() const = 0;

    virtual void setStepperProperty( libecs::StringCref    aStepperID,
				     libecs::StringCref    aPropertyName,
				     libecs::PolymorphCref aValue ) = 0;

    virtual const libecs::Polymorph
    getStepperProperty( libecs::StringCref aStepperID,
			libecs::StringCref aPropertyName ) const = 0;

    virtual void createEntity( libecs::StringCref   aClassname, 
			       libecs::StringCref   aFullIDString,
			       libecs::StringCref   aName ) = 0;

    virtual const bool 
    isEntityExist( libecs::StringCref  aFullIDString ) const = 0;

    virtual void setEntityProperty( libecs::StringCref    aFullPNString,
				    libecs::PolymorphCref aValue ) = 0;

    virtual const libecs::Polymorph
    getEntityProperty( libecs::StringCref aFullPNString ) const = 0;


    virtual void createLogger( libecs::StringCref aFullPNString ) = 0;

    virtual const libecs::Polymorph getLoggerList() const = 0;

    virtual const libecs::DataPointVectorRCPtr 
    getLoggerData( libecs::StringCref aFullPNString ) const = 0;

    virtual const libecs::DataPointVectorRCPtr
    getLoggerData( libecs::StringCref aFullPNString, 
		   libecs::RealCref aStartTime, 
		   libecs::RealCref anEndTime ) const = 0;

    virtual const libecs::DataPointVectorRCPtr
    getLoggerData( libecs::StringCref aFullPNString,
		   libecs::RealCref aStartTime, libecs::RealCref anEndTime, 
		   libecs::RealCref interval ) const = 0;

    virtual const libecs::Real 
    getLoggerStartTime( libecs::StringCref aFullPNString ) const = 0;

    virtual const libecs::Real 
    getLoggerEndTime( libecs::StringCref aFullPNString ) const = 0;

    virtual void 
    setLoggerMinimumInterval( libecs::StringCref aFullPNString, 
			      libecs::RealCref anInterval ) = 0;

    virtual const libecs::Real 
    getLoggerMinimumInterval( libecs::StringCref aFullPNString ) const = 0;

    virtual const libecs::Int 
    getLoggerSize( libecs::StringCref aFullPNString ) const = 0;

    virtual void step() = 0;

    virtual void initialize() = 0;

    virtual const libecs::Real getCurrentTime() const = 0;

    virtual void run() = 0;

    virtual void run( libecs::Real aDuration ) = 0;

    virtual void stop() = 0;

    virtual void setPendingEventChecker( PendingEventCheckerPtr 
					 aPendingEventChecker ) = 0;

    virtual void setEventHandler( EventHandlerPtr anEventHandler ) = 0;

  };   //end of class Simulator

  /** @} */ //end of libemc_module 

} // namespace libemc

#endif   /* ___SIMULATOR_IMPLEMENTATION_H___ */

