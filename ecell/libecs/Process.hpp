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

#ifndef ___PROCESS_H___
#define ___PROCESS_H___

#include "AssocVector.h"

#include "libecs.hpp"
#include "Entity.hpp"
#include "Connection.hpp"

namespace libecs
{

  /** @addtogroup entities
   *@{
   */

  /** @file */


  DECLARE_ASSOCVECTOR( String, Connection, std::less< const String >, 
		       ConnectionMap  );

  /**
     Process class is used to represent chemical and other phenonema which 
     may or may not result in change in value of one or more Variables.

  */

  class Process 
    : 
    public Entity
  {

  public: 


    class PriorityCompare
    {
    public:
      bool operator()( ProcessPtr aLhs, ProcessPtr aRhs ) const
      {
	return ( aLhs->getPriority() < aRhs->getPriority() );
      }
    };


    /** 
	A function type that returns a pointer to Process.  

	Every Process class must have this type of a function which returns
	an instance for the ProcessMaker.
    */

    typedef ProcessPtr (* AllocatorFuncPtr )();


  public:

    Process();
    virtual ~Process();

    virtual const EntityType getEntityType() const
    {
      return EntityType( EntityType::PROCESS );
    }

    virtual void initialize();
    
    virtual void process() = 0;
    
    
    /**
       Set activity variable.  This must be set at every turn and takes
       [number of molecule that this process yields] / [deltaT].
       However, public activity() method returns it as number of molecule
       per a second, not per deltaT.

       @param activity [number of molecule that this yields] / [deltaT].
       @see getActivity(), getActivityPerSecond()
    */

    void setActivity( RealCref anActivity ) 
    { 
      theActivity = anActivity; 
    }

    const Real getActivity() const
    {
      return theActivity;
    }

    void setConnection( PolymorphCref aValue );

    void setConnectionList( PolymorphCref );

    const Polymorph getConnectionList() const;

    void registerConnection( StringCref aName, FullIDCref aFullID,
			   const Int aCoefficient );

    void registerConnection( StringCref aName, VariablePtr aVariable, 
			   const Int aCoefficient );

    /**
       Get Connection by tag name.

       @param aConnectionName
       @return a Connection
       @see Connection
    */

    Connection getConnection( StringCref aConnectionName );

    /**
       @return a const reference to the connection map
    */
    ConnectionMapCref getConnectionMap() const
    {
      return theConnectionMap;
    }

    void setPriority( IntCref aValue )
    {
      thePriority = aValue;
    }

    const Int getPriority() const
    {
      return thePriority;
    }

  protected:

    void makeSlots();

  protected:

    ConnectionMap theConnectionMap;

  private:

    Real        theActivity;
    Int         thePriority;

  };





  /*@}*/

} // namespace libecs

#endif /* ___PROCESS_H___ */

/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/