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

#include "System.hpp"
#include "FullID.hpp"
#include "PropertySlotMaker.hpp"

#include "Entity.hpp"


namespace libecs
{

  Entity::Entity()
    : 
    theSuperSystem( NULLPTR ),
    theID( "" ),
    theName( "" ) 
  {
    CREATE_PROPERTYSLOT_SET_GET( String, Name,   Entity );
    CREATE_PROPERTYSLOT        ( String, FullID, 
				 NULLPTR, &Entity::getFullIDString );
  }


  Entity::~Entity()
  {
    ; // do nothing
  }

  const FullID Entity::getFullID() const
  {
    return FullID( getEntityType(), getSystemPath(), getID() );
  }

  const String Entity::getFullIDString() const
  {
    return getFullID().getString();
  }

  const SystemPath Entity::getSystemPath() const
  {
    SystemPtr aSystemPtr( getSuperSystem() );
    SystemPath aSystemPath( aSystemPtr->getSystemPath() );
    aSystemPath.push_back( aSystemPtr->getID() );
    return aSystemPath;
  }

} // namespace libecs

/*
  Do not modify
  $Author$
  $Revision$
  $Date$
  $Locker$
*/
