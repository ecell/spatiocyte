//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-CELL Simulation Environment package
//
//                Copyright (C) 1996-2000 Keio University
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
// written by Masayuki Okayama <smash@e-cell.org> at
// E-CELL Project, Lab. for Bioinformatics, Keio University.
//


#ifndef __EMC_DATAPOINT_HPP
#define __EMC_DATAPOINT_HPP



#include "libecs/libecs.hpp"
#include "libecs/Defs.hpp"
#include "libecs/CoreLinuxCompatibility.hpp"
#include "libecs/DataPoint.hpp"

namespace libemc
{

  /** @defgroup libemc_module The Libemc Module 
   * This is the libemc module 
   * @{ 
   */ 
  
  using namespace libecs;

  class EmcDataPoint
  {
  public:

    EmcDataPoint( void )
    {
      ; // do nothing
    }

    EmcDataPoint( const DataPoint& dp )
      :
      theDataPoint( &dp )
    {
      ; // do nothing
    }

    virtual ~EmcDataPoint()
    {
      ;
    }

    void setDataPoint( const DataPoint& dp )
    {
      theDataPoint = &dp;
    }

    const double& getTime() const
    {
      return theDataPoint->getTime();
    }

    const UVariable& getValue() const
    {
      return theDataPoint->getValue();
    }

    const DataPoint& getDataPoint() const
    {
      return *theDataPoint;
    }


  private:
    const DataPoint* theDataPoint;
  };

  /** @} */ //end of libemc_module 
}


#endif
