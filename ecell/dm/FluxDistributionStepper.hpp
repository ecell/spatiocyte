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

#ifndef __FLUXDISTRIBUTIONSTEPPER_HPP
#define __FLUXDISTRIBUTIONSTEPPER_HPP

#include "libecs/libecs.hpp"

#include "libecs/VariableProxy.hpp"
#include "libecs/Stepper.hpp"
#include "libecs/DifferentialStepper.hpp"

USE_LIBECS;

LIBECS_DM_CLASS( FluxDistributionStepper, DifferentialStepper )
{  

 public:
  
  LIBECS_DM_OBJECT( FluxDistributionStepper, DifferentialStepper )
    {
      INHERIT_PROPERTIES( Stepper );
    }
  
  FluxDistributionStepper();
  ~FluxDistributionStepper() {}
  
  virtual void initialize();

  virtual void interrupt( StepperPtr const aCaller )
    {
      integrate( aCaller->getCurrentTime() );
      step();
      log();
    }
  
  virtual void step()
    {
      clear();
      process();

      for( UnsignedInt c( 0 ); c < getReadOnlyVariableOffset(); ++c )
	{
	  theVelocityBuffer[ c ] = theVariableVector[ c ]->getVelocity();
	}
    }
  
};

#endif /* __FLUXDISTRIBUTIONSTEPPER_HPP */
