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

#ifndef __RULESRMPROCESS_HPP
#define __RULESRMPROCESS_HPP

#include "SRMProcess.hpp"

namespace libecs
{

  DECLARE_CLASS( RuleSRMProcess );

  class RuleSRMProcess
    :
    public SRMProcess
  {

  public:

    class IsRuleProcess
      : 
      public std::unary_function<ProcessPtr,bool>
    {
    public:
      result_type operator()( const argument_type aProcessPtr ) const
      {
	if( dynamic_cast<RuleSRMProcessPtr> ( aProcessPtr ) != NULLPTR ) 
	  {
	    return true;
	  }
	else
	  {
	    return false;
	  }
      }
    };



    RuleSRMProcess()
    {
      ; // do nothing
    }

    virtual ~RuleSRMProcess()
    {
      ; // do nothing
    }


  };

} // namespace libecs


#endif /* __RULESRMPROCESS_HPP */