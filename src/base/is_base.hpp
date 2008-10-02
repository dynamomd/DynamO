/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef IS_Base_H
#define IS_Base_H

#include "is_stream_op.hpp"
#include "constants.hpp"

namespace DYNAMO
{
  struct SimData;

  class Base_Class
  {
  public:
    Base_Class(const char *aName, const char *aColor):
      name(aName),color(aColor) {};
    
    Stream_Operator I_cout() const
    {
      return ((std::cout << Stream_Operator(name,color)) << "\n");
    }

    Stream_Operator I_cerr() const
    {
      return ((std::cout << Stream_Operator(name,IC_red)) << "\n");
    }
    
  protected:
    const char* name;
    const char* color;
  };

  class SimBase: public Base_Class
  {
  public:
    SimBase(SimData* const& SD,const char *aName, const char *aColor):
      Base_Class(aName,aColor),
      Sim(SD)    
    {};
    
  protected:
    SimData* Sim;
  };
  
  class SimBase_const: public Base_Class
  {
  public:
    SimBase_const(const SimData* const& SD, const char *aName, const char *aColor):
      Base_Class(aName,aColor),
      Sim(SD)
    {};

  protected:
    const SimData* Sim;
  };

}
#endif
