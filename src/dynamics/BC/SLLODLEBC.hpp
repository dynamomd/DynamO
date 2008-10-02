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


#ifndef SLLOD_LE_BC_H
#define SLLOD_LE_BC_H

#include "LEBC.hpp"

/*
template<class T>
class SLLOD_LE_BC: public LEBC<T>
{
 public:
  inline SLLOD_LE_BC():LEBC<T>() 
    { 
      LEBC<T>::I_cout() << "SLLOD variant Lee's Edwards BC loaded"; 
    };

  virtual ~SLLOD_LE_BC() {};
  
  inline virtual BC* Clone () const
    {
      return new SLLOD_LE_BC(*this);
    };

  inline virtual void setPBC(CVector<> &pos, CVector<> &vel) const 
    {
      LEBC<T>::setPBC(pos);
    };
};

*/
#endif
