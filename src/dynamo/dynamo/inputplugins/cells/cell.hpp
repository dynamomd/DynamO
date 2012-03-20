/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#pragma once
#include <vector>
#include <boost/scoped_ptr.hpp>

namespace dynamo {
  class UCell
  {
  public:
    UCell(UCell* nP): uc(nP) {}

    virtual ~UCell() {}

    virtual void initialise() { uc->initialise(); }

    virtual std::vector<Vector  > placeObjects(const Vector & ) = 0;  
  
    boost::scoped_ptr<UCell> uc;
  };

  /*! \brief A simple terminator, used to place a particle at this
      point*/
  class UParticle: public UCell
  {
  public:
    UParticle(): UCell(NULL) {}

    //Terminate initialisation
    virtual void initialise() {}

    virtual std::vector<Vector  > placeObjects(const Vector & center)
    {
      std::vector<Vector  > retval;
      retval.push_back(center);
      return retval;
    }
  };

  /*! \brief A simple terminator, used to place a particle at this
      point*/
  class UList: public UCell
  {
  public:
    UList(const std::vector<Vector>& list, double scale, UCell* nextCell): 
      UCell(nextCell), 
      _list(list) 
    {
      //Center the list of positions
      
      Vector center(0,0,0);
      BOOST_FOREACH(const Vector& vec, _list)
	center += vec;
      
      center /= _list.size();

      BOOST_FOREACH(Vector& vec, _list)
	vec = scale * (vec - center);
    }

    virtual std::vector<Vector  > placeObjects(const Vector & center)
    {
      std::vector<Vector> retval;

      BOOST_FOREACH(const Vector& vec, _list)
	{
	  const std::vector<Vector>& newsites = uc->placeObjects(vec + center);
	  retval.insert(retval.end(), newsites.begin(), newsites.end());
	}

      return retval;
    }

    std::vector<Vector> _list;
  };
}
