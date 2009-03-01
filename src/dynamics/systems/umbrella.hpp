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

#ifndef CSUmbrella_HPP
#define CSUmbrella_HPP

#include "system.hpp"
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include "../../extcode/include/boost/random/normal_distribution.hpp"
#include "../../base/is_simdata.hpp"
#include "../ranges/1range.hpp"
#include "../../datatypes/pluginpointer.hpp"

class CSUmbrella: public CSystem
{
public:
  CSUmbrella(const XMLNode& XML, DYNAMO::SimData*);

  CSUmbrella(DYNAMO::SimData*, Iflt, Iflt, Iflt, std::string, CRange*, CRange*);
  
  virtual CSystem* Clone() const { return new CSUmbrella(*this); }

  virtual void runEvent() const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  void particlesUpdated(const CNParticleData&);

  void recalculateTime();

  Iflt a,b,delU;
  int ulevelcenter;
  mutable int ulevel;

  smrtPlugPtr<CRange> range1;
  smrtPlugPtr<CRange> range2;
};

#endif
