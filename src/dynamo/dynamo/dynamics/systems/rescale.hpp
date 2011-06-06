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
#include "system.hpp"
#include "../../base/is_simdata.hpp"
#include "../ranges/1range.hpp"
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <magnet/cloneptr.hpp>

class CSysRescale: public System
{
public:
  CSysRescale(const magnet::xml::Node& XML, dynamo::SimData*);
  CSysRescale(dynamo::SimData*, size_t frequency, std::string name);

  virtual System* Clone() const { return new CSysRescale(*this); }

  virtual void runEvent() const;

  virtual void initialise(size_t);

  virtual void operator<<(const magnet::xml::Node&);

  void checker(const NEventData&);
  
  inline const long double& getScaleFactor() const {return scaleFactor; }

protected:
  virtual void outputXML(xml::XmlStream&) const;

  size_t _frequency;

  mutable long double scaleFactor;

  mutable long double LastTime, RealTime;
  
};
