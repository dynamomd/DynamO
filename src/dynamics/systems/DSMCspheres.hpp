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

#ifndef CSDSMCSpheres_HPP
#define CSDSMCSpheres_HPP

#include "system.hpp"
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include "../../extcode/include/boost/random/normal_distribution.hpp"
#include "../../base/is_simdata.hpp"
#include "../ranges/1range.hpp"
#include "../../datatypes/pluginpointer.hpp"

class CSDSMCSpheres: public CSystem
{
public:
  CSDSMCSpheres(const XMLNode& XML, DYNAMO::SimData*);

  CSDSMCSpheres(DYNAMO::SimData*, Iflt, Iflt, Iflt, Iflt, std::string, CRange*, CRange*);
  
  virtual CSystem* Clone() const { return new CSDSMCSpheres(*this); }

  virtual void stream(Iflt);

  virtual void runEvent() const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  mutable boost::variate_generator<DYNAMO::baseRNG&, 
				   boost::uniform_real<> > uniformRand;  

  Iflt tstep;
  Iflt chi;
  Iflt d2;
  Iflt diameter;
  mutable Iflt maxprob;
  Iflt e;
  Iflt factor;

  smrtPlugPtr<CRange> range1;
  smrtPlugPtr<CRange> range2;
};

#endif
