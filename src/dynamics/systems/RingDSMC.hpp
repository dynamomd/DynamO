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

#ifndef CSRingDSMC_HPP
#define CSRingDSMC_HPP

#include "system.hpp"
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include "../../extcode/include/boost/random/normal_distribution.hpp"
#include "../../base/is_simdata.hpp"
#include "../ranges/1range.hpp"
#include "../../datatypes/pluginpointer.hpp"

class CSRingDSMC: public CSystem
{
public:
  CSRingDSMC(const XMLNode& XML, DYNAMO::SimData*);

  CSRingDSMC(DYNAMO::SimData*, Iflt, Iflt, Iflt, Iflt, Iflt, std::string, CRange*);
  
  virtual CSystem* Clone() const { return new CSRingDSMC(*this); }

  virtual void runEvent() const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  mutable boost::variate_generator<DYNAMO::baseRNG&, 
				   boost::uniform_real<> > uniformRand;  

  Iflt tstep;
  Iflt chi12, chi13;
  Iflt d2;
  Iflt diameter;
  mutable Iflt maxprob12;
  mutable Iflt maxprob13;
  Iflt e;
  Iflt factor12;
  Iflt factor13;

  smrtPlugPtr<CRange> range1;
};

#endif
