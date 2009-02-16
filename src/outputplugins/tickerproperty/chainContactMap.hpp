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

#ifndef COPCContactMap_H
#define COPCContactMap_H

#include "ticker.hpp"
#include "../../datatypes/histogram.hpp"
#include <boost/shared_array.hpp>

class CTChain;
class CRange;

class COPCContactMap: public COPTicker
{
 public:
  COPCContactMap(const DYNAMO::SimData*, const XMLNode&);

  virtual COutputPlugin *Clone() const
  { return new COPCContactMap(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();

  virtual void changeSystem(COutputPlugin*);

  void temperatureRescale(const Iflt&) {}

  virtual void output(xmlw::XmlStream&);
  
 protected:

  struct Cdata
  {
    Cdata(const CTChain*, unsigned long);

    const CTChain* chainPtr;
    boost::shared_array<unsigned long> array;
    unsigned long counter;
    unsigned long chainlength;
  };

  std::list<Cdata> chains;

};

#endif
