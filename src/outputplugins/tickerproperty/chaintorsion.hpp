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

#ifndef COPCTorsion_H
#define COPCTorsion_H

#include "ticker.hpp"
#include "../../datatypes/histogram.hpp"

class CTChain;

class COPCTorsion: public COPTicker
{
 public:
  COPCTorsion(const DYNAMO::SimData*, const XMLNode&);

  virtual COutputPlugin *Clone() const
  { return new COPCTorsion(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();

  virtual void changeSystem(COutputPlugin*);

  virtual void output(xmlw::XmlStream&);
  
 protected:

  struct CTCdata
  {
    const CTChain* chainPtr;
    C1DHistogram gammaMol;
    C1DHistogram gammaSys;
    C1DHistogram f;
    CTCdata(const CTChain* ptr, Iflt binwidth1, Iflt binwidth2, Iflt binwidth3):
      chainPtr(ptr), gammaMol(binwidth1),
      gammaSys(binwidth2),f(binwidth3)
    {}

  };

  std::list<CTCdata> chains;
};

#endif
