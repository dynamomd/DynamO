/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef OPRGyration_H
#define OPRGyration_H

#include "ticker.hpp"
#include "../../datatypes/histogram.hpp"
class CRange;

class CTChain;

class OPRGyration: public OPTicker
{
 public:
  OPRGyration(const DYNAMO::SimData*, const XMLNode&);

  virtual OutputPlugin *Clone() const
  { return new OPRGyration(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();

  virtual void changeSystem(OutputPlugin*);

  virtual void output(xml::XmlStream&);

  struct molGyrationDat
  {
    Vector  EigenVal;
    Vector  EigenVec[3];
    Vector  MassCentre;
  };
  
  static molGyrationDat getGyrationEigenSystem(const ClonePtr<CRange>&, const DYNAMO::SimData*);

  static Vector  NematicOrderParameter(const std::list<Vector  >&);
  static Iflt CubaticOrderParameter(const std::list<Vector  >&);

  virtual void operator<<(const XMLNode&);
  
 protected:

  struct CTCdata
  {
    const CTChain* chainPtr;
    std::vector<C1DHistogram> gyrationRadii;
    std::vector<C1DHistogram> nematicOrder;
    C1DHistogram cubaticOrder;    

    CTCdata(const CTChain* ptr, Iflt binwidth1, Iflt binwidth2, Iflt binwidth3):
      chainPtr(ptr),
      cubaticOrder(binwidth3)
    {
      for (size_t i = 0; i < NDIM; i++)
	{
	  gyrationRadii.push_back(C1DHistogram(binwidth1));
	  nematicOrder.push_back(C1DHistogram(binwidth2));
	}
    }

  };

  std::list<CTCdata> chains;

  Iflt binwidth1, binwidth2, binwidth3;  
};

#endif
