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

#ifndef COPRGyration_H
#define COPRGyration_H

#include "ticker.hpp"
#include "../../datatypes/histogram.hpp"
class CRange;

class CTChain;

class COPRGyration: public COPTicker
{
 public:
  COPRGyration(const DYNAMO::SimData*, const XMLNode&);

  virtual COutputPlugin *Clone() const
  { return new COPRGyration(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();

  virtual void changeSystem(COutputPlugin*);

  virtual void output(xmlw::XmlStream&);

  struct molGyrationDat
  {
    CVector<> EigenVal;
    CVector<> EigenVec[3];
    CVector<> MassCentre;
  };
  
  static molGyrationDat getGyrationEigenSystem(const smrtPlugPtr<CRange>&, const DYNAMO::SimData*);

  static CVector<> NematicOrderParameter(const std::list<CVector<> >&);
  static Iflt CubaticOrderParameter(const std::list<CVector<> >&);
  
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
      for (int i = 0; i < NDIM; i++)
	{
	  gyrationRadii.push_back(C1DHistogram(binwidth1));
	  nematicOrder.push_back(C1DHistogram(binwidth2));
	}
    }

  };

  std::list<CTCdata> chains;
  
};

#endif
