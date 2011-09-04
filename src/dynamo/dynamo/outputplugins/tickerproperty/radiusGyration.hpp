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

#ifdef DYNAMO_GSL

#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <dynamo/datatypes/histogram.hpp>

class CRange;

namespace dynamo {
  class CTChain;

  class OPRGyration: public OPTicker
  {
  public:
    OPRGyration(const dynamo::SimData*, const magnet::xml::Node&);

    virtual OutputPlugin *Clone() const
    { return new OPRGyration(*this); }

    virtual void initialise();

    virtual void stream(double) {}

    virtual void ticker();

    virtual void changeSystem(OutputPlugin*);

    virtual void output(magnet::xml::XmlStream&);

    struct molGyrationDat
    {
      Vector  EigenVal;
      Vector  EigenVec[3];
      Vector  MassCentre;
    };
  
    static molGyrationDat getGyrationEigenSystem(const std::tr1::shared_ptr<CRange>&, const dynamo::SimData*);

    static Vector  NematicOrderParameter(const std::list<Vector  >&);
    static double CubaticOrderParameter(const std::list<Vector  >&);

    virtual void operator<<(const magnet::xml::Node&);
  
  protected:

    struct CTCdata
    {
      const CTChain* chainPtr;
      std::vector<C1DHistogram> gyrationRadii;
      std::vector<C1DHistogram> nematicOrder;
      C1DHistogram cubaticOrder;    

      CTCdata(const CTChain* ptr, double binwidth1, double binwidth2, double binwidth3):
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

    double binwidth1, binwidth2, binwidth3;  
  };
}

#endif
