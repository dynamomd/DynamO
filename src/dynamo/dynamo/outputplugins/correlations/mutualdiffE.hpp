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
#include <dynamo/outputplugins/outputplugin.hpp>
#include <boost/circular_buffer.hpp>

namespace dynamo {
  class OPMutualDiffusionE: public OutputPlugin
  {
  public:
    OPMutualDiffusionE(const dynamo::SimData*, const magnet::xml::Node&);
  
    virtual void operator<<(const magnet::xml::Node&);

    virtual OutputPlugin* Clone() const { return new OPMutualDiffusionE(*this); }

    virtual void stream(const double);

    virtual void eventUpdate(const GlobalEvent&, const NEventData&);

    virtual void eventUpdate(const LocalEvent&, const NEventData&);

    virtual void eventUpdate(const System&, const NEventData&, const double&);
  
    virtual void eventUpdate(const IntEvent&, const PairEventData&);

    virtual double rescaleFactor();

    virtual void output(magnet::xml::XmlStream&);

    virtual void initialise();
 
    std::list<Vector  > getAvgAcc() const;
  
  protected:  
    virtual void updateDelG(const PairEventData&);

    virtual void updateDelG(const ParticleEventData&);

    virtual void updateDelG(const NEventData&);
    
    virtual void newG();
  
    virtual void accPass();

    double getdt();
    
    boost::circular_buffer<Vector  > G1;
    boost::circular_buffer<Vector  > G2;
    std::vector<Vector  > accG;
    size_t count;
    double dt, currentdt;

    Vector  delGsp1, delGsp2, Gsp1, Gsp2;

    size_t species1;
    size_t species2;
  
    Vector  sysMom;

    double massFracSp1;
    double massFracSp2;

    size_t CorrelatorLength;
    size_t currCorrLen;
    bool notReady;
    std::string species1name, species2name;

  };
}
