/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <vector>

namespace dynamo {
  class Particle;
  class OPSHCrystal: public OPTicker
  {
  public:
    OPSHCrystal(const dynamo::Simulation*, const magnet::xml::Node&);

    virtual void initialise();

    virtual void stream(double) {}

    virtual void ticker();
  
    virtual void output(magnet::xml::XmlStream&);

    virtual void operator<<(const magnet::xml::Node&);

  protected:

    std::complex<double> localq(const Particle& part, int l, int m);

    //! Cut-off radius 
    double rg;
    size_t maxl;
    size_t nblistID;
    long count;
  
    std::vector<std::vector<std::complex<double> > > globalcoeff;

    struct sphericalsum
    {
      sphericalsum(const dynamo::Simulation * const, 
		   const double&, const size_t&);
    
      void operator()(const Particle&, const size_t&) const;
    
      void clear();

      const dynamo::Simulation* const Sim;
      const double rg;
      const size_t maxl;
      mutable size_t count;
      mutable std::vector<std::vector<std::complex<double> > > coeffsum;
    };
  };
}
