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

#include <dynamo/interactions/squarewell.hpp>

namespace dynamo {
  class IThinThread: public ISquareWell
  {
  public:
    IThinThread(const magnet::xml::Node&, dynamo::Simulation*);
  
    void operator<<(const magnet::xml::Node&);

    virtual bool validateState(const Particle& p1, const Particle& p2, bool textoutput = true) const;

    /*! \brief This capture test returns false as (initially) there are no bridges.
     */
    virtual size_t captureTest(const Particle&, const Particle&) const { return false; }

    virtual IntEvent getEvent(const Particle&, const Particle&) const;
  
    virtual void runEvent(Particle&, Particle&, const IntEvent&) const;
  
    virtual void outputXML(magnet::xml::XmlStream&) const;
  };
}
