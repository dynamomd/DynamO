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
#include <dynamo/base.hpp>
#include <dynamo/1particleEventData.hpp>
#include <magnet/math/vector.hpp>

namespace magnet { namespace xml { class XmlStream; } }

namespace dynamo {
  class Particle;
  class Event;

  namespace boundary {
    class Object : public dynamo::SimBase {
    public:
      Object(Simulation* const SD, const std::string aName);

      virtual bool validateState(const Particle& part, bool textoutput = true) const = 0;
      virtual Event getLinearEvent(const Particle& part, const double diameter, const Vector origin, const Vector velocity) const = 0;

      virtual Event getOscillatingEvent(const Particle& part, const double diameter, const Vector origin, const Vector amplitude, const double freq, const double t_shift) const {
	M_throw() << "Not implemented";
      }
    
      virtual Vector getContactNormal(const Particle&, const Event&) const = 0;
  
      virtual void outputXML(magnet::xml::XmlStream&) const = 0;
    
      static shared_ptr<Object> getClass(const magnet::xml::Node&, dynamo::Simulation*);
    };

    class PlanarWall : public Object {
    public:
      PlanarWall(const magnet::xml::Node&, dynamo::Simulation*);

      virtual bool validateState(const Particle& part, bool textoutput = true) const;

      virtual Event getLinearEvent(const Particle& part, const double diameter, const Vector origin, const Vector velocity) const;
      
      virtual Vector getContactNormal(const Particle&, const Event&) const;

      virtual void outputXML(magnet::xml::XmlStream&) const;

    protected:
      Vector _position;
      Vector _normal;
    };
  }
}
