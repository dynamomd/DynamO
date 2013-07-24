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
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/particle.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <vector>

namespace dynamo {
  class IDRangeUnion: public IDRange
  {
  public:
    IDRangeUnion(const magnet::xml::Node& XML, const dynamo::Simulation* Sim) 
    { 
      for (magnet::xml::Node node = XML.fastGetNode("IDRange"); node.valid(); ++node)
	{
	  shared_ptr<IDRange> ptr(IDRange::getClass(node, Sim));
	  ranges.push_back(ptr);
	}
    }

    IDRangeUnion() {}

    virtual bool isInRange(const Particle &part) const
    {
      for (const shared_ptr<IDRange>& r : ranges)
	if (r->isInRange(part)) return true;

      return false;
    }

    virtual unsigned long size() const 
    { 
      size_t size(0);
      for (const shared_ptr<IDRange>& r : ranges)
	size += r->size();
      return size;
    }

    virtual unsigned long operator[](unsigned long i) const 
    { 
      for (const shared_ptr<IDRange>& r : ranges)
	{
	  if (i < r->size())
	    return r->operator[](i);
	  i -= r->size();
	}
      M_throw() << "Out of range access";
    }

    virtual unsigned long at(unsigned long i) const 
    { 
      return operator[](i);
    }

  protected:

    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Type") << "Union";
      for (const shared_ptr<IDRange>& r : ranges)
	XML << *r;
    }

    std::vector<shared_ptr<IDRange> > ranges;
  };
}
