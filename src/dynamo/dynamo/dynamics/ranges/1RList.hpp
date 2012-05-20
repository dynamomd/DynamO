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
#include <dynamo/dynamics/ranges/1range.hpp>
#include <dynamo/simulation/particle.hpp>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>
#include <vector>

namespace dynamo {
  class RList: public Range
  {
  public:
    RList(const magnet::xml::Node& XML) 
    { operator<<(XML); }

    RList() {}

    virtual bool isInRange(const Particle &part) const
    {
      BOOST_FOREACH(const unsigned long ID, IDs)
	if (part.getID() == ID)
	  return true;
      return false;
    }

    //The data output classes
    virtual void operator<<(const magnet::xml::Node& XML)
    {
      if (strcmp(XML.getAttribute("Range"),"List"))
	M_throw() << "Attempting to load RList from non list";
      try {
    
	for (magnet::xml::Node node = XML.fastGetNode("ID"); node.valid(); ++node)
	  IDs.push_back(node.getAttribute("val").as<unsigned long>());
      }
      catch (boost::bad_lexical_cast &)
	{
	  M_throw() << "Failed a lexical cast in RList";
	}
    }
    
    virtual unsigned long size() const { return IDs.size(); };

    virtual iterator begin() const { return Range::iterator(0, this); }

    virtual iterator end() const { return Range::iterator(IDs.size(), this); }

    virtual unsigned long operator[](unsigned long i) const { return IDs[i]; }

    virtual unsigned long at(unsigned long i) const { return IDs.at(i); }

  protected:
    virtual const unsigned long& getIteratorID(const unsigned long &i) const { return IDs[i]; }

    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Range") << "List";
      BOOST_FOREACH(unsigned long ID, IDs)
	XML << magnet::xml::tag("ID") << magnet::xml::attr("val") << ID << magnet::xml::endtag("ID");
    }

    std::vector<unsigned long> IDs;
  };
}
