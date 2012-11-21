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
#include <boost/foreach.hpp>
#include <vector>

namespace dynamo {
  class IDRangeList: public IDRange
  {
  public:
    IDRangeList(const magnet::xml::Node& XML) 
    { 
      try {
	for (magnet::xml::Node node = XML.fastGetNode("ID"); node.valid(); ++node)
	  IDs.push_back(node.getAttribute("val").as<size_t>());
      }
      catch (boost::bad_lexical_cast &)
	{
	  M_throw() << "Failed a lexical cast in IDRangeList";
	}
    }

    template<class T>
    IDRangeList(const std::vector<T>& data):
      IDs(data)
    {}

    IDRangeList() {}

    virtual bool isInRange(const Particle &part) const
    {
      BOOST_FOREACH(const unsigned long ID, IDs)
	if (part.getID() == ID)
	  return true;
      return false;
    }

    std::vector<size_t>& getContainer() { return IDs; }
    
    virtual unsigned long size() const { return IDs.size(); }

    virtual unsigned long operator[](unsigned long i) const { return IDs[i]; }

    virtual unsigned long at(unsigned long i) const { return IDs.at(i); }

  protected:

    virtual void outputXML(magnet::xml::XmlStream& XML) const
    {
      XML << magnet::xml::attr("Type") << "List";
      BOOST_FOREACH(unsigned long ID, IDs)
	XML << magnet::xml::tag("ID") << magnet::xml::attr("val") << ID << magnet::xml::endtag("ID");
    }

    std::vector<size_t> IDs;
  };
}
