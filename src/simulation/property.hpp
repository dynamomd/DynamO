/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#include <vector>
#include <string>
#include <algorithm>
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>

class Property
{
public:
  Property(std::string name, size_t N):
    _name(name),
    _storage(N)
  {}
  
  inline double& getProperty(size_t ID) { return _storage[ID]; }
  inline const double& getProperty(size_t ID) const { return _storage[ID]; }
 
  bool operator==(const std::string& str) const { return _name == str; }

  void XMLOutput(xml::XmlStream& XML, const size_t ID) const
  {
    XML << xml::attr(_name) << _storage[ID];
  }

protected:
  std::string _name;
  std::vector<double> _storage;
};

class PropertyStore : public std::vector<Property>
{
public:
  const Property& getProperty(const std::string& name) const
  {
    const_iterator it = std::find(begin(), end(), name);
    
    if (it == end()) M_throw() << "Could not find the property";
    
    return *it;
  }

  void XMLOutput(xml::XmlStream& XML, const size_t ID) const
  {
    for (const_iterator iPtr = begin(); iPtr != end(); ++iPtr)
      iPtr->XMLOutput(XML, ID);
  }
  
private:
};
