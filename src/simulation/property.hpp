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
#include <magnet/thread/refPtr.hpp>


/*! \brief A container to hold a property every particle in the system
 * has.
 *
 * 
 */
class Property : public std::vector<double>
{
public:
  Property(std::string name, size_t N):
    std::vector<double>(N),
    _name(name)
  {}
  
  inline double& getProperty(size_t ID) { return operator[](ID); }
  inline const double& getProperty(size_t ID) const { return operator[](ID); }
 
  bool operator==(const std::string& str) const { return _name == str; }

  void XMLOutput(xml::XmlStream& XML, const size_t ID) const
  { XML << xml::attr(_name) << operator[](ID); }

protected:
  std::string _name;
};

class PropertyStore : public std::vector<magnet::thread::RefPtr<Property> >
{
public:
  const magnet::thread::RefPtr<Property>& getProperty(const std::string& name) const
  {
    const_iterator it = std::find(begin(), end(), name);    
    if (it == end()) M_throw() << "Could not find the property";    
    return *it;
  }

  void XMLOutput(xml::XmlStream& XML, const size_t ID) const
  {
    for (const_iterator iPtr = begin(); iPtr != end(); ++iPtr)
      (*iPtr)->XMLOutput(XML, ID);
  }
  
private:
};


class PropertyHandle
{
public:  
  virtual double getProperty(size_t ID) = 0;
};

class PropertyReference : public PropertyHandle
{
public:  
  virtual double getProperty(size_t ID) 
  { return _ref->getProperty(ID); }

protected:
  const magnet::thread::RefPtr<Property> _ref;
};
