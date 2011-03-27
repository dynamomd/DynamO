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

protected:
  std::string _name;
  std::vector<double> _storage;
};

class PropertyStore
{
public:

  Property& getProperty(const std::string& name)
  {
    std::vector<Property>::iterator it 
      = std::find(_properties.begin(), _properties.end(), name);

    if (it == _properties.end())
      M_throw() << "Could not find the property";

    return *it;
  }
  
private:
  std::vector<Property> _properties;
};
