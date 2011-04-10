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
#include <magnet/xmlreader.hpp>
#include <magnet/thread/refPtr.hpp>


/*! A interface class which allows other classes to access a property
 * of a particle.  These properties are looked up by a name, and the
 * value extracted using the ID of a particle. Some properties are
 * just a single fixed value, their name is their value (see
 * NumericProperty). Others are more complicated and use look-up
 * tables or functions. These are usually defined in the PropertyStore
 * and PropertyHandles are used to access them.
 */
class Property
{
public:
  //! Fetch the value of this property for a particle with a certain ID
  inline virtual const double& getProperty(size_t ID) const = 0;
  
  //! Fetch the name of this property
  inline virtual std::string getName() const = 0;
};

class NumericProperty: public Property
{
public:
  inline NumericProperty(double val): _val(val) {}
  
  inline virtual const double& getProperty(size_t ID) const { return _val; }
  inline virtual std::string getName() const { return boost::lexical_cast<std::string>(_val); }
private:
  double _val;
};

/*! This class stores the properties of the particles loaded from the
 * configuration file and hands out reference counting pointers to the
 * properties to other classes when they're requested by name.
 */
class PropertyStore: private std::vector<magnet::thread::RefPtr<Property> >
{
  typedef magnet::thread::RefPtr<Property> Value;
  typedef std::vector<Value> Base;
  
public:
  typedef Base::const_iterator const_iterator;

  /*! Request a handle to a property using an xml attribute containing
    the properties name.  If the name is a numeric type, the look-up
    in the property store will fail but a one-time NumericProperty is
    create. You may then have lines in the configuration file like so

    <Interaction Elasticity="0.9" ... For a fixed value or

    <Interaction Elasticity="e" ... For a lookup of the particle property "e"

    \param name An Attribute containing either the name or the value of a property.
    \return A reference to the property requested or an instance of NumericProperty.
  */
  inline magnet::thread::RefPtr<Property> getProperty(const magnet::xml::Attribute& name)
  {
    //Try name based lookup first
    for (const_iterator iPtr = Base::begin(); iPtr != Base::end(); ++iPtr)
      if ((*iPtr)->getName() == name.getValue())
	return *iPtr;

    //Try name-is-the-value lookup
    try { 
      double value = boost::lexical_cast<double>(name);
      return Value(new NumericProperty(value));
    } catch (boost::bad_lexical_cast&)
      {
	M_throw() << "Could not find the property named by " << name.getPath();
      }
  }

  /*! Method which loads the properties from the XML configuration file. 
    \param node A xml Node at the root DYNAMOconfig Node of the config file.
   */
  inline void loadProperties(const magnet::xml::Node& node)
  {
    if (!node.getNode("Properties").valid()) return;

    for (magnet::xml::Node propNode = node.getNode("Properties").getNode("Property");
	 propNode.valid(); ++propNode)
      M_throw() << "Unsupported Property type, " << propNode.getAttribute("Type");
  }
private:
};
