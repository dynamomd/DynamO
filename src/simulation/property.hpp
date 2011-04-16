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
#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/thread/refPtr.hpp>
#include <magnet/units.hpp>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

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
  typedef magnet::units::Units Units;

  Property(Units units): 
    _units(units) {}

  //! Fetch the value of this property for a particle with a certain ID
  inline virtual const double& getProperty(size_t ID) const { M_throw() << "Unimplemented"; }

  //! Fetch the maximum value of this property
  inline virtual const double& getMaxValue() const { M_throw() << "Unimplemented"; }

  //! This is called whenever a unit is rescaled.
  //!
  //! This function must check the _units of the property and raise
  //! the rescale factor to the correct power.
  //! \param dim The unit that is being rescaled [(L)ength, (T)ime, (M)ass].
  //! \param rescale The factor to rescale the unit by.
  inline virtual const void rescaleUnit(const Units::Dimension dim, 
					const double rescale)
  { M_throw() << "Unimplemented"; }
  
  //! Fetch the name of this property
  inline virtual std::string getName() const { M_throw() << "Unimplemented"; }

  //! Fetch the units of this property
  inline const Units& getUnits() const { return _units; }

  //! Helper to write out derived classes
  friend xml::XmlStream operator<<(xml::XmlStream& XML, const Property& prop)
  { prop.outputXML(XML); return XML; }

protected:
  virtual void outputXML(xml::XmlStream& XML) const = 0;

  //! The Units of the property.
  magnet::units::Units _units;
};

class NumericProperty: public Property
{
public:
  inline NumericProperty(double val, const Property::Units& units):
    Property(units), _val(val) {}
  
  inline virtual const double& getProperty(size_t ID) const { return _val; }
  inline virtual std::string getName() const { return boost::lexical_cast<std::string>(_val); }
  inline virtual const double& getMaxValue() const { return _val; }

  //! \sa Property::rescaleUnit
  inline virtual const void rescaleUnit(const Units::Dimension dim, 
					const double rescale)
  { _val *= std::pow(rescale, _units.getUnitsPower(dim));  }

private:
  virtual void outputXML(xml::XmlStream& XML) const {}

  double _val;
};

/*! This class stores the properties of the particles loaded from the
 * configuration file and hands out reference counting pointers to the
 * properties to other classes when they're requested by name.
 */
class PropertyStore
{
  typedef magnet::thread::RefPtr<Property> Value;
  typedef std::vector<Value> Container;
  
  //!Contains the NumericProperty's that are defined by their
  //!name. These are only stored in the PropertyStore for unit
  //!rescaling.
  Container _numericProperties;

  //!Contains the properties that are looked up by their name.
  Container _namedProperties;

  typedef Container::iterator iterator;

public:
  typedef Container::const_iterator const_iterator;

  /*! Request a handle to a property using a string containing the
    properties name.  If the name is a string representation of a
    numeric type, the look-up in the property store will fail but a
    one-time NumericProperty is created. You may then have lines in
    the configuration file like so

    For a fixed value
    <Interaction Elasticity="0.9" ... 

    or for a lookup in the property store
    <Interaction Elasticity="e" ... For a lookup of the particle property "e"

    \param name An Attribute containing either the name or the value of a property.
    \return A reference to the property requested or an instance of NumericProperty.
  */
  inline magnet::thread::RefPtr<Property> getProperty(const std::string& name,
						      const Property::Units& units)
  {
    try { return getPropertyBase(name, units); }
    catch (boost::bad_lexical_cast&)
      { M_throw() << "Could not find the property named by " << name; }
  }

  /*! Request a handle to a property using an xml attribute containing
    the properties name. See getProperty(const std::string& name) for
    usage info. */
  inline magnet::thread::RefPtr<Property> getProperty(const magnet::xml::Attribute& name,
						      const Property::Units& units)
  {
    try { return getPropertyBase(name.getValue(), units); }
    catch (boost::bad_lexical_cast&)
      { M_throw() << "Could not find the property named by " << name.getPath(); }
  }

  /*! Request a handle to a property, but this specialization always
      returns a new instance of NumericProperty.
      \sa getProperty(const std::string& name)
  */
  inline magnet::thread::RefPtr<Property> getProperty(const double& name, 
						      const Property::Units& units)
  { return Value(new NumericProperty(name, units)); }

  /*! Method which loads the properties from the XML configuration file. 
    \param node A xml Node at the root DYNAMOconfig Node of the config file.
   */
  inline PropertyStore& operator<<(const magnet::xml::Node& node)
  {
    if (!node.getNode("Properties").valid()) return *this;

    for (magnet::xml::Node propNode = node.getNode("Properties").getNode("Property");
	 propNode.valid(); ++propNode)
      M_throw() << "Unsupported Property type, " << propNode.getAttribute("Type");

    return *this;
  }

  inline friend xml::XmlStream operator<<(xml::XmlStream& XML, const PropertyStore& propStore)
  {
    XML << xml::tag("Properties");

    for (const_iterator iPtr = propStore._namedProperties.begin(); 
	 iPtr != propStore._namedProperties.end(); ++iPtr)
      XML << (*(*iPtr));

    XML << xml::endtag("Properties");

    return XML;
  }

  //! Function to rescale the units of all Property-s.
  //!
  //! \param dim The unit that is being rescaled [(L)ength, (T)ime, (M)ass].
  //! \param rescale The factor to rescale the unit by.
  inline const void rescaleUnit(const Property::Units::Dimension dim, 
				const double rescale)
  {  
    for (iterator iPtr = _namedProperties.begin(); 
	 iPtr != _namedProperties.end(); ++iPtr)
      (*iPtr)->rescaleUnit(dim, rescale);

    for (iterator iPtr = _namedProperties.begin(); 
	 iPtr != _namedProperties.end(); ++iPtr)
      (*iPtr)->rescaleUnit(dim, rescale);
  }


private:

  inline magnet::thread::RefPtr<Property> getPropertyBase(const std::string name,
							  const Property::Units& units)
  {
    //Try name based lookup first
    for (const_iterator iPtr = _namedProperties.begin(); 
	 iPtr != _namedProperties.end(); ++iPtr)
      if ((*iPtr)->getName() == name)
	if ((*iPtr)->getUnits() == units)
	  return *iPtr;
	else
	  M_throw() << "Property \"" << name << "\" found with units of " 
		    << std::string((*iPtr)->getUnits())
		    << ", but the requested property has units of " 
		    << std::string(units);
    
    //Try name-is-the-value lookup, if this fails a
    //boost::bad_lexical_cast& will be thrown and must be caught by
    //the caller. 
    return Value(new NumericProperty(boost::lexical_cast<double>(name), units));
  }
};
