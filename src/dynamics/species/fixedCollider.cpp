/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "fixedCollider.hpp"
#include <boost/foreach.hpp>
#include "../../base/is_simdata.hpp"
#include <boost/tokenizer.hpp>

void 
SpFixedCollider::initialise()
{
  Species::initialise();

  BOOST_FOREACH(size_t ID, *range)
    Sim->particleList[ID].clearState(Particle::DYNAMIC);
}

void 
SpFixedCollider::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
  
  try {
    mass = 0;
    spName = XML.getAttribute("Name");
    intName = XML.getAttribute("IntName");

    if (XML.isAttributeSet("Color"))
      {
	typedef boost::tokenizer<boost::char_separator<char> >
	  Tokenizer;
	
	boost::char_separator<char> colorSep(",");
	
	std::string data(XML.getAttribute("Color"));

	Tokenizer tokens(data, colorSep);
	Tokenizer::iterator value_iter = tokens.begin();

	if (value_iter == tokens.end())
	  throw std::runtime_error("Malformed color in species");
	_constColor[0] = boost::lexical_cast<int>(*value_iter);
	
	if (++value_iter == tokens.end())
	  throw std::runtime_error("Malformed color in species");
	_constColor[1] = boost::lexical_cast<int>(*value_iter);
	
	if (++value_iter == tokens.end())
	  throw std::runtime_error("Malformed color in species");
	_constColor[2] = boost::lexical_cast<int>(*value_iter);
	
	if (++value_iter != tokens.end())
	  throw std::runtime_error("Malformed color in species");
	
	_constColor[3] = 255;
	
	_colorMode = CONSTANT;
      }    
  } 
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in SpFixedCollider";
    }
}

void 
SpFixedCollider::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Name") << spName
      << xml::attr("IntName") << intName
      << xml::attr("Type") << "FixedCollider";

  if (_colorMode == CONSTANT)
    {
      std::string colorval = boost::lexical_cast<std::string>(_constColor[0] + 0);
      colorval += ",";
      colorval += boost::lexical_cast<std::string>(_constColor[1] + 0);
      colorval += ",";
      colorval += boost::lexical_cast<std::string>(_constColor[2] + 0);
      XML << xml::attr("Color") << colorval;
    }
  
  XML << range;
}
