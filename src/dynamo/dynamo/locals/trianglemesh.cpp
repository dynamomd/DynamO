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

#include <dynamo/locals/trianglemesh.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/locals/localEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>


namespace dynamo {
  LTriangleMesh::LTriangleMesh(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
    Local(tmp, "LocalWall")
  { operator<<(XML); }

  LocalEvent 
  LTriangleMesh::getEvent(const Particle& part) const
  {
#ifdef ISSS_DEBUG
    if (!Sim->dynamics->isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    size_t triangleid = 0; //The id of the triangle for which the event is for
    double diam = 0.5 * _diameter->getProperty(part.getID());

    std::pair<double, size_t> tmin(HUGE_VAL, 0); //Default to no collision

    for (size_t id(0); id < _elements.size(); ++id)
      {
	std::pair<double, size_t> t = Sim->dynamics->getSphereTriangleEvent(part,
				  _vertices[_elements[id].get<0>()],
				  _vertices[_elements[id].get<1>()],
				  _vertices[_elements[id].get<2>()],
				  diam);
	if (t < tmin) { tmin = t; triangleid = id; }
      }

    return LocalEvent(part, tmin.first, WALL, *this, 8 * triangleid + tmin.second);
  }

  void
  LTriangleMesh::runEvent(Particle& part, const LocalEvent& iEvent) const
  { 
    ++Sim->eventCount;
  
    const size_t triangleID = iEvent.getExtraData() / Dynamics::T_COUNT;
    const size_t trianglepart = iEvent.getExtraData() % Dynamics::T_COUNT;

    const TriangleElements& elem = _elements[triangleID];
  
    const Vector& A(_vertices[elem.get<0>()]);
    const Vector& B(_vertices[elem.get<1>()]);
    const Vector& C(_vertices[elem.get<2>()]);
    
    //Run the collision and catch the data
    Vector normal;
    switch (trianglepart)
      {
      case Dynamics::T_FACE:
	{
	  normal = (B - A) ^ (C - B);
	  normal /= normal.nrm();
	  break;
	}
      case Dynamics::T_A_CORNER: 
	{ 
	  normal = part.getPosition() - A; 
	  Sim->BCs->applyBC(normal);
	  normal /= normal.nrm();
	  break; 
	}
      case Dynamics::T_B_CORNER: 
	{ 
	  normal = part.getPosition() - B; 
	  Sim->BCs->applyBC(normal);
	  normal /= normal.nrm();
	  break; 
	}
      case Dynamics::T_C_CORNER: 
	{ 
	  normal = part.getPosition() - C; 
	  Sim->BCs->applyBC(normal);
	  normal /= normal.nrm();
	  break; 
	}
      case Dynamics::T_AB_EDGE:
	{
	  Vector edge = B - A;
	  edge /= edge.nrm();

	  normal = part.getPosition() - A;
	  Sim->BCs->applyBC(normal);
	  normal -= (normal | edge) * edge;	
	  normal /= normal.nrm();
	  break; 	
	}
      case Dynamics::T_AC_EDGE:
	{
	  Vector edge = C - A;
	  edge /= edge.nrm();

	  normal = part.getPosition() - A;
	  Sim->BCs->applyBC(normal);
	  normal -= (normal | edge) * edge;	
	  normal /= normal.nrm();
	  break; 	
	}
      case Dynamics::T_BC_EDGE:
	{
	  Vector edge = B - C;
	  edge /= edge.nrm();

	  normal = part.getPosition() - C;
	  Sim->BCs->applyBC(normal);
	  normal -= (normal | edge) * edge;	
	  normal /= normal.nrm();
	  break; 	
	}
      default:
	M_throw() << "Unhandled triangle sphere intersection type encountered";
      }

    NEventData EDat(Sim->dynamics->runPlaneEvent(part, normal, _e->getProperty(part.getID()), 0.0));

    Sim->signalParticleUpdate(EDat);

    //Now we're past the event update the scheduler and plugins
    Sim->ptrScheduler->fullUpdate(part);
  
    BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
  }

  void 
  LTriangleMesh::operator<<(const magnet::xml::Node& XML)
  {
    range = shared_ptr<IDRange>(IDRange::getClass(XML.getNode("IDRange"), Sim));
  
    try {
      _diameter = Sim->_properties.getProperty(XML.getAttribute("Diameter"),
					       Property::Units::Length());
      _e = Sim->_properties.getProperty(XML.getAttribute("Elasticity"),
					Property::Units::Dimensionless());

      localName = XML.getAttribute("Name");

      {//Load the vertex coordinates
	std::istringstream is(XML.getNode("Vertices").getValue());
	is.exceptions(std::ostringstream::badbit | std::ostringstream::failbit);
	is.peek(); //Set the eof flag if needed
	Vector tmp;
	while (!is.eof())
	  {
	    is >> tmp[0];
	    if (is.eof()) M_throw() << "The vertex coordinates is not a multiple of 3";

	    is >> tmp[1];
	    if (is.eof()) M_throw() << "The vertex coordinates is not a multiple of 3";

	    is >> tmp[2];	  
	    _vertices.push_back(tmp * Sim->units.unitLength());
	  }
      }

      {//Load the triangle elements
	std::istringstream is(XML.getNode("Elements").getValue());
	is.exceptions(std::ostringstream::badbit | std::ostringstream::failbit);
	is.peek(); //Set the eof flag if needed

	TriangleElements tmp;
	while (!is.eof())
	  {
	    is >> tmp.get<0>();
	    if (is.eof()) M_throw() << "The triangle elements are not a multiple of 3";
	  
	    is >> tmp.get<1>();;
	    if (is.eof()) M_throw() << "The triangle elements are not a multiple of 3";

	    is >> tmp.get<2>();;

	    if ((tmp.get<0>() >= _vertices.size()) 
		|| (tmp.get<1>() >= _vertices.size()) 
		|| (tmp.get<2>() >= _vertices.size()))
	      M_throw() << "Triangle " << _elements.size() << " has an out of range vertex ID";

	    Vector normal
	      = (_vertices[tmp.get<1>()] - _vertices[tmp.get<0>()])
	      ^ (_vertices[tmp.get<2>()] - _vertices[tmp.get<1>()]);

	    if (normal.nrm() == 0) 
	      M_throw() << "Triangle " << _elements.size() << " has a zero normal!";


	    _elements.push_back(tmp);
	  }
      }

    }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in LTriangleMesh";
      }
  }

  void 
  LTriangleMesh::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "TriangleMesh" 
	<< magnet::xml::attr("Name") << localName
	<< magnet::xml::attr("Elasticity") << _e->getName()
	<< magnet::xml::attr("Diameter") << _diameter->getName()
	<< range;

    XML << magnet::xml::tag("Vertices") << magnet::xml::chardata();
    BOOST_FOREACH(Vector vert, _vertices)
      XML << vert[0] / Sim->units.unitLength() << " " 
	  << vert[1] / Sim->units.unitLength() << " "
	  << vert[2] / Sim->units.unitLength() << "\n";
    XML << magnet::xml::endtag("Vertices");

    XML << magnet::xml::tag("Elements") << magnet::xml::chardata();
    BOOST_FOREACH(TriangleElements elements, _elements)
      XML << elements.get<0>() << " " 
	  << elements.get<1>() << " "
	  << elements.get<2>() << "\n";
    XML << magnet::xml::endtag("Elements");
  }

#ifdef DYNAMO_visualizer

  shared_ptr<coil::RenderObj>
  LTriangleMesh::getCoilRenderObj() const
  {
    const double lengthRescale = 1 / Sim->primaryCellSize.maxElement();

    if (!_renderObj)
      {
	std::vector<float> verts;
	verts.reserve(3 * _vertices.size());
	BOOST_FOREACH(const Vector& v, _vertices)
	  {
	    verts.push_back(v[0] * lengthRescale);
	    verts.push_back(v[1] * lengthRescale);
	    verts.push_back(v[2] * lengthRescale);
	  }

	std::vector<GLuint> elems;
	elems.reserve(3 * _elements.size());
	BOOST_FOREACH(const TriangleElements& e, _elements)
	  {
	    elems.push_back(e.get<0>());
	    elems.push_back(e.get<1>());
	    elems.push_back(e.get<2>());
	  }
      
	_renderObj.reset(new coil::RTriangleMesh(getName(), verts, elems));
      }
  
    return std::tr1::static_pointer_cast<coil::RenderObj>(_renderObj);
  }
#endif
}
