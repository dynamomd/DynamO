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

#include <dynamo/outputplugins/contactmap.hpp>
#include <dynamo/interactions/captures.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace {
  ::std::size_t
  hash_combine(const ::std::size_t hash1, const ::std::size_t hash2)
  {
    return hash1 ^ (hash2 + 0x9e3779b9 + (hash1 << 6) + (hash1 >> 2));
  };  
}

namespace dynamo {
  namespace detail {
    ::std::size_t
    OPContactMapPairHash::operator()(const std::pair<std::size_t, std::size_t>& entry) const
    {
      return hash_combine(entry.first, entry.second);
    };
  }

  OPContactMap::OPContactMap(const dynamo::Simulation* t1, const magnet::xml::Node& XML):
    OutputPlugin(t1,"ContactMap", 1) //This plugin should be updated and initialised after the misc plugin
  { operator<<(XML); }

  void 
  OPContactMap::initialise() 
  {
    _collected_maps.clear();
    _map_links.clear();
    _next_map_id = 0;
    _weight = 0;
    _total_weight = 0;
  
    _interaction = std::dynamic_pointer_cast<ICapture>(Sim->interactions[_interaction_name]);

    if (!_interaction)
      M_throw() << "Could not cast \"" << _interaction_name << "\" to an ICapture type to build the contact map";
    
    _current_map = _collected_maps.insert(CollectedMapType::value_type(*_interaction, MapData(Sim->calcInternalEnergy(), _next_map_id++))).first;
  }

  void OPContactMap::stream(double dt) { _weight += dt; }

  void 
  OPContactMap::flush()
  {
    //Cannot create new maps here, as flush may happen when the output plugins are invalid
    MapData& data = _current_map->second;
    data._weight += _weight;
    _total_weight += _weight;
    _weight = 0;
  }

  void 
  OPContactMap::operator<<(const magnet::xml::Node& XML) {
    if (!XML.hasAttribute("Interaction"))
      M_throw() << "You must specify an Interaction name for ContactMap tracking using the Interaction option";
    
    _interaction_name = XML.getAttribute("Interaction").as<std::string>();
  }

  void 
  OPContactMap::eventUpdate(const IntEvent &event, const PairEventData &eventdata) 
  {
    stream(event.getdt());

    if (event.getInteractionID() == _interaction->getID())
      if ((event.getType() == STEP_IN) || (event.getType() == STEP_OUT))
	mapChanged(true);
  }

  void 
  OPContactMap::mapChanged(bool addLink) {
    flush();
    size_t oldMapID(_current_map->second._id);
    
    //Try and find the current map in the collected maps
    _current_map = _collected_maps.find(*_interaction);
    if (_current_map == _collected_maps.end())
      //Insert the new map
      _current_map = _collected_maps.insert(CollectedMapType::value_type(*_interaction, MapData(Sim->getOutputPlugin<OPMisc>()->getConfigurationalU(), _next_map_id++))).first;
    
    //Add the link	    
    if (addLink)
      ++(_map_links[std::make_pair(oldMapID, _current_map->second._id)]);
  }

  void 
  OPContactMap::replicaExchange(OutputPlugin& otherplugin)
  {
    OPContactMap& other_map = static_cast<OPContactMap&>(otherplugin);

    //Swap over the sim pointers
    std::swap(Sim, other_map.Sim);
    std::swap(_interaction, other_map._interaction);

    //Now let each plugin know the map has changed
    mapChanged(false);
    other_map.mapChanged(false);
  }

  void 
  OPContactMap::eventUpdate(const GlobalEvent &event, const NEventData&) 
  { stream(event.getdt()); }

  void 
  OPContactMap::eventUpdate(const LocalEvent &event, const NEventData&) 
  { stream(event.getdt()); }

  void 
  OPContactMap::eventUpdate(const System&, const NEventData&, const double& dt)
  { stream(dt); }

  void
  OPContactMap::periodicOutput()
  {
    I_Pcout() << ", Maps " << _collected_maps.size() << ", links " << _map_links.size();
  }

  void 
  OPContactMap::output(magnet::xml::XmlStream& XML)
  {
    dout << "Writing out " << _collected_maps.size() << " Contact maps with "
	 << _map_links.size() << " links" << std::endl;

    XML << magnet::xml::tag("ContactMap")
	<< magnet::xml::tag("Maps")
      	<< magnet::xml::attr("Count") << _collected_maps.size();
    ;
    
    for (const CollectedMapType::value_type& entry : _collected_maps)
      {
	XML << magnet::xml::tag("Map")
	    << magnet::xml::attr("ID") << entry.second._id
	    << magnet::xml::attr("Energy") << entry.second._energy / Sim->units.unitEnergy()
	    << magnet::xml::attr("Weight") << entry.second._weight / _total_weight;
	
	for (const ICapture::value_type& ids : entry.first)
	  XML << magnet::xml::tag("Contact")
	      << magnet::xml::attr("ID1") << ids.first.first
	      << magnet::xml::attr("ID2") << ids.first.second
	      << magnet::xml::attr("State") << ids.second
	      << magnet::xml::endtag("Contact");
	
	XML << magnet::xml::endtag("Map");
      }

    XML << magnet::xml::endtag("Maps")
	<< magnet::xml::tag("Links")
      	<< magnet::xml::attr("Count") << _map_links.size();

    for (const LinksMapType::value_type& entry : _map_links)
      XML << magnet::xml::tag("Link")
	  << magnet::xml::attr("Source") << entry.first.first
	  << magnet::xml::attr("Target") << entry.first.second
	  << magnet::xml::attr("Occurrences") << entry.second
	  << magnet::xml::endtag("Link");

    XML << magnet::xml::endtag("Links")
	<< magnet::xml::endtag("ContactMap");
  }
}
