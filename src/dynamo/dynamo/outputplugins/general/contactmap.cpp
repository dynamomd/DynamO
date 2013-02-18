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

#include <dynamo/outputplugins/general/contactmap.hpp>
#include <dynamo/interactions/captures.hpp>
#include <dynamo/outputplugins/0partproperty/misc.hpp>
#include <dynamo/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  namespace detail {
    ::std::size_t
    OPContactMapPairHash::operator()(const std::pair<size_t, size_t>& pair) const
    {
      return pair.first ^ (pair.second + 0x9e3779b9 + (pair.first << 6) + (pair.first >> 2));
    };

    ::std::size_t
    OPContactMapHash::operator()(const std::vector<std::pair<size_t, size_t> >& map) const
    {
      OPContactMapPairHash pairhash;
      ::std::size_t hash(0);
      typedef ::std::pair< ::std::size_t, ::std::size_t> Entry;
      BOOST_FOREACH(const Entry& entry, map)
	hash = pairhash(std::make_pair(pairhash(std::make_pair(hash, entry.first)), entry.second));
      return hash;
    }
  }

  OPContactMap::OPContactMap(const dynamo::Simulation* t1, 
				   const magnet::xml::Node& XML):
    OutputPlugin(t1,"ContactMap", 1) //This plugin should be updated and initialised after the misc plugin
  { operator<<(XML); }

  void 
  OPContactMap::initialise() 
  {
    _current_map.clear();
    _collected_maps.clear();
    _map_links.clear();
    _next_map_id = 0;
    _weight = 0;
    _total_weight = 0;
  
    BOOST_FOREACH(const shared_ptr<Interaction>& interaction, Sim->interactions)
      {
	shared_ptr<ISingleCapture> capture_interaction = std::tr1::dynamic_pointer_cast<ISingleCapture>(interaction);
	typedef std::pair<size_t, size_t> key;
	if (capture_interaction)
	  _current_map.insert(capture_interaction->getMap().begin(), capture_interaction->getMap().end());
      }
    
    MapKey key(_current_map.begin(), _current_map.end());
    _collected_maps.insert(std::pair<const MapKey, MapData>(key, MapData(Sim->calcInternalEnergy(), _next_map_id++)));
  }

  void OPContactMap::stream(double dt) { _weight += dt; }

  void 
  OPContactMap::flush() 
  {
    //Cannot create new maps here, as flush may happen when the output plugins are invalid
    MapKey key(_current_map.begin(), _current_map.end());
    MapData& data = _collected_maps[key];
    data._weight += _weight;
    _total_weight += _weight;
    _weight = 0;
    //If you change anything here, you'll have to change the function
    //below too:
    //eventUpdate(const IntEvent &event, const PairEventData &data)
  }

  void OPContactMap::operator<<(const magnet::xml::Node& XML) {}

  void 
  OPContactMap::eventUpdate(const IntEvent &event, const PairEventData &eventdata) 
  {
    stream(event.getdt());

    if ((event.getType() == STEP_IN) || (event.getType() == STEP_OUT))
      {
	shared_ptr<ISingleCapture> capture_interaction = std::tr1::dynamic_pointer_cast<ISingleCapture>(Sim->interactions[event.getInteractionID()]);
	if (capture_interaction)
	  {
	    //Cache the old map data, and flush the entry
	    MapKey oldMapKey(_current_map.begin(), _current_map.end());
	    MapData& olddata = _collected_maps[oldMapKey];
	    olddata._weight += _weight;
	    _total_weight += _weight;
	    _weight = 0;
	    size_t oldMapID(olddata._id);

	    //Update the map
	    size_t minID = std::min(event.getParticle1ID(), event.getParticle2ID());
	    size_t maxID = std::max(event.getParticle1ID(), event.getParticle2ID());
	    if (capture_interaction->isCaptured(minID, maxID))
	      _current_map.insert(std::make_pair(minID, maxID));
	    else
	      _current_map.erase(std::make_pair(minID, maxID));

	    //Check if the new map is already in the list.  If the
	    //current map is not, insert it, and initialise its ID
	    MapKey newkey(_current_map.begin(), _current_map.end());
	    std::tr1::unordered_map<MapKey, MapData,  detail::OPContactMapHash>::iterator _map_it
	      = _collected_maps.find(newkey);
	    if (_map_it == _collected_maps.end())
	      _map_it = _collected_maps.insert(std::pair<const MapKey, MapData>(newkey, MapData(Sim->getOutputPlugin<OPMisc>()->getConfigurationalU(), _next_map_id++))).first;
	    
	    //Add the link	    
	    ++(_map_links[std::make_pair(oldMapID, _map_it->second._id)]);
	  }
      }
  }

  void 
  OPContactMap::changeSystem(OutputPlugin* otherplugin)
  {
    OPContactMap& other_map = static_cast<OPContactMap&>(*otherplugin);

    //First, flush both maps, this handles the weights.
    flush();
    other_map.flush();
        
    MapKey key1(_current_map.begin(), _current_map.end());
    MapKey key2(other_map._current_map.begin(), other_map._current_map.end());

    //Check that this plugin has the map from the other in its collection
    if (_collected_maps.find(key2) == _collected_maps.end())
      _collected_maps.insert(std::make_pair(key2, MapData(other_map._collected_maps[key2]._energy, _next_map_id++)));
    //And vice versa
    if (other_map._collected_maps.find(key1) == other_map._collected_maps.end())
      other_map._collected_maps.insert(std::make_pair(key1, MapData(_collected_maps[key1]._energy, other_map._next_map_id++)));

    //Now swap the current contact maps
    std::swap(_current_map, other_map._current_map);

    //Now swap over the sim pointers
    std::swap(Sim, other_map.Sim);
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
  OPContactMap::output(magnet::xml::XmlStream& XML)
  {
    flush();
    XML << magnet::xml::tag("ContactMap")
	<< magnet::xml::tag("Maps")
      	<< magnet::xml::attr("Count") << _collected_maps.size();
      ;

    typedef std::pair<MapKey, MapData> MapDataType;
    BOOST_FOREACH(const MapDataType& entry, _collected_maps)
      {
	XML << magnet::xml::tag("Map")
	    << magnet::xml::attr("ID") << entry.second._id
	    << magnet::xml::attr("Energy") << entry.second._energy / Sim->units.unitEnergy()
	    << magnet::xml::attr("Weight") << entry.second._weight / _total_weight;
	
	typedef std::pair<size_t, size_t> IDPair;
	BOOST_FOREACH(const IDPair& ids, entry.first)
	  XML << magnet::xml::tag("Contact")
	      << magnet::xml::attr("ID1") << ids.first
	      << magnet::xml::attr("ID2") << ids.second
	      << magnet::xml::endtag("Contact");
	
	XML << magnet::xml::endtag("Map");
      }

    XML << magnet::xml::endtag("Maps")
	<< magnet::xml::tag("Links")
      	<< magnet::xml::attr("Count") << _map_links.size();


    typedef std::pair<const std::pair<size_t, size_t>, size_t> EdgeDataType;
    BOOST_FOREACH(const EdgeDataType& entry, _map_links)
      XML << magnet::xml::tag("Link")
	  << magnet::xml::attr("Source") << entry.first.first
	  << magnet::xml::attr("Target") << entry.first.second
	  << magnet::xml::attr("Occurances") << entry.second
	  << magnet::xml::endtag("Link");

    XML << magnet::xml::endtag("Links")
	<< magnet::xml::endtag("ContactMap");
  }
}
