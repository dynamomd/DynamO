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
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/interactions/captures.hpp>
#include <map>
#include <vector>
#include <tr1/unordered_map>

namespace dynamo {
  namespace detail {
    struct OPContactMapHash
    {
      ::std::size_t operator()(const std::vector<ICapture::value_type>& map) const;
    };

    struct OPContactMapValueHash
    {
      ::std::size_t operator()(const ICapture::value_type& pair) const;
    };

    struct OPContactMapPairHash
    {
      ::std::size_t operator()(const std::pair<std::size_t, std::size_t> & pair) const;
    };
  }
  class OPContactMap: public OutputPlugin
  {
  public:
    OPContactMap(const dynamo::Simulation*, const magnet::xml::Node&);

    ~OPContactMap() {}

    virtual void eventUpdate(const IntEvent&, const PairEventData&);
    virtual void eventUpdate(const GlobalEvent&, const NEventData&);
    virtual void eventUpdate(const LocalEvent&, const NEventData&);
    virtual void eventUpdate(const System&, const NEventData&, const double&);

    virtual void initialise();

    virtual void output(magnet::xml::XmlStream&);

    virtual void operator<<(const magnet::xml::Node&);

    virtual void changeSystem(OutputPlugin*);

  private:
    void stream(double);
    void flush();
    
    double _weight;
    double _total_weight;
    /*! \brief A sorted listing of all the captured pairs in the
     system
    */
    std::map<ICapture::key_type, ICapture::mapped_type> _current_map;
    size_t _next_map_id;

    struct MapData
    {
      MapData(double energy=0, size_t id = 0): _weight(0), _energy(energy), _id(id) {}
      double _weight;
      double _energy;
      size_t _id;
    };

    typedef std::vector<ICapture::value_type> MapKey;
    /*! \brief A hash table storing the histogram of the contact maps.
      
      The key of this map is a sorted list of the captured pairs in
      the system. This sorting is implicitly carried out by the
      _current_map.
     */
    std::tr1::unordered_map<MapKey, MapData,  detail::OPContactMapHash> _collected_maps;
    std::tr1::unordered_map<std::pair<size_t, size_t>, size_t, detail::OPContactMapPairHash> _map_links;
  };
}
