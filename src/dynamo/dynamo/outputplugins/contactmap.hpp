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
#include <dynamo/interactions/captures.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <map>
#include <unordered_map>
#include <vector>

namespace dynamo {
namespace detail {
struct OPContactMapPairHash {
  size_t operator()(const std::pair<std::size_t, std::size_t> &entry) const;
};
} // namespace detail

class OPContactMap : public OutputPlugin {
public:
  OPContactMap(const dynamo::Simulation *, const magnet::xml::Node &);

  ~OPContactMap() {}

  virtual void eventUpdate(const Event &, const NEventData &);

  virtual void initialise();

  virtual void output(magnet::xml::XmlStream &);

  virtual void operator<<(const magnet::xml::Node &);

  virtual void replicaExchange(OutputPlugin &);

  void periodicOutput();

private:
  void stream(double);
  void flush();

  void mapChanged(bool addLink);

  double _weight;
  double _total_weight;
  /*! \brief A sorted listing of all the captured pairs in the
   system
  */
  size_t _next_map_id;

  struct MapData {
    MapData(double discovery_time = 0, double energy = 0, size_t id = 0)
        : _weight(0), _energy(energy), _discovery_time(discovery_time),
          _id(id) {}
    double _weight;
    double _energy;
    double _discovery_time;
    size_t _id;
  };

  typedef std::unordered_map<detail::CaptureMapKey, MapData,
                             detail::CaptureMapKeyHash>
      CollectedMapType;
  typedef std::unordered_map<std::pair<size_t, size_t>, size_t,
                             detail::OPContactMapPairHash>
      LinksMapType;
  /*! \brief A hash table storing the histogram of the contact maps.

    The key of this map is a sorted list of the captured pairs in
    the system. This sorting is implicitly carried out by the
    _current_map.
   */
  CollectedMapType _collected_maps;
  CollectedMapType::iterator _current_map;
  LinksMapType _map_links;
  std::string _interaction_name;
  std::shared_ptr<ICapture> _interaction;
};
} // namespace dynamo
