#include <dynamo/particle.hpp>

namespace dynamo {
//! \brief Operator to write out an XML representation of a Particle.
magnet::xml::XmlStream &operator<<(magnet::xml::XmlStream &XML,
                                   const Particle &particle) {
  XML << magnet::xml::attr("ID") << particle._ID;

  if (!particle.testState(Particle::DYNAMIC))
    XML << magnet::xml::attr("Static") << "Static";

  XML << magnet::xml::tag("P") << (particle._pos) << magnet::xml::endtag("P")
      << magnet::xml::tag("V") << (particle._vel) << magnet::xml::endtag("V");

  return XML;
}

Particle::Particle(const Vector &position, const Vector &velocity,
                   const unsigned long &nID)
    : _pos(position), _peculiarTime(0.0), _vel(velocity), _ID(nID),
      _state(DEFAULT) {}

Particle::Particle(const magnet::xml::Node &XML, unsigned long nID)
    : _peculiarTime(0.0), _ID(nID), _state(DEFAULT) {
  if (XML.hasAttribute("Static"))
    clearState(DYNAMIC);

  _pos << XML.getNode("P");
  _vel << XML.getNode("V");
}
} // namespace dynamo