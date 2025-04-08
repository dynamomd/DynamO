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

#include <boost/algorithm/string/predicate.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/include.hpp>
#include <dynamo/outputplugins/tickerproperty/structureImage.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
OPStructureImaging::OPStructureImaging(const dynamo::Simulation *tmp,
                                       const magnet::xml::Node &XML)
    : OPTicker(tmp, "StructureImaging"), id(0), imageCount(500) {
  operator<<(XML);
}

void OPStructureImaging::operator<<(const magnet::xml::Node &XML) {
  if (!XML.hasAttribute("Structure"))
    M_throw() << "You must specify the name of the structure to monitor for "
                 "StructureImaging";

  structureName = XML.getAttribute("Structure");

  if (XML.hasAttribute("MaxImages"))
    imageCount = XML.getAttribute("MaxImages").as<size_t>();
}

void OPStructureImaging::initialise() {
  dout << "Initialising Structure imaging with a max of " << imageCount
       << " snapshots" << std::endl;

  id = Sim->topology.size();

  for (const shared_ptr<Topology> &ptr : Sim->topology)
    if (boost::iequals(structureName, ptr->getName()))
      id = ptr->getID();

  if (id == Sim->topology.size())
    M_throw() << "Could not find a structure named " << structureName
              << " in the simulation";

  imagelist.clear();
  ticker();
}

//  void
//  OPStructureImaging::replicaExchange(OutputPlugin& nplug)
//  {
//    std::swap(Sim, static_cast<OPStructureImaging&>(nplug).Sim);
//  }

void OPStructureImaging::ticker() {
  if (imageCount != 0) {
    --imageCount;
    printImage();
  }
}

void OPStructureImaging::printImage() {
  for (const shared_ptr<IDRange> &prange : Sim->topology[id]->getMolecules()) {
    std::vector<Vector> atomDescription;

    Vector lastpos(Sim->particles[*prange->begin()].getPosition());

    Vector masspos{0, 0, 0};

    double sysMass(0.0);

    Vector sumrij{0, 0, 0};

    for (const size_t &pid : *prange) {
      // This is all to make sure we walk along the structure
      const Particle &part(Sim->particles[pid]);
      Vector rij = part.getPosition() - lastpos;
      lastpos = part.getPosition();
      Sim->BCs->applyBC(rij);

      sumrij += rij;

      double pmass = Sim->species[part]->getMass(pid);
      sysMass += pmass;
      masspos += sumrij * pmass;

      atomDescription.push_back(sumrij);
    }

    masspos /= sysMass;

    for (Vector &rijpos : atomDescription)
      rijpos -= masspos;

    // We use a trick to not copy the vector, just swap it over
    imagelist.push_back(std::vector<Vector>());
    std::swap(imagelist.back(), atomDescription);
  }
}

void OPStructureImaging::output(magnet::xml::XmlStream &XML) {
  XML << magnet::xml::tag("StructureImages") << magnet::xml::attr("version")
      << 2;

  for (const std::vector<Vector> &vec : imagelist) {
    XML << magnet::xml::tag("Image");

    size_t id(0);

    for (const Vector &vec2 : vec) {
      XML << magnet::xml::tag("Atom") << magnet::xml::attr("ID") << id++
          << (vec2 / Sim->units.unitLength()) << magnet::xml::endtag("Atom");
    }

    XML << magnet::xml::endtag("Image");
  }

  XML << magnet::xml::endtag("StructureImages");
}
} // namespace dynamo
