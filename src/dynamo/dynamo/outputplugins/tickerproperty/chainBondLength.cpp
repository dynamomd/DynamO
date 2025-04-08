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

#include <dynamo/include.hpp>
#include <dynamo/outputplugins/tickerproperty/chainBondLength.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/topology/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <vector>

namespace dynamo {
OPChainBondLength::Cdata::Cdata(size_t ID, size_t CL) : chainID(ID) {
  BondLengths.resize(CL - 1, magnet::math::Histogram<>(0.0001));
}

OPChainBondLength::OPChainBondLength(const dynamo::Simulation *tmp,
                                     const magnet::xml::Node &)
    : OPTicker(tmp, "ChainBondLength") {}

void OPChainBondLength::initialise() {
  for (const shared_ptr<Topology> &plugPtr : Sim->topology)
    if (std::dynamic_pointer_cast<TChain>(plugPtr))
      chains.push_back(
          Cdata(plugPtr->getID(), plugPtr->getMolecules().front()->size()));
}

void OPChainBondLength::ticker() {
  for (Cdata &dat : chains)
    for (const shared_ptr<IDRange> &range :
         Sim->topology[dat.chainID]->getMolecules())
      if (range->size() > 2)
        // Walk the polymer
        for (size_t j = 0; j < range->size() - 1; ++j)
          dat.BondLengths[j].addVal(
              (Sim->particles[(*range)[j + 1]].getPosition() -
               Sim->particles[(*range)[j]].getPosition())
                  .nrm());
}

void OPChainBondLength::output(magnet::xml::XmlStream &XML) {
  XML << magnet::xml::tag("BondAngleLength");

  for (Cdata &dat : chains) {
    XML << magnet::xml::tag("Chain") << magnet::xml::attr("Name")
        << Sim->topology[dat.chainID]->getName();

    size_t Nc = Sim->topology[dat.chainID]->getMolecules().front()->size() - 1;

    for (size_t i = 0; i < Nc; ++i)
      dat.BondLengths[i].outputHistogram(XML, 1.0 / Sim->units.unitLength());

    XML << magnet::xml::endtag("Chain");
  }

  XML << magnet::xml::endtag("BondAngleLength");
}
} // namespace dynamo
