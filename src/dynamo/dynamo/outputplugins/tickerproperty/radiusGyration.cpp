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

#include <cmath>
#include <dynamo/include.hpp>
#include <dynamo/outputplugins/tickerproperty/radiusGyration.hpp>
#include <dynamo/ranges/IDRange.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/topology/include.hpp>
#include <magnet/math/matrix.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>
#include <vector>

namespace dynamo {
OPRGyration::OPRGyration(const dynamo::Simulation *tmp,
                         const magnet::xml::Node &XML)
    : OPTicker(tmp, "GyrationRadius"), _binWidthGyration(0.01),
      _binWidthNematic(0.001) {
  operator<<(XML);
}

void OPRGyration::operator<<(const magnet::xml::Node &XML) {
  if (XML.hasAttribute("BinWidthGyration"))
    _binWidthGyration = XML.getAttribute("BinWidthGyration").as<double>();

  if (XML.hasAttribute("BinWidthGyration"))
    _binWidthNematic = XML.getAttribute("BinWidthGyration").as<double>();
}

void OPRGyration::initialise() {
  for (const shared_ptr<Topology> &plugPtr : Sim->topology)
    if (std::dynamic_pointer_cast<TChain>(plugPtr))
      chains.push_back(CTCdata(static_cast<const TChain *>(plugPtr.get()),
                               _binWidthGyration * Sim->units.unitArea(),
                               _binWidthNematic));
}

void OPRGyration::replicaExchange(OutputPlugin &plug) {
  std::list<CTCdata>::iterator iPtr1 = chains.begin(),
                               iPtr2 = static_cast<OPRGyration &>(plug)
                                           .chains.begin();

#ifdef DYNAMO_DEBUG
  if (chains.size() != static_cast<OPRGyration &>(plug).chains.size())
    M_throw() << "Size mismatch when exchanging!";
#endif

  while (iPtr1 != chains.end()) {
#ifdef DYNAMO_DEBUG
    if (iPtr1->chainPtr->getName() != iPtr2->chainPtr->getName())
      M_throw() << "Name mismatch while replexing!";
#endif
    std::swap(iPtr1->gyrationRadii, iPtr2->gyrationRadii);
    std::swap(iPtr1->nematicOrder, iPtr2->nematicOrder);
    ++iPtr1;
    ++iPtr2;
  }

  std::swap(_binWidthGyration, _binWidthGyration);
  std::swap(_binWidthNematic, _binWidthNematic);
}

OPRGyration::molGyrationDat
OPRGyration::getGyrationEigenSystem(const shared_ptr<IDRange> &range,
                                    const dynamo::Simulation *Sim) {
  // Determine the centre of mass. Watch for periodic images
  Vector tmpVec;

  molGyrationDat retVal;
  retVal.MassCentre = Vector{0, 0, 0};

  double totmass = Sim->species[Sim->particles[*(range->begin())]]->getMass(
      *(range->begin()));
  // Walk along the chain
  Vector origin_position = Vector{0, 0, 0};
  Matrix inertiaTensor;

  for (IDRange::iterator iPtr = range->begin() + 1; iPtr != range->end();
       iPtr++) {
    Vector currRelPos = Sim->particles[*iPtr].getPosition() -
                        Sim->particles[*(iPtr - 1)].getPosition();
    Sim->BCs->applyBC(currRelPos);

    const Vector unfolded_pos = currRelPos + origin_position;

    const double mass = Sim->species[Sim->particles[*iPtr]]->getMass(*iPtr);

    retVal.MassCentre += origin_position * mass;
    inertiaTensor +=
        mass * ((unfolded_pos * unfolded_pos) * Matrix::identity() -
                Dyadic(unfolded_pos, unfolded_pos));
    totmass += mass;
    origin_position = unfolded_pos;
  }

  retVal.MassCentre /= totmass;
  retVal.MassCentre += Sim->particles[*(range->begin())].getPosition();

  std::pair<std::array<Vector, 3>, std::array<double, 3>> result =
      magnet::math::symmetric_eigen_decomposition(inertiaTensor / totmass);

  for (size_t i = 0; i < NDIM; i++) {
    retVal.EigenVal[i] = result.second[i] / range->size();

    // EigenVec Components
    for (size_t j = 0; j < NDIM; j++)
      retVal.EigenVec[i][j] = result.first[i][j];
  }
  return retVal;
}

Vector OPRGyration::NematicOrderParameter(const std::list<Vector> &molAxis) {
  Matrix Q;

  for (const Vector &vec : molAxis)
    for (size_t i = 0; i < NDIM; i++)
      for (size_t j = i; j < NDIM; j++)
        Q(i, j) += (3.0 * vec[i] * vec[j]) - (i == j ? 1 : 0);

  double Factor = 1.0 / (2.0 * molAxis.size());

  for (size_t i = 0; i < NDIM; i++)
    for (size_t j = i; j < NDIM; j++)
      Q(i, j) *= Factor;

  // Copy over the triangle matrix
  for (size_t i = 0; i < NDIM - 1; i++)
    for (size_t j = i + 1; j < NDIM; j++)
      Q(j, i) = Q(i, j);

  std::pair<std::array<Vector, 3>, std::array<double, 3>> result =
      magnet::math::symmetric_eigen_decomposition(Q);

  return Vector{result.second[0], result.second[1], result.second[2]};
}

void OPRGyration::ticker() {
  for (CTCdata &dat : chains) {
    std::list<Vector> molAxis;

    for (const shared_ptr<IDRange> &range : dat.chainPtr->getMolecules()) {
      molGyrationDat vals = getGyrationEigenSystem(range, Sim);
      // Take the largest eigenvector as the molecular axis
      molAxis.push_back(vals.EigenVec[NDIM - 1]);
      // Now add the radius of gyration
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
        dat.gyrationRadii[iDim].addVal(vals.EigenVal[iDim]);
    }

    Vector EigenVal = NematicOrderParameter(molAxis);

    for (size_t i = 0; i < NDIM; i++)
      if (std::isnormal(EigenVal[i]))
        dat.nematicOrder[i].addVal(EigenVal[i]);
  }
}

void OPRGyration::output(magnet::xml::XmlStream &XML) {
  XML << magnet::xml::tag("ChainGyration");

  for (CTCdata &dat : chains) {
    XML << magnet::xml::tag("Chain") << magnet::xml::attr("Name")
        << dat.chainPtr->getName().c_str() << magnet::xml::tag("GyrationRadii");

    for (size_t i = 0; i < NDIM; i++)
      dat.gyrationRadii.at(i).outputHistogram(XML, 1.0 / Sim->units.unitArea());

    XML << magnet::xml::endtag("GyrationRadii")
        << magnet::xml::tag("NematicOrderParameter");

    std::list<Vector> molAxis;

    for (const shared_ptr<IDRange> &range : dat.chainPtr->getMolecules())
      molAxis.push_back(getGyrationEigenSystem(range, Sim).EigenVec[NDIM - 1]);

    Vector EigenVal = NematicOrderParameter(molAxis);

    for (size_t i = 0; i < NDIM; i++)
      if (std::isnormal(EigenVal[i])) {
        char lett[2] = {char('x' + i), '\0'};

        XML << magnet::xml::attr(lett) << EigenVal[i];
      }

    for (size_t i = 0; i < NDIM; i++)
      dat.nematicOrder.at(i).outputHistogram(XML, 1.0);

    XML << magnet::xml::endtag("NematicOrderParameter")
        << magnet::xml::endtag("Chain");
  }

  XML << magnet::xml::endtag("ChainGyration");
}
} // namespace dynamo
