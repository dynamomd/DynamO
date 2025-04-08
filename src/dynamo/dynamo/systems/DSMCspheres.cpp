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

#include <dynamo/systems/DSMCspheres.hpp>

#include <dynamo/BC/BC.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
SysDSMCSpheres::SysDSMCSpheres(const magnet::xml::Node &XML,
                               dynamo::Simulation *tmp)
    : System(tmp), maxprob(0.0) {
  dt = std::numeric_limits<float>::infinity();
  SysDSMCSpheres::operator<<(XML);
  type = DSMC;
}

SysDSMCSpheres::SysDSMCSpheres(dynamo::Simulation *nSim, double nd,
                               double ntstp, double nChi, double ne,
                               std::string nName, IDRange *r1, IDRange *r2)
    : System(nSim), tstep(ntstp), chi(nChi), d2(nd * nd), diameter(nd),
      maxprob(0.0), e(ne), range1(r1), range2(r2) {
  sysName = nName;
  type = DSMC;
}

NEventData SysDSMCSpheres::runEvent() {
  dt = tstep;
  std::normal_distribution<> norm_sampler;
  std::uniform_real_distribution<> uniform_sampler;
  std::uniform_int_distribution<size_t> id1sampler(0, range1->size() - 1);
  std::uniform_int_distribution<size_t> id2sampler(0, range2->size() - 1);

  // Find the likely maximum number of interacting pairs. The
  // addition of the random variable is a neat way to randomly pick
  // an extra pair to, on average, pick the correct number of
  // fractional pairs (thanks Severin!)
  const size_t nmax = static_cast<size_t>(0.5 * maxprob * range1->size() +
                                          uniform_sampler(Sim->ranGenerator));

  NEventData retval;

  for (size_t n = 0; n < nmax; ++n) {
    Particle &p1(
        Sim->particles[*(range1->begin() + id1sampler(Sim->ranGenerator))]);

    size_t p2id = *(range2->begin() + id2sampler(Sim->ranGenerator));

    // Find another particle which is not p1
    while (p2id == p1.getID())
      p2id = *(range2->begin() + id2sampler(Sim->ranGenerator));

    Particle &p2(Sim->particles[p2id]);

    Sim->dynamics->updateParticlePair(p1, p2);

    Vector rij;
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      rij[iDim] = norm_sampler(Sim->ranGenerator);

    // This is the extra diameter term missing from the "factor" variable
    rij *= diameter / rij.nrm();

    if (Sim->dynamics->DSMCSpheresTest(p1, p2, maxprob, factor, rij)) {
      ++Sim->eventCount;
      retval.L2partChanges.push_back(
          PairEventData(Sim->dynamics->DSMCSpheresRun(p1, p2, e, rij)));
    }
  }
  return retval;
}

void SysDSMCSpheres::initialise(size_t nID) {
  ID = nID;
  dt = tstep;

  // An extra factor of diameter is missing here, which is used to
  // give the vector rij below and in runEvent the "correct"
  // magnitude. This is incase applyBC in the dynamics (and any
  // other function) expects realistic rij values.
  factor = 4.0 * range2->size() * diameter * M_PI * chi * tstep /
           Sim->getSimVolume();

  if (maxprob == 0.0) {
    std::normal_distribution<> norm_sampler;
    std::uniform_int_distribution<size_t> id1sampler(0, range1->size() - 1);
    std::uniform_int_distribution<size_t> id2sampler(0, range2->size() - 1);

    // Just do some quick testing to get an estimate
    for (size_t n = 0; n < 1000; ++n) {
      Particle &p1(
          Sim->particles[*(range1->begin() + id1sampler(Sim->ranGenerator))]);

      size_t p2id = *(range2->begin() + id2sampler(Sim->ranGenerator));

      while (p2id == p1.getID())
        p2id = *(range2->begin() + id2sampler(Sim->ranGenerator));

      Particle &p2(Sim->particles[p2id]);

      Sim->dynamics->updateParticlePair(p1, p2);

      Vector rij;
      for (size_t iDim(0); iDim < NDIM; ++iDim)
        rij[iDim] = norm_sampler(Sim->ranGenerator);

      rij *= diameter / rij.nrm();

      Sim->dynamics->DSMCSpheresTest(p1, p2, maxprob, factor, rij);
    }
  }

  if (maxprob > 0.5)
    derr << "MaxProbability is " << maxprob << "\nNpairs per step is "
         << 0.5 * range1->size() * maxprob << std::endl;
  else
    dout << "MaxProbability is " << maxprob << "\nNpairs per step is "
         << 0.5 * range1->size() * maxprob << std::endl;

  if (0.5 * range1->size() * maxprob < 2.0)
    derr << "This probability is low" << std::endl;
}

void SysDSMCSpheres::operator<<(const magnet::xml::Node &XML) {
  tstep = XML.getAttribute("tStep").as<double>() * Sim->units.unitTime();
  chi = XML.getAttribute("Chi").as<double>();
  sysName = XML.getAttribute("Name");
  diameter =
      XML.getAttribute("Diameter").as<double>() * Sim->units.unitLength();
  e = XML.getAttribute("Inelasticity").as<double>();
  d2 = diameter * diameter;
  magnet::xml::Node subRangeXML = XML.getNode("IDRange");
  range1 = shared_ptr<IDRange>(IDRange::getClass(subRangeXML, Sim));
  ++subRangeXML;
  range2 = shared_ptr<IDRange>(IDRange::getClass(subRangeXML, Sim));
  if (XML.hasAttribute("MaxProbability"))
    maxprob = XML.getAttribute("MaxProbability").as<double>();
}

void SysDSMCSpheres::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::tag("System") << magnet::xml::attr("Type")
      << "DSMCSpheres" << magnet::xml::attr("tStep")
      << tstep / Sim->units.unitTime() << magnet::xml::attr("Chi") << chi
      << magnet::xml::attr("Diameter") << diameter / Sim->units.unitLength()
      << magnet::xml::attr("Inelasticity") << e << magnet::xml::attr("Name")
      << sysName << magnet::xml::attr("MaxProbability") << maxprob << range1
      << range2 << magnet::xml::endtag("System");
}
} // namespace dynamo
