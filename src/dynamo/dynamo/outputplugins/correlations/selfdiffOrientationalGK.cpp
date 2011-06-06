/*  dynamo:- Event driven molecular dynamics simulator
    http://www.marcusbannerman.co.uk/dynamo
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

#include "selfdiffOrientationalGK.hpp"
#include "../../dynamics/liouvillean/SLLOD.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include <magnet/math/matrix.hpp>

OPSelfDiffusionOrientationalGK::OPSelfDiffusionOrientationalGK(const dynamo::SimData* tmp,
							       const magnet::xml::Node& XML):
  OutputPlugin(tmp, "SelfDiffusionOrientationalGK", 60), //Note the sort order set later
  count(0),
  dt(0),
  currentdt(0.0),
  CorrelatorLength(100),
  currCorrLen(0),
  notReady(true)
{
  operator<<(XML);
}

void
OPSelfDiffusionOrientationalGK::initialise()
{
  Sim->getOutputPlugin<OPMisc>();

  dt = getdt();

  if (!(Sim->dynamics.getLiouvillean().hasOrientationData()))
    M_throw() << "There is no orientation data available.";

  G.resize(Sim->N, boost::circular_buffer<VUpair> (CorrelatorLength, VUpair(Vector(0,0,0), Vector(0,0,0) )));

  accG2_parallel.resize(Sim->dynamics.getSpecies().size());
  accG2_perp.resize(Sim->dynamics.getSpecies().size());

  BOOST_FOREACH(std::vector<double>& listref_parallel, accG2_parallel)
  {
    listref_parallel.resize(CorrelatorLength, 0);
  }

  BOOST_FOREACH(std::vector<double>& listref_perp, accG2_perp)
  {
    listref_perp.resize(CorrelatorLength, 0);
  }

  I_cout() << "dt set to " << dt / Sim->dynamics.units().unitTime();
}

void
OPSelfDiffusionOrientationalGK::operator<<(const magnet::xml::Node& XML)
{
  try
  {
    CorrelatorLength = XML.getAttribute("Length").as<size_t>(100);
    if (XML.getAttribute("dt").valid())
      dt = XML.getAttribute("dt").as<double>() * Sim->dynamics.units().unitTime();

    if (XML.getAttribute("t").valid())
      dt = XML.getAttribute("t").as<double>() * Sim->dynamics.units().unitTime() 
	/ CorrelatorLength;
  }
  catch (boost::bad_lexical_cast &)
  {
    M_throw() << "Failed a lexical cast in OPSelfDiffusionOrientationalGK";
  }
}

void
OPSelfDiffusionOrientationalGK::eventUpdate(const GlobalEvent& iEvent, const NEventData& PDat)
{
  //Move the time forward
  currentdt += iEvent.getdt();

  //Now test if we've gone over the step time
  while (currentdt >= dt)
    {
      currentdt -= dt;
      newG(PDat);
    }
}

void
OPSelfDiffusionOrientationalGK::eventUpdate(const LocalEvent& iEvent, const NEventData& PDat)
{
  //Move the time forward
  currentdt += iEvent.getdt();

  //Now test if we've gone over the step time
  while (currentdt >= dt)
    {
      currentdt -= dt;
      newG(PDat);
    }
}

void
OPSelfDiffusionOrientationalGK::eventUpdate(const System&, const NEventData& PDat, const double& edt)
{
  //Move the time forward
  currentdt += edt;

  //Now test if we've gone over the step time
  while (currentdt >= dt)
    {
      currentdt -= dt;
      newG(PDat);
    }
}

void
OPSelfDiffusionOrientationalGK::eventUpdate(const IntEvent& iEvent, const PairEventData& PDat)
{
  //Move the time forward
  currentdt += iEvent.getdt();

  //Now test if we've gone over the step time
  while (currentdt >= dt)
    {
      currentdt -= dt;
      newG(PDat);
    }
}

void
OPSelfDiffusionOrientationalGK::newG(const ParticleEventData& PDat)
{
  if (Sim->dynamics.liouvilleanTypeTest<LSLLOD>())
  {
    Sim->dynamics.getLiouvillean().updateAllParticles();
  }

  for (size_t i = 0; i < Sim->N; ++i)
  {
    const Liouvillean::rotData& rdat(Sim->dynamics.getLiouvillean().getRotData(Sim->particleList[i]));
    G[i].push_front(VUpair(Sim->particleList[i].getVelocity(), rdat.orientation));
  }

  //Now correct the fact that the wrong velocity and orientation have been pushed
  const Liouvillean::rotData& fetchpart(Sim->dynamics.getLiouvillean()
					.getRotData(Sim->particleList[PDat.getParticle().getID()]));
  G[PDat.getParticle().getID()].front() = VUpair(PDat.getOldVel(), fetchpart.orientation);

  //This ensures the list gets to accumilator size
  if (notReady)
  {
    if (++currCorrLen != CorrelatorLength)
    {
      return;
    }

      notReady = false;
  }

  accPass();
}

void
OPSelfDiffusionOrientationalGK::newG(const PairEventData& PDat)
{
  for (size_t i = 0; i < Sim->N; ++i)
  {
    const Liouvillean::rotData& rdat(Sim->dynamics.getLiouvillean().getRotData(Sim->particleList[i]));
    G[i].push_front(VUpair(Sim->particleList[i].getVelocity(), rdat.orientation));
  }

  //Now correct the fact that the wrong velocities and orientations have been pushed
  const Liouvillean::rotData& fetch1(Sim->dynamics.getLiouvillean()
				     .getRotData(Sim->particleList[PDat.particle1_.getParticle().getID()]));
  const Liouvillean::rotData& fetch2(Sim->dynamics.getLiouvillean()
				     .getRotData(Sim->particleList[PDat.particle2_.getParticle().getID()]));

  G[PDat.particle1_.getParticle().getID()].front() = VUpair(PDat.particle1_.getOldVel(), fetch1.orientation);
  G[PDat.particle2_.getParticle().getID()].front() = VUpair(PDat.particle2_.getOldVel(), fetch2.orientation);

  //This ensures the list gets to accumilator size
  if (notReady)
  {
    if (++currCorrLen != CorrelatorLength)
    {
      return;
    }

    notReady = false;
  }

  accPass();
}

void
OPSelfDiffusionOrientationalGK::newG(const NEventData& PDat)
{
  //This ensures the list stays at accumilator size
  for (size_t i = 0; i < Sim->N; ++i)
  {
    const Liouvillean::rotData& rdat(Sim->dynamics.getLiouvillean().getRotData(Sim->particleList[i]));
    G[i].push_front(VUpair(Sim->particleList[i].getVelocity(), rdat.orientation));
  }

  //Go back and fix the pushes
  BOOST_FOREACH(const ParticleEventData&PDat2, PDat.L1partChanges)
  {
    const Liouvillean::rotData& fetchpart(Sim->dynamics.getLiouvillean()
					  .getRotData(Sim->particleList[PDat2.getParticle().getID()]));
    G[PDat2.getParticle().getID()].front() = VUpair(PDat2.getOldVel(), fetchpart.orientation);
  }

  BOOST_FOREACH(const PairEventData& PDat2, PDat.L2partChanges)
  {
    const Liouvillean::rotData& fetch1(Sim->dynamics.getLiouvillean().getRotData(Sim->particleList[PDat2.particle1_.getParticle().getID()]));
    const Liouvillean::rotData& fetch2(Sim->dynamics.getLiouvillean().getRotData(Sim->particleList[PDat2.particle2_.getParticle().getID()]));

    G[PDat2.particle1_.getParticle().getID()].front() = VUpair(PDat2.particle1_.getOldVel(), fetch1.orientation);
    G[PDat2.particle2_.getParticle().getID()].front() = VUpair(PDat2.particle2_.getOldVel(), fetch2.orientation);
  }

  //This ensures the list gets to accumilator size
  if (notReady)
  {
    if (++currCorrLen != CorrelatorLength)
    {
      return;
    }

    notReady = false;
  }

  accPass();
}

void
OPSelfDiffusionOrientationalGK::output(xml::XmlStream& XML)
{
  double factor = Sim->dynamics.units().unitTime() / (Sim->dynamics.units().unitDiffusion() * count);

  // Counting over the PERPENDICULAR vector - both should be same size anyway
  for (size_t i = 0; i < accG2_perp.size(); ++i)
  {
    // Common headers
    double specCount = Sim->dynamics.getSpecies()[i]->getCount();

    XML << xml::tag("Correlator")
	<< xml::attr("name") << "SelfDiffusionOrientationalGK"
	<< xml::attr("species") << Sim->dynamics.getSpecies()[i]->getName()
	<< xml::attr("size") << accG2_perp.size()
	<< xml::attr("dt") << dt / Sim->dynamics.units().unitTime()
	<< xml::attr("LengthInMFT") << dt * accG2_perp[i].size() / Sim->getOutputPlugin<OPMisc>()->getMFT()
	<< xml::attr("simFactor") << factor / specCount
	<< xml::attr("SampleCount") << count;


    // Perpendicular section
    double acc_perp = 0.5*(accG2_perp[i].front() + accG2_perp[i].back());

    for (size_t j = 1; j < accG2_perp[i].size() - 1; ++j)
      { acc_perp += accG2_perp[i][j];}

    acc_perp *= factor * dt / (Sim->dynamics.units().unitTime() * specCount);

    XML << xml::tag("Component")
	<< xml::attr("Type") << "Perpendicular"
	<< xml::tag("Integral")
	<< xml::attr("value") << acc_perp
	<< xml::endtag("Integral")
	<< xml::chardata();

    for (size_t j = 0; j < accG2_perp[i].size(); ++j)
    {
      XML << j * dt / Sim->dynamics.units().unitTime()
	  << "\t" << accG2_perp[i][j] * factor / specCount << "\n";
    }

    XML << xml::endtag("Component");

    // Parallel section
    double acc_parallel = 0.5*(accG2_parallel[i].front() + accG2_parallel[i].back());

    for (size_t j = 1; j < accG2_parallel[i].size() - 1; ++j)
    {
      acc_parallel += accG2_parallel[i][j];
    }

    acc_parallel *= factor * dt / (Sim->dynamics.units().unitTime() * specCount);

    XML << xml::tag("Component")
	<< xml::attr("Type") << "Parallel"
	<< xml::tag("Integral")
	<< xml::attr("value") << acc_parallel
	<< xml::endtag("Integral")
	<< xml::chardata();

    for (size_t j = 0; j < accG2_parallel[i].size(); ++j)
    {
      XML << j * dt / Sim->dynamics.units().unitTime()
	  << "\t" << accG2_parallel[i][j] * factor / specCount << "\n";
    }

    XML << xml::endtag("Component")
	<< xml::endtag("Correlator");
    }
}

double
OPSelfDiffusionOrientationalGK::getdt()
{
  //Get the simulation temperature
  if (dt == 0.0)
  {
    if (Sim->lastRunMFT != 0.0)
    {
      return Sim->lastRunMFT * 50.0 / CorrelatorLength;
    }
    else
    {
      return 10.0 / (((double) CorrelatorLength)*sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
    }
  }
  else
  {
    return dt;
  }
}

void
OPSelfDiffusionOrientationalGK::accPass()
{
  ++count;

  BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
  {
    BOOST_FOREACH(const size_t& ID, *spec->getRange())
    {
      for (size_t j = 0; j < CorrelatorLength; ++j)
      {
	// Parallel = <[v(t).u(0)][v(0).u(0)]>
	accG2_parallel[spec->getID()][j] +=  ((G[ID].front().first |  G[ID][j].second) * (G[ID][j].first | G[ID][j].second));

	// Get a unit matrix
	Matrix unit, dyad;
	unit.one();
	dyad = Dyadic(G[ID][j].second, G[ID][j].second);

	// Perpendicular = <v(t).[I - u(0)u(0)]v(0)>
	accG2_perp[spec->getID()][j] += (G[ID].front().first | ((unit - dyad) * G[ID][j].first));

	//for (size_t iDim(0); iDim < NDIM; ++iDim)
	//{
	//}
      }
    }
  }
}
