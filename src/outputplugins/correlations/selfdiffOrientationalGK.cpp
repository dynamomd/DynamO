/*  DYNAMO:- Event driven molecular dynamics simulator
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../../dynamics/liouvillean/OrientationL.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"

OPSelfDiffusionOrientationalGK::OPSelfDiffusionOrientationalGK(const DYNAMO::SimData* tmp,const XMLNode& XML):
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

  if (dynamic_cast<const LNOrientation*>(&Sim->dynamics.getLiouvillean()) == NULL)
  {
    D_throw() << "Species does not specify an orientation";
  }

  G.resize(Sim->lN, boost::circular_buffer<VUpair> (CorrelatorLength, VUpair(Vector(0,0,0), Vector(0,0,0) )));

  accG2_parallel.resize(Sim->dynamics.getSpecies().size());
  accG2_perp.resize(Sim->dynamics.getSpecies().size());

  BOOST_FOREACH(std::vector<Iflt>& listref_parallel, accG2_parallel)
  {
    listref_parallel.resize(CorrelatorLength, 0);
  }

  BOOST_FOREACH(std::vector<Iflt>& listref_perp, accG2_perp)
  {
    listref_perp.resize(CorrelatorLength, 0);
  }

  I_cout() << "dt set to " << dt / Sim->dynamics.units().unitTime();
}

void
OPSelfDiffusionOrientationalGK::operator<<(const XMLNode& XML)
{
  try
  {
    if (XML.isAttributeSet("Length"))
    {
      CorrelatorLength = boost::lexical_cast<unsigned int> (XML.getAttribute("Length"));
    }

    if (XML.isAttributeSet("dt"))
    {
      dt = Sim->dynamics.units().unitTime() *
      boost::lexical_cast<Iflt>(XML.getAttribute("dt"));
    }

    if (XML.isAttributeSet("t"))
    {
      dt = Sim->dynamics.units().unitTime() * boost::lexical_cast<Iflt> (XML.getAttribute("t"))/CorrelatorLength;
    }
  }

  catch (boost::bad_lexical_cast &)
  {
    D_throw() << "Failed a lexical cast in OPSelfDiffusionOrientationalGK";
  }
}

void
OPSelfDiffusionOrientationalGK::eventUpdate(const CGlobEvent& iEvent, const CNParticleData& PDat)
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
OPSelfDiffusionOrientationalGK::eventUpdate(const CLocalEvent& iEvent, const CNParticleData& PDat)
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
OPSelfDiffusionOrientationalGK::eventUpdate(const CSystem&, const CNParticleData& PDat, const Iflt& edt)
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
OPSelfDiffusionOrientationalGK::eventUpdate(const CIntEvent& iEvent, const C2ParticleData& PDat)
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
OPSelfDiffusionOrientationalGK::newG(const C1ParticleData& PDat)
{
  if (Sim->dynamics.liouvilleanTypeTest<LSLLOD>())
  {
    Sim->dynamics.getLiouvillean().updateAllParticles();
  }

  for (size_t i = 0; i < Sim->lN; ++i)
  {
    const LNOrientation::rotData& rdat(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getRotData(Sim->vParticleList[i]));
    G[i].push_front(VUpair(Sim->vParticleList[i].getVelocity(), rdat.orientation));
  }

  //Now correct the fact that the wrong velocity and orientation have been pushed
  const LNOrientation::rotData& fetchpart(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getRotData(Sim->vParticleList[PDat.getParticle().getID()]));
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
OPSelfDiffusionOrientationalGK::newG(const C2ParticleData& PDat)
{
  for (size_t i = 0; i < Sim->lN; ++i)
  {
    const LNOrientation::rotData& rdat(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getRotData(Sim->vParticleList[i]));
    G[i].push_front(VUpair(Sim->vParticleList[i].getVelocity(), rdat.orientation));
  }

  //Now correct the fact that the wrong velocities and orientations have been pushed
  const LNOrientation::rotData& fetch1(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getRotData(Sim->vParticleList[PDat.particle1_.getParticle().getID()]));
  const LNOrientation::rotData& fetch2(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getRotData(Sim->vParticleList[PDat.particle2_.getParticle().getID()]));

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
OPSelfDiffusionOrientationalGK::newG(const CNParticleData& PDat)
{
  //This ensures the list stays at accumilator size
  for (size_t i = 0; i < Sim->lN; ++i)
  {
    const LNOrientation::rotData& rdat(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getRotData(Sim->vParticleList[i]));
    G[i].push_front(VUpair(Sim->vParticleList[i].getVelocity(), rdat.orientation));
  }

  //Go back and fix the pushes
  BOOST_FOREACH(const C1ParticleData&PDat2, PDat.L1partChanges)
  {
    const LNOrientation::rotData& fetchpart(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getRotData(Sim->vParticleList[PDat2.getParticle().getID()]));
    G[PDat2.getParticle().getID()].front() = VUpair(PDat2.getOldVel(), fetchpart.orientation);
  }

  BOOST_FOREACH(const C2ParticleData& PDat2, PDat.L2partChanges)
  {
    const LNOrientation::rotData& fetch1(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getRotData(Sim->vParticleList[PDat2.particle1_.getParticle().getID()]));
    const LNOrientation::rotData& fetch2(static_cast<const LNOrientation&> (Sim->dynamics.getLiouvillean()).getRotData(Sim->vParticleList[PDat2.particle2_.getParticle().getID()]));

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
OPSelfDiffusionOrientationalGK::output(xmlw::XmlStream& XML)
{
  Iflt factor = Sim->dynamics.units().unitTime() / (Sim->dynamics.units().unitDiffusion() * count);

  // Counting over the PERPENDICULAR vector - both should be same size anyway
  for (size_t i = 0; i < accG2_perp.size(); ++i)
  {
    // Common headers
    Iflt specCount = Sim->dynamics.getSpecies()[i]->getCount();

    XML << xmlw::tag("Correlator")
	<< xmlw::attr("name") << "SelfDiffusionOrientationalGK"
	<< xmlw::attr("species") << Sim->dynamics.getSpecies()[i]->getName()
	<< xmlw::attr("size") << accG2_perp.size()
	<< xmlw::attr("dt") << dt / Sim->dynamics.units().unitTime()
	<< xmlw::attr("LengthInMFT") << dt * accG2_perp[i].size() / Sim->getOutputPlugin<OPMisc>()->getMFT()
	<< xmlw::attr("simFactor") << factor / specCount
	<< xmlw::attr("SampleCount") << count;


    // Perpendicular section
    Iflt acc_perp = 0.5*(accG2_perp[i].front() + accG2_perp[i].back());

    for (size_t j = 1; j < accG2_perp[i].size() - 1; ++j)
      { acc_perp += accG2_perp[i][j];}

    acc_perp *= factor * dt / (Sim->dynamics.units().unitTime() * specCount);

    XML << xmlw::tag("Component")
	<< xmlw::attr("Type") << "Perpendicular"
	<< xmlw::tag("Integral")
	<< xmlw::attr("value") << acc_perp
	<< xmlw::endtag("Integral")
	<< xmlw::chardata();

    for (size_t j = 0; j < accG2_perp[i].size(); ++j)
    {
      XML << j * dt / Sim->dynamics.units().unitTime()
	  << "\t" << accG2_perp[i][j] * factor / specCount << "\n";
    }

    XML << xmlw::endtag("Component");

    // Parallel section
    Iflt acc_parallel = 0.5*(accG2_parallel[i].front() + accG2_parallel[i].back());

    for (size_t j = 1; j < accG2_parallel[i].size() - 1; ++j)
    {
      acc_parallel += accG2_parallel[i][j];
    }

    acc_parallel *= factor * dt / (Sim->dynamics.units().unitTime() * specCount);

    XML << xmlw::tag("Component")
	<< xmlw::attr("Type") << "Parallel"
	<< xmlw::tag("Integral")
	<< xmlw::attr("value") << acc_parallel
	<< xmlw::endtag("Integral")
	<< xmlw::chardata();

    for (size_t j = 0; j < accG2_parallel[i].size(); ++j)
    {
      XML << j * dt / Sim->dynamics.units().unitTime()
	  << "\t" << accG2_parallel[i][j] * factor / specCount << "\n";
    }

    XML << xmlw::endtag("Component")
	<< xmlw::endtag("Correlator");
    }
}

Iflt
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
      return 10.0 / (((Iflt) CorrelatorLength)*sqrt(Sim->dynamics.getLiouvillean().getkT()) * CorrelatorLength);
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

  BOOST_FOREACH(const smrtPlugPtr<CSpecies>& spec, Sim->dynamics.getSpecies())
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
