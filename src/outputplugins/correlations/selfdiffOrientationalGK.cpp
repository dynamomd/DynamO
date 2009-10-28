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

  G_parallel.resize(Sim->lN, boost::circular_buffer<Vector  >(CorrelatorLength, Vector(0,0,0)));
  G_perp.resize(Sim->lN, boost::circular_buffer<Vector  >(CorrelatorLength, Vector(0,0,0)));

  accG2_parallel.resize(Sim->dynamics.getSpecies().size());
  accG2_perp.resize(Sim->dynamics.getSpecies().size());

  BOOST_FOREACH(std::vector<Vector  >& listref_parallel, accG2_parallel)
  {
    listref_parallel.resize(CorrelatorLength, Vector (0,0,0));
  }

  BOOST_FOREACH(std::vector<Vector  >& listref_perp, accG2_perp)
  {
    listref_perp.resize(CorrelatorLength, Vector (0,0,0));
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
			CorrelatorLength = boost::lexical_cast<unsigned int>
				(XML.getAttribute("Length"));
		}

		if (XML.isAttributeSet("dt"))
		{
			dt = Sim->dynamics.units().unitTime() *
				boost::lexical_cast<Iflt>(XML.getAttribute("dt"));
		}

		if (XML.isAttributeSet("t"))
		{
			dt = Sim->dynamics.units().unitTime() *
				boost::lexical_cast<Iflt>
				(XML.getAttribute("t"))/CorrelatorLength;
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
	if (Sim->dynamics.liouvilleanTypeTest<CLSLLOD>())
	{
		Sim->dynamics.getLiouvillean().updateAllParticles();
	}

	for (size_t i = 0; i < Sim->lN; ++i)
	{
		G_perp[i].push_front(Sim->vParticleList[i].getVelocity());
		G_parallel[i].push_front(Sim->vParticleList[i].getVelocity());
	}

	//Now correct the fact that the wrong velocity has been pushed
	G_perp[PDat.getParticle().getID()].front() = PDat.getOldVel();
	G_parallel[PDat.getParticle().getID()].front() = PDat.getOldVel();

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
		G_perp[i].push_front(Sim->vParticleList[i].getVelocity());
		G_parallel[i].push_front(Sim->vParticleList[i].getVelocity());
	}

	//Now correct the fact that the wrong velocity has been pushed
	G_perp[PDat.particle1_.getParticle().getID()].front() = PDat.particle1_.getOldVel();
	G_perp[PDat.particle2_.getParticle().getID()].front() = PDat.particle2_.getOldVel();

	G_parallel[PDat.particle1_.getParticle().getID()].front() = PDat.particle1_.getOldVel();
	G_parallel[PDat.particle2_.getParticle().getID()].front() = PDat.particle2_.getOldVel();

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
		G_perp[i].push_front(Sim->vParticleList[i].getVelocity());
		G_parallel[i].push_front(Sim->vParticleList[i].getVelocity());
	}

	//Go back and fix the pushes
	BOOST_FOREACH(const C1ParticleData&PDat2, PDat.L1partChanges)
	{
		G_perp[PDat2.getParticle().getID()].front() = PDat2.getOldVel();
		G_parallel[PDat2.getParticle().getID()].front() = PDat2.getOldVel();
	}

	BOOST_FOREACH(const C2ParticleData& PDat2, PDat.L2partChanges)
    {
		G_perp[PDat2.particle1_.getParticle().getID()].front() = PDat2.particle1_.getOldVel();
		G_perp[PDat2.particle2_.getParticle().getID()].front() = PDat2.particle2_.getOldVel();

		G_parallel[PDat2.particle1_.getParticle().getID()].front() = PDat2.particle1_.getOldVel();
		G_parallel[PDat2.particle2_.getParticle().getID()].front() = PDat2.particle2_.getOldVel();
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
	Iflt factor = Sim->dynamics.units().unitTime()
		/ (Sim->dynamics.units().unitDiffusion() * count);

	// Run through the perpendicular components
	for (size_t i = 0; i < accG2_perp.size(); ++i)
    {
		Iflt specCount = Sim->dynamics.getSpecies()[i]->getCount();

		Vector acc_perp = 0.5*(accG2_perp[i].front() + accG2_perp[i].back());

		for (size_t j = 1; j < accG2_perp[i].size() - 1; ++j)
		{
			acc_perp += accG2_perp[i][j];
		}

		acc_perp *= factor * dt / (Sim->dynamics.units().unitTime() * specCount);

		XML << xmlw::tag("Correlator")
			<< xmlw::attr("name") << "SelfDiffusionPerpendicularGK"
			<< xmlw::attr("species") << Sim->dynamics.getSpecies()[i]->getName()
			<< xmlw::attr("size") << accG2_perp.size()
			<< xmlw::attr("dt") << dt / Sim->dynamics.units().unitTime()
			<< xmlw::attr("LengthInMFT") << dt * accG2_perp[i].size()
				/ Sim->getOutputPlugin<OPMisc>()->getMFT()
			<< xmlw::attr("simFactor") << factor / specCount
			<< xmlw::attr("SampleCount") << count
			<< xmlw::tag("Integral") << acc_perp
			<< xmlw::endtag("Integral")
			<< xmlw::chardata();

		for (size_t j = 0; j < accG2_perp[i].size(); ++j)
		{
			XML << j * dt / Sim->dynamics.units().unitTime();
			for (size_t iDim = 0; iDim < NDIM; iDim++)
				XML << "\t" << accG2_perp[i][j][iDim] * factor / specCount;
			XML << "\n";
		}

		XML << xmlw::endtag("Correlator");
    }

	// Run through the parallel components
	for (size_t i = 0; i < accG2_parallel.size(); ++i)
    {
		Iflt specCount = Sim->dynamics.getSpecies()[i]->getCount();

		Vector acc_parallel = 0.5*(accG2_parallel[i].front() + accG2_parallel[i].back());

		for (size_t j = 1; j < accG2_parallel[i].size() - 1; ++j)
		{
			acc_parallel += accG2_parallel[i][j];
		}

		acc_parallel *= factor * dt / (Sim->dynamics.units().unitTime() * specCount);

		XML << xmlw::tag("Correlator")
			<< xmlw::attr("name") << "SelfDiffusionParallelGK"
			<< xmlw::attr("species") << Sim->dynamics.getSpecies()[i]->getName()
			<< xmlw::attr("size") << accG2_parallel.size()
			<< xmlw::attr("dt") << dt / Sim->dynamics.units().unitTime()
			<< xmlw::attr("LengthInMFT") << dt * accG2_parallel[i].size()
				/ Sim->getOutputPlugin<OPMisc>()->getMFT()
			<< xmlw::attr("simFactor") << factor / specCount
			<< xmlw::attr("SampleCount") << count
			<< xmlw::tag("Integral") << acc_parallel
			<< xmlw::endtag("Integral")
			<< xmlw::chardata();

		for (size_t j = 0; j < accG2_parallel[i].size(); ++j)
		{
			XML << j * dt / Sim->dynamics.units().unitTime();
			for (size_t iDim = 0; iDim < NDIM; iDim++)
				XML << "\t" << accG2_parallel[i][j][iDim] * factor / specCount;
			XML << "\n";
		}

		XML << xmlw::endtag("Correlator");
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
				for (size_t iDim(0); iDim < NDIM; ++iDim)
				{
					accG2_parallel[spec->getID()][j][iDim] +=  G_parallel[ID].front()[iDim] * G_parallel[ID][j][iDim];
					accG2_perp[spec->getID()][j][iDim] +=  G_perp[ID].front()[iDim] * G_perp[ID][j][iDim];
				}
			}
		}
	}
}
