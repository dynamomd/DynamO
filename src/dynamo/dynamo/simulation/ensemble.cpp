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

#include <dynamo/simulation/ensemble.hpp>
#include <dynamo/dynamics/systems/ghost.hpp>
#include <dynamo/dynamics/liouvillean/CompressionL.hpp>
#include <dynamo/dynamics/BC/LEBC.hpp>
#include <dynamo/outputplugins/1partproperty/uenergy.hpp>
#include <dynamo/dynamics/liouvillean/NewtonMCL.hpp>

#include <magnet/exception.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  Ensemble* 
  Ensemble::getClass(const magnet::xml::Node& XML, const dynamo::SimData* Sim)
  {
    if (!strcmp(XML.getAttribute("Type"), "NVT"))
      return new EnsembleNVT(Sim);
    else if (!strcmp(XML.getAttribute("Type"), "NVE"))
      return new EnsembleNVE(Sim);
    else if (!strcmp(XML.getAttribute("Type"), "NVShear"))
      return new EnsembleNVShear(Sim);
    else if (!strcmp(XML.getAttribute("Type"), "NECompression"))
      return new EnsembleNECompression(Sim);
    else if (!strcmp(XML.getAttribute("Type"), "NTCompression"))
      return new EnsembleNTCompression(Sim);
    else
      M_throw() << "Cannot correctly identify the ensemble";
  }

  double 
  Ensemble::exchangeProbability(const Ensemble&) const
  { M_throw() << "Exchange move not written for this Ensemble"; }

  xml::XmlStream& operator<<(xml::XmlStream& XML, 
			      const Ensemble& g)
  {
    XML << xml::tag("Ensemble")
	<< xml::attr("Type") << g.getName()
	<< xml::endtag("Ensemble");
    return XML;
  }

  void
  EnsembleNVE::initialise()
  {
    EnsembleVals[0] = Sim->particleList.size();
    EnsembleVals[1] = Sim->primaryCellSize[0] * Sim->primaryCellSize[1] * Sim->primaryCellSize[2];
    EnsembleVals[2] = Sim->dynamics.calcInternalEnergy() + Sim->dynamics.getLiouvillean().getSystemKineticEnergy();

    dout << "NVE Ensemble initialised\nN=" << EnsembleVals[0]
	     << "\nV=" << EnsembleVals[1] / Sim->dynamics.units().unitVolume()
	     << "\nE=" << EnsembleVals[2] / Sim->dynamics.units().unitEnergy() << std::endl;
  }

  boost::array<double,3> 
  EnsembleNVE::getReducedEnsembleVals() const
  {
    boost::array<double,3> retval;
    retval[0] = EnsembleVals[0];
    retval[1] = EnsembleVals[1] / Sim->dynamics.units().unitVolume();
    retval[2] = EnsembleVals[2] / Sim->dynamics.units().unitEnergy();

    return retval;
  }

  void
  EnsembleNVT::initialise()
  {
    EnsembleVals[0] = Sim->particleList.size();
    EnsembleVals[1] = Sim->dynamics.units().unitVolume();

    try {
      thermostat = Sim->dynamics.getSystem("Thermostat").get_ptr();
    } catch (std::exception &)
      {
	M_throw() << "Could not find the Thermostat in NVT system";
      }
    
    //Only one kind of thermostat so far!
    if (dynamic_cast<const CSysGhost*>(thermostat) == NULL)
      {
	M_throw() << "Could not upcast thermostat to Andersens";
      }    
    
    EnsembleVals[2] = static_cast<const CSysGhost*>(thermostat)->getTemperature();
    
    dout << "NVT Ensemble initialised\nN=" << EnsembleVals[0]
	     << "\nV=" << EnsembleVals[1] / Sim->dynamics.units().unitVolume()
	     << "\nT=" << EnsembleVals[2] / Sim->dynamics.units().unitEnergy() << std::endl;
  }

  boost::array<double,3> 
  EnsembleNVT::getReducedEnsembleVals() const
  {
    boost::array<double,3> retval;
    retval[0] = EnsembleVals[0];
    retval[1] = EnsembleVals[1] / Sim->dynamics.units().unitVolume();
    retval[2] = EnsembleVals[2] / Sim->dynamics.units().unitEnergy();

    return retval;
  }

  double 
  EnsembleNVT::exchangeProbability(const Ensemble& oE) const
  {
#ifdef DYNAMO_DEBUG
    try {
      dynamic_cast<const EnsembleNVT&>(oE);
    } catch (std::bad_cast)
      {
	M_throw() << "The ensembles types differ";
      }    
#endif

    //Must use static cast to allow access to protected members
    //This is -\Delta in the Sugita_Okamoto paper
    
    const EnsembleNVT& ensemble2(static_cast<const EnsembleNVT&>(oE));

    double beta1 = 1 / EnsembleVals[2];
    double E1 = Sim->getOutputPlugin<OPUEnergy>()->getSimU();
    double beta2 = 1 / ensemble2.getEnsembleVals()[2];
    double E2 = ensemble2.Sim->getOutputPlugin<OPUEnergy>()->getSimU();
    
    double factor = (E1 - E2) * (beta1 - beta2);
    
    if (Sim->dynamics.liouvilleanTypeTest<LNewtonianMC>())
      {
	factor -= static_cast<const LNewtonianMC&>(Sim->dynamics.getLiouvillean()).W(E1);
	factor += static_cast<const LNewtonianMC&>(Sim->dynamics.getLiouvillean()).W(E2);
      }

    if (ensemble2.Sim->dynamics.liouvilleanTypeTest<LNewtonianMC>())
      {
	factor -= static_cast<const LNewtonianMC&>(ensemble2.Sim->dynamics.getLiouvillean()).W(E2);
	factor += static_cast<const LNewtonianMC&>(ensemble2.Sim->dynamics.getLiouvillean()).W(E1);
      }

    return std::exp(factor);
  }

  void
  EnsembleNVShear::initialise()
  {
    EnsembleVals[0] = Sim->particleList.size();
    EnsembleVals[1] = Sim->primaryCellSize[0] * Sim->primaryCellSize[1] * Sim->primaryCellSize[2];
    EnsembleVals[2] = CLEBC::shearRate();

    dout << "NVShear Ensemble initialised\nN=" << EnsembleVals[0]
	     << "\nV=" << EnsembleVals[1] / Sim->dynamics.units().unitVolume()
	     << "\nGamma=" << EnsembleVals[2] * Sim->dynamics.units().unitTime() << std::endl;
  }

  boost::array<double,3> 
  EnsembleNVShear::getReducedEnsembleVals() const
  {
    boost::array<double,3> retval;
    retval[0] = EnsembleVals[0];
    retval[1] = EnsembleVals[1] / Sim->dynamics.units().unitVolume();
    retval[2] = EnsembleVals[2] * Sim->dynamics.units().unitTime();

    return retval;
  }

  void
  EnsembleNECompression::initialise()
  {
    EnsembleVals[0] = Sim->particleList.size();
    EnsembleVals[1] = Sim->dynamics.calcInternalEnergy() 
      + Sim->dynamics.getLiouvillean().getSystemKineticEnergy();
    
    try {
      EnsembleVals[2] = dynamic_cast<const LCompression&>
	(Sim->dynamics.getLiouvillean()).getGrowthRate();
    }
    catch (std::exception&)
      {
	M_throw() << "Compression ensemble requires the use of compression liouvillean";
      }

    dout << "NECompression Ensemble initialised\nN=" << EnsembleVals[0]
	     << "\nE=" << EnsembleVals[1] / Sim->dynamics.units().unitEnergy()
	     << "\nGamma=" << EnsembleVals[2] * Sim->dynamics.units().unitTime() << std::endl;
  }

  boost::array<double,3> 
  EnsembleNECompression::getReducedEnsembleVals() const
  {
    boost::array<double,3> retval;
    retval[0] = EnsembleVals[0];
    retval[1] = EnsembleVals[1] / Sim->dynamics.units().unitEnergy();
    retval[2] = EnsembleVals[2] * Sim->dynamics.units().unitTime();

    return retval;
  }

  void
  EnsembleNTCompression::initialise()
  {
    EnsembleVals[0] = Sim->particleList.size();

    try {
      thermostat = Sim->dynamics.getSystem("Thermostat").get_ptr();
    } catch (std::exception&)
      {
	M_throw() << "Could not find the Thermostat in NVT system";
      }
    
    //Only one kind of thermostat so far!
    if (dynamic_cast<const CSysGhost*>(thermostat) == NULL)
      {
	M_throw() << "Could not upcast thermostat to Andersens";
      }
    
    EnsembleVals[1] = static_cast<const CSysGhost*>
      (thermostat)->getTemperature();
    
    try {
      EnsembleVals[2] = dynamic_cast<const LCompression&>
	(Sim->dynamics.getLiouvillean()).getGrowthRate();
    }
    catch (std::exception&)
      {
	M_throw() << "Compression ensemble requires the use of compression liouvillean";
      }

    dout << "NTCompression Ensemble initialised\nN=" << EnsembleVals[0]
	     << "\nT=" << EnsembleVals[1] / Sim->dynamics.units().unitEnergy()
	     << "\nGamma=" << EnsembleVals[2] * Sim->dynamics.units().unitTime() << std::endl;
  }

  boost::array<double,3> 
  EnsembleNTCompression::getReducedEnsembleVals() const
  {
    boost::array<double,3> retval;
    retval[0] = EnsembleVals[0];
    retval[1] = EnsembleVals[1] / Sim->dynamics.units().unitEnergy();
    retval[2] = EnsembleVals[2] * Sim->dynamics.units().unitTime();

    return retval;
  }

}
