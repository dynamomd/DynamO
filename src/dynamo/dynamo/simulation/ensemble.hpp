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
#pragma once

#include <dynamo/base.hpp>
#include <boost/array.hpp>
#include <string>

namespace magnet { namespace xml { class Node; class XmlStream; } }

class System;

namespace dynamo {
  class SimData;

  /*! \brief This class specifies the simulation ensemble that the
   * simulation is being performed in.
   *
   * This is the abstract base class defining the Ensemble interface.
   */
  class Ensemble : public SimBase_const
  {
  public:
    Ensemble(const SimData* const& SD, const char *aName):
      SimBase_const(SD, aName) {}
    
    virtual ~Ensemble() {}

    //! Used to load an xml tag corresponding to an Ensemble and
    //! generate the correct derived Ensemble type.
    //! \param XML The XML node containing the type of the Ensemble
    //! \param Sim A pointer to the simulation data this Ensemble is valid for.
    static Ensemble* getClass(const magnet::xml::Node& XML, 
			       const dynamo::SimData* Sim);

    //! Helper function which generates an Ensemble XML tag.
    //! \sa getName
    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const Ensemble&);
    
    //! This function returns a string corresponding to the type of
    //! Ensemble. This name is used to load the Ensemble again from
    //! the config file.
    virtual std::string getName() const = 0;

    //! Called to generate and store the Ensemble variables.
    virtual void initialise() = 0;

    //! Returns an array containing the control values of the Ensemble
    //! (e.g., NVE) in the units of the output.
    //! \sa getEnsembleVals
    virtual boost::array<double,3> getReducedEnsembleVals() const = 0;
    
    //! Swaps the underlying ensemble control values.
    //! \sa EReplicaExchangeSimulation
    virtual void swap(Ensemble& rhs) { std::swap(EnsembleVals, rhs.EnsembleVals); }

    //! Calculates the probability of carrying out a replica exchange
    //! move between this Ensemble and another.
    virtual double exchangeProbability(const Ensemble&) const;
    
    //! Returns an array containing the ensemble values in simulation units.
    //! \sa getReducedEnsembleVals
    const boost::array<double,3>& getEnsembleVals() const { return EnsembleVals; }

  protected:
    boost::array<double,3> EnsembleVals;
  };

  //! \brief An Ensemble where N (no. of particles), V (simulation volume),
  //! and E (total energy) are held constant.
  class EnsembleNVE           : public Ensemble 
  {
  public:
    EnsembleNVE(const dynamo::SimData* SD): 
      Ensemble(SD, "EnsembleNVE") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NVE"); }
  };

  //! \brief An Ensemble where N (no. of particles), V (simulation volume),
  //! and T (temperature) are held constant.
  //! 
  //! This class also stores a pointer to the thermostat used to hold
  //! the temperature constant
  class EnsembleNVT           : public Ensemble 
  {
  public:
    EnsembleNVT(const dynamo::SimData* SD): 
      Ensemble(SD, "EnsembleNVT") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NVT"); }

    virtual double exchangeProbability(const Ensemble&) const;

  protected:
    const System* thermostat;
  };

  //! \brief An Ensemble where N (no. of particles), V (simulation volume),
  //! and shear rate are held constant.
  class EnsembleNVShear       : public Ensemble 
  {
  public:
    EnsembleNVShear(const dynamo::SimData* SD): 
      Ensemble(SD, "EnsembleNVShear") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NVShear"); }
  };

  //! \brief An Ensemble where N (no. of particles), E (total energy),
  //! and isotropic compression rate are held constant.
  class EnsembleNECompression : public Ensemble 
  {
  public:
    EnsembleNECompression(const dynamo::SimData* SD): 
      Ensemble(SD, "EnsembleNECompression") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NECompression"); }
  };

  //! \brief An Ensemble where N (no. of particles), T (Temperature),
  //! and isotropic compression rate are held constant.
  //! 
  //! This class also stores a pointer to the thermostat used to hold
  //! the temperature constant.
  class EnsembleNTCompression : public Ensemble 
  {
  public:
    EnsembleNTCompression(const dynamo::SimData* SD): 
      Ensemble(SD, "EnsembleNTCompression") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NTCompression"); }

  protected:
    const System* thermostat;
  };
}
