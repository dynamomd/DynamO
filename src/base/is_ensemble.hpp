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
/*! \file is_ensemble.hpp
 * \brief Contains the definition of the Ensemble class.
 */
#pragma once

#include <boost/array.hpp>
#include <string>
#include "constants.hpp"
#include "is_base.hpp"

namespace xml { class XmlStream; }
namespace magnet { namespace xml { class Node; } }

class System;

namespace dynamo {
  class SimData;

  /*! \brief This class specifies the simulation ensemble that the
   * simulation is being performed in.
   *
   *
   */
  class Ensemble : public SimBase_const
  {
  public:
    Ensemble(const SimData* const& SD, const char *aName):
      SimBase_const(SD, aName, IC_blue) {}
    
    virtual ~Ensemble() {}

    static Ensemble* getClass(const magnet::xml::Node& XML, 
			       const dynamo::SimData* Sim);

    friend xml::XmlStream& operator<<(xml::XmlStream&, const Ensemble&);
    
    virtual std::string getName() const = 0;

    virtual void initialise() = 0;

    virtual boost::array<double,3> getReducedEnsembleVals() const = 0;
    
    virtual void exchange(Ensemble& rhs) { std::swap(EnsembleVals, rhs.EnsembleVals); }

    virtual double exchangeProbability(const Ensemble&) const;
    
    const boost::array<double,3>& getEnsembleVals() const { return EnsembleVals; }

  protected:
    boost::array<double,3> EnsembleVals;
  };
  
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
