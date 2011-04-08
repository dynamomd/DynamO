/*  DYNAMO:- Event driven molecular dynamics simulator 
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
 * \brief Contains the definition of the CEnsemble class.
 */
#pragma once

#include <boost/array.hpp>
#include <string>
#include "constants.hpp"
#include "is_base.hpp"

namespace xml { class XmlStream; }
namespace magnet { namespace xml { class Node; } }

class System;

namespace DYNAMO {
  class SimData;

  /*! \brief This class specifies the simulation ensemble that the
   * simulation is being performed in.
   *
   *
   */
  class CEnsemble : public SimBase_const
  {
  public:
    CEnsemble(const SimData* const& SD, const char *aName):
      SimBase_const(SD, aName, IC_blue) {}
    
    virtual ~CEnsemble() {}

    static CEnsemble* getClass(const magnet::xml::Node& XML, 
			       const DYNAMO::SimData* Sim);

    friend xml::XmlStream& operator<<(xml::XmlStream&, const CEnsemble&);
    
    virtual std::string getName() const = 0;

    virtual void initialise() = 0;

    virtual boost::array<double,3> getReducedEnsembleVals() const = 0;
    
    virtual void exchange(CEnsemble& rhs) { std::swap(EnsembleVals, rhs.EnsembleVals); }

    virtual double exchangeProbability(const CEnsemble&) const;
    
    const boost::array<double,3>& getEnsembleVals() const { return EnsembleVals; }

  protected:
    boost::array<double,3> EnsembleVals;
  };
  
  class CENVE           : public CEnsemble 
  {
  public:
    CENVE(const DYNAMO::SimData* SD): 
      CEnsemble(SD, "CENVE") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NVE"); }
  };

  class CENVT           : public CEnsemble 
  {
  public:
    CENVT(const DYNAMO::SimData* SD): 
      CEnsemble(SD, "CENVT") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NVT"); }

    virtual double exchangeProbability(const CEnsemble&) const;

  protected:
    const System* thermostat;
  };

  class CENVShear       : public CEnsemble 
  {
  public:
    CENVShear(const DYNAMO::SimData* SD): 
      CEnsemble(SD, "CENVShear") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NVShear"); }
  };

  class CENECompression : public CEnsemble 
  {
  public:
    CENECompression(const DYNAMO::SimData* SD): 
      CEnsemble(SD, "CENECompression") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NECompression"); }
  };

  class CENTCompression : public CEnsemble 
  {
  public:
    CENTCompression(const DYNAMO::SimData* SD): 
      CEnsemble(SD, "CENTCompression") {}

    virtual void initialise();

    virtual boost::array<double,3> getReducedEnsembleVals() const;

    virtual std::string getName() const
    { return std::string("NTCompression"); }

  protected:
    const System* thermostat;
  };
}
