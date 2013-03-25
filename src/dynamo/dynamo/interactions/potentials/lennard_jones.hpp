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

#pragma once
#include <dynamo/interactions/potentials/potential.hpp>
#include <limits>

namespace dynamo {
    /*! \brief The Lennard-Jones potential.
    
      This class implements an automatically stepped potential where
      each step is specified through one of many algorithms.
   */
  class PotentialLennardJones :public Potential {
  public:
    /*! \brief An enum of the types of step energy algorithms available.*/
    enum UMode {
      MIDPOINT,
      LEFT,
      RIGHT,
      VOLUME,
      VIRIAL
    };
    
    /*! \brief An enum of types of step positionining algorithms available.*/
    enum RMode {
      DELTAR,
      DELTAU
    };

    PotentialLennardJones(double sigma, double epsilon, double cutoff, UMode umode, RMode rmode, double attractivesteps, double kT=1);

    PotentialLennardJones(const magnet::xml::Node& XML) {
      operator<<(XML);
    }
      
    virtual std::size_t steps() const;

    virtual void operator<<(const magnet::xml::Node&);
  
    virtual double hard_core_diameter() const {
      return 0;
    }

    virtual double render_diameter() const {
      return _sigma;
    }

    double U(double) const;
    double U_uncut(double) const;
    double minimum() const;
    
  protected:

    virtual void calculateToStep(size_t) const;

    virtual void outputXML(magnet::xml::XmlStream&) const;

    double B2func(double) const;

    double _sigma;
    double _epsilon;
    double _cutoff;
    double _kT;

    /*! \brief The number of steps in the attractive section of the potential.
      
      You may have a fractional number of steps in the potential.
     */
    double _attractiveSteps;

    //! \brief The active step energy algorithm.
    UMode _U_mode;
    //! \brief The active step position algorithm.
    RMode _R_mode;
  };
}
