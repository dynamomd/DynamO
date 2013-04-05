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
#include <vector>
#include <dynamo/base.hpp>

namespace magnet { namespace xml { class Node; class XmlStream; } }

namespace dynamo {
  /*! \brief The base class for any stepped potential.
    
    This class represents a general stepped potential. Each "step" is
    represented by a pairing of distance and energy. Depending on the
    nature of the potential, this energy may correspond to the left or
    right of the discontinuity (see direction()). In short, if a
    particle is on a step with an ID of zero, it has an interaction
    energy of zero.

    This class also implements a cache, to allow fast lookup of previously
    accessed steps, as some calculated potentials are expensive to
    compute.
   */
  class Potential {
  public:
    
    typedef std::pair<double, double> value_type;

    /*! \brief Accessor to give a value_type containing the
        discontinuity location and energy change.

	If the step is not in the step cache, this function calculates
	new steps and adds them to the cache until the requested step
	ID is located.
     */
    virtual value_type operator[](const std::size_t step_id) const {
#ifdef DYNAMO_DEBUG
      if (step_id >= steps()) M_throw() << "Out of range access";
#endif 

      if (step_id >= cached_steps()) calculateToStep(step_id);

      return value_type(_r_cache[step_id], _u_cache[step_id]);
    }
    
    /*!\brief Return the maximum number of steps in the potential.
     */
    virtual std::size_t steps() const = 0;

    /*!\brief Return the hard core diameter of the potential, or zero
       if there is no hard core.
     */
    virtual double hard_core_diameter() const = 0; 

    /*!\brief Return the diameter which should be used to render/draw
       the particle as a sphere.
     */
    virtual double render_diameter() const = 0;
    
    /*!\brief Return how many steps of the potential have already been
      calculated and cached.
    */
    std::size_t cached_steps() const {
      return std::min(_r_cache.size(), _u_cache.size());
    }

    static shared_ptr<Potential> getClass(const magnet::xml::Node&);

    /*! \brief Loads the Potential from an XML node in a
        configuration file.
     */
    virtual void operator<<(const magnet::xml::Node&) = 0;
  
    /*! \brief A helper function that calls outputXML to write out the
        parameters of this Potential to a config file.
     */
    friend magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream&, const Potential&);

    /*! \brief Determine which step in the potential the passed radius
        corresponds to.
    */
    size_t calculateStepID(const double r) const {
      size_t retval(0);
      if (direction())
	for (; (retval < steps()) && (r > operator[](retval).first); ++retval) {}
      else
	for (; (retval < steps()) && (r < operator[](retval).first); ++retval) {}
      return retval;
    }

    /*! \brief Return a pair with the min-max bounds of the potential
        step ID given.
     */
    std::pair<double, double> getStepBounds(size_t ID) const {
      if (direction())
	return std::pair<double, double>((ID == 0) ? 0: operator[](ID-1).first, operator[](ID).first);
      else
	return std::pair<double, double>(operator[](ID).first, ((ID == 0) ? HUGE_VAL : operator[](ID - 1).first));
    }
    
    /*! \brief Returns false if the discontinuity energies specify the
	step energy to its left or true if it is the energy to its
	right.
	
	If a potential limits to infinity at zero separation
	(e.g. Lennard-Jones), then it is natural to start the stepping
	ID's at the cut-off radius and have them increase as r
	approaches zero. This implies that the way in which each step
	is stored is through a discontinuity position (r) and the
	value of the energy to its left (r^-). The reverse is true if
	the potential is zero at r=0, and diverges as r increases,
	each discontinuity should be stored as its position (r) and
	the energy to its right (r^+).
     */
    virtual bool direction() const { return false;}

  protected:
    virtual void calculateToStep(size_t) const = 0;
    virtual void outputXML(magnet::xml::XmlStream&) const = 0;

    mutable std::vector<double> _r_cache;
    mutable std::vector<double> _u_cache;
  };

  /*! \brief A manually stepped potential.
    
    This class implements a stepped potential where each step is
    specified manually.
   */
  class PotentialStepped :public Potential {
  public:
    PotentialStepped(const magnet::xml::Node& XML) {
      operator<<(XML);
    }
    
    PotentialStepped(std::vector<std::pair<double, double> >);
  
    virtual std::size_t steps() const {
      return _r_cache.size();
    }

    virtual void operator<<(const magnet::xml::Node&);
  
    virtual double hard_core_diameter() const {
      for (size_t i(0); i < _u_cache.size(); ++i)
	if (std::isinf(_u_cache[i]))
	  return _r_cache[i];
      
      return 0;
    }
    
    virtual double render_diameter() const {
      double hard_core_d = hard_core_diameter();

      if (hard_core_d != 0) return hard_core_d;
      //If there is no hard core, just return the innermost step diameter
      return _r_cache.back();
    }

  protected:
    virtual void calculateToStep(size_t) const {
      M_throw() << "Cannot calculate new steps for this potential!";
    }

    virtual void outputXML(magnet::xml::XmlStream&) const;
  };
}
