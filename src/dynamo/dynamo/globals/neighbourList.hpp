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
#include <dynamo/globals/global.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/ranges/IDRangeList.hpp>
#include <magnet/function/delegate.hpp>
#include <magnet/math/vector.hpp>
#include <memory>
#include <vector>

namespace dynamo {
  /*! \brief A base class for Global events which implement a neighbour list.
    
    This is the interface for neighbour lists, which are used to
    optimise the look up of particles in the neighbourhood of a given
    \ref Particle.
   
    This class also defines callback's that can be registered so that
    other parts of DynamO can be updated when a particle changes
    neighbours.
   */
  class GNeighbourList: public Global, magnet::Tracked
  {
  public:
    GNeighbourList(dynamo::Simulation* a, const char *b): 
      Global(a, b),
      _initialised(false),
      _maxInteractionRange(0),
      isUsedInScheduler(false),
      lambda(0.9)
    {}

    virtual void getParticleNeighbours(const Particle&, std::vector<size_t>&) const = 0;
    virtual void getParticleNeighbours(const Vector&, std::vector<size_t>&) const = 0;

    /*! \brief This returns the maximum interaction length this
      neighbourlist supports.
      
      Due the to neighbourlists using integer numbers of cells, they
      end up supporting an interaction range larger than 
      \ref getMaxInteractionRange().
    */
    virtual double
    getMaxSupportedInteractionLength() const = 0;

    virtual void reinitialise()
    {
      if (!_maxInteractionRange)
	_maxInteractionRange = Sim->getLongestInteraction();

      _initialised = true;
    }
    
    void markAsUsedInScheduler() { isUsedInScheduler = true; }

    void setCellOverlap(bool overlap) 
    {
      if (overlap)
	lambda = 0.9; 
      else
	lambda = 0.001;
    }

    /*! \brief Set the minimum range this neighbourlist is to support.
      
      This is the minimum as neighbourlists usually must support a
      slightly larger distance.
      
      \sa getMaxSupportedInteractionLength()
     */
    void setMaxInteractionRange(double range)
    {
      _maxInteractionRange = range;
      if (_initialised) reinitialise();
    }

    /*! \brief Returns the requested minimum supported interaction
        range.
     */
    double getMaxInteractionRange() const
    { return _maxInteractionRange; }

    mutable magnet::Signal<void(const Particle&, const size_t&)> _sigNewNeighbour;
    mutable magnet::Signal<void(const Particle&, const size_t&)> _sigCellChange;
    mutable magnet::Signal<void()> _sigReInitialise;

  protected:
    bool _initialised;
    double _maxInteractionRange;

    GNeighbourList(const GNeighbourList&);

    virtual void outputXML(magnet::xml::XmlStream&) const = 0;

    bool isUsedInScheduler;
    double lambda; 
  };
}

