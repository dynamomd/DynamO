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

#pragma once

#include <boost/rational.hpp>

//!A macro factory for generating the Units to std::string and std::string
//!to Units conversions.
#define UNITS_FACTORY(F) \
  F(Dimensionless) \
  F(Length) \
  F(Area) \
  F(Volume) \
  F(Time) \
  F(Mass) \
  F(Velocity) \
  F(Momentum) \
  F(Energy) \
  F(Diffusion) \
  F(MutualDiffusion) \
  F(Thermalconductivity) \
  F(ThermalDiffusion) \
  F(Viscosity) \
  F(Pressure)
  

namespace magnet {
  namespace units {
    //! \brief This class is used for checking the units of a value at
    //! runtime.  
    //!
    //! This class can be used as a base class to store the units of a
    //! value and to check them at runtime. This runtime requirement means
    //! many of the popular template techniques are unsuitable.
    class Units
    {
      //! The type used to implement the value of each units dimension.
      typedef boost::rational<int> Value;
    public:
      //! Constructor, allowing RAII of the units with any types
      //! supported by Value.
      template<typename T1, typename T2, typename T3>
      inline Units(T1 l, T2 t, T3 m)
      { _unitPowers[L] = l; _unitPowers[T] = t; _unitPowers[M] = m; }
      
      //!Enumeration of the unit dimensions
      enum {L=0, //!> Length units.
	    T=1, //!> Time units.
	    M=2, //!> Mass units.
	    END=3};
      
      //!Allows the generation of new units when two units scales are
      //!multiplied.
      inline Units operator*(const Units& ou) const
      {
	Units retval;
	for (size_t i(0); i < END; ++i)
	  retval._unitPowers[i] = _unitPowers[i] + ou._unitPowers[i];
	return retval;
      }
      
      //!Allows the generation of new units when two units scales are
      //!divided.
      inline Units operator/(const Units& ou) const
      {
	Units retval;
	for (size_t i(0); i < END; ++i)
	  retval._unitPowers[i] = _unitPowers[i] - ou._unitPowers[i];
	return retval;
      }
      
      //!Comparison operator
      inline bool operator==(const Units& ou) const
      {
	for (size_t i(0); i < END; ++i)
	  if (_unitPowers[i] != ou._unitPowers[i]) return false;
	return true;
      }
      
      //The definitions of many units of interest

      //!Returns a handy no-units Units
      static inline Units Dimensionless() { return Units(0,0,0); }

      //!This is Boltzmanns constant which gives units to the
      //!temperature. By setting this to zero units we are essentially
      //!saying temperature is in units of energy. This is only here
      //!to make the expressions below complete
      static inline Units kB() { return Dimensionless(); }

      //!Returns a Units for length
      static inline Units Length() { return Units(1,0,0); }
      //!Returns a Units for area
      static inline Units Area() { return Length() * Length(); }
      //!Returns a Units for volume
      static inline Units Volume() { return Area() * Length(); }
      //!Returns a Units for time
      static inline Units Time() { return Units(0,1,0); }
      //!Returns a Units for mass
      static inline Units Mass() { return Units(0,0,1); }
      //!Returns a Units for velocity
      static inline Units Velocity() { return Length() / Time(); }
      //!Returns a Units for Momentum
      static inline Units Momentum() { return Velocity() * Mass(); }
      //!Returns a Units for energy
      static inline Units Energy() { return Velocity() * Velocity() * Mass(); }

      //!Returns a Units for diffustion
      static inline Units Diffusion() { return Area() / Time(); }
      //!Returns a Units for mutual diffustion
      static inline Units MutualDiffusion() { return Mass() * Time() / Volume(); }
      //!Returns a Units for thermal conductivity
      static inline Units ThermalConductivity() { return kB() / (Time() * Length()); }
      //!Returns a Units for thermal diffusion
      static inline Units ThermalDiffusion() { return Mass() / (Time() * Length()); }
      //!Returns a Units for viscosity
      static inline Units Viscosity() { return  Dimensionless() / (Time() * Length()); }

      //!Returns a Units for pressue
      static inline Units Pressure() { return  Mass() / (Time() * Time() * Length()); }
      
    private:

      inline Units() {}
      Value _unitPowers[END];
    };
  }
}

#undef UNITS_FACTORY
