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
#include <magnet/stream/console_specials.hpp>
#include <magnet/stream/formattedostream.hpp>
#include <magnet/exception.hpp>
#include <memory>
#include <iostream>
#include <iomanip>
#include <limits>
#include <functional>

namespace dynamo
{
  using std::shared_ptr;

  class Simulation;

  /*! \brief Provides some basic IO functionality to a derived class.
    
    This class is the Base class for most of the classes in
    DynamO. Its purpose is to provide some helpful functionality,
    such as formatted screen output.
   */
  class Base
  {
  public:
    /*! \brief Initialises the Base class.
     
      \param aName The name of the class.
     */
    Base(const std::string name):
      dout(std::cout), derr(std::cerr)
    {
      setOutputPrefix(name);
      //Reasonable precision for output
      dout << std::setprecision(std::numeric_limits<float>::digits10);
      derr << std::setprecision(std::numeric_limits<float>::digits10);
    }

  protected:
    void setOutputPrefix(const std::string& prefix)
    {
      dout.setPrefix(colorCode(prefix) + prefix + ": " + magnet::console::reset());
      derr.setPrefix(magnet::console::bold() + magnet::console::red_fg() + prefix + ": " + magnet::console::reset());
    }

    /*! \brief A std::cout style output stream. 
     
      This member is meant as a replacement to std::cout, as it
      provides automatic formatting of the output.
     
      \note Before any output will appear on the screen, the stream
      must be flushed. The most convenient way of doing this is to
      always end your output with a std::endl like so:
      \code dout << "An example output" << std::endl; \endcode
     */
    mutable magnet::stream::FormattedOStream dout;

    /*! \brief See \ref dout for more information. */
    mutable magnet::stream::FormattedOStream derr;

    /*! \brief This constructor is only available for virtual
      inheritance. The concrete derived class must call the other
      constructor.
    */
    Base() { M_throw() << "Calling the default constructor!"; }

  private:
    /*! \brief Generate a random console text-color command based off
      a string.
     
      This function is used to automatically pick a color for the
      formatted output of a class, by using a hash of the classes
      name.
     */
    std::string colorCode(std::string str)
    {
      switch (std::hash<std::string>()(str) % 9)
	{
	case 0: return magnet::console::cyan_fg();
	case 1: return magnet::console::purple_fg();
	case 2: return magnet::console::blue_fg();
	case 3: return magnet::console::yellow_fg();
	case 4: return magnet::console::green_fg();
	case 5: return magnet::console::bold() + magnet::console::green_fg();
	case 6: return magnet::console::bold() + magnet::console::blue_fg();
	case 7: return magnet::console::bold() + magnet::console::purple_fg();
	case 8: return magnet::console::bold() + magnet::console::cyan_fg();
	}
      return "";
    }
  };

  /*! \brief A Base class which contains a writable pointer to a
    Simulation structure.
   
    This class must be able to change the Simulation struct it points
    to.
   */
  class SimBase: public Base
  {
  public:
    /*! \brief Constructor
     
      \param SD Pointer to the Simulation class.
      \param aName The name of the class deriving from this class.
     */
    SimBase(Simulation* const SD, const std::string aName):
      Base(aName),
      Sim(SD)
    {};
    
    Simulation* const getSimPointer() { return Sim; }

  protected:
    /*! \brief This constructor is only available for virtual
      inheritance. The concrete derived class must call the other
      constructor.
    */
    SimBase();

    /*! \brief A pointer to a Simulation, which allows the simulation
        to be altered.*/
    Simulation* const Sim;
  };
  
  /*! \brief Similar to the SimBase class except it contains a const
    pointer to a Simulation class.
   
    This class must be able to change the Simulation struct it points to.
   */
  class SimBase_const: public Base
  {
  public:
    /*! \brief Constructor
     
      \param SD Const pointer to the Simulation struct
      \param aName The name of the class deriving from this
      \param aColor The colour of the output from this class.
     */
    SimBase_const(const Simulation* const SD, const std::string aName):
      Base(aName),
      Sim(SD)
    {};

    /*! \brief A pointer to a Simulation, which does not allow the
      simulation to be altered.*/
    const Simulation * const  getSimPointer() { return Sim; }

  protected:
    /*! \brief This constructor is only available for virtual
     inheritance. The concrete derived class must call the other
     constructor.
    */
    SimBase_const();

    /*! \brief A const pointer to a Simulation class.*/
    const Simulation * const Sim;
  };

}
