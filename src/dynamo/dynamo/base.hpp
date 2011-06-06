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
#include <magnet/exception.hpp>
#include <iostream>

namespace dynamo
{
  class SimData;

  /*! \brief Provides some basic IO functionality to a derived class.
   * 
   * This class is the Base class for most of the classes in
   * DynamO. Its purpose is to provide some helpful functionality,
   * such as formatted screen output.
   */
  class Base
  {
  public:
    /*! \brief Initialises the Base class.
     *
     * \param aName The name of the class.
     */
    Base(const std::string aName):
      name(aName) {};
    
    /*! \brief A private stream to format the standard output stream. */
    std::ostream& I_cout() const
    {
      return (std::cout << "\n");
    }

    /*! \brief A private stream to format the standard error stream. */
    std::ostream& I_cerr() const
    {
      return (std::cerr << "\n");
    }
    
  protected:
    //! This constructor is only available for virtual
    //! inheritance. The concrete derived class must call the other
    //! constructor.
    Base() { M_throw() << "Calling the default constructor!"; }

    /*! \brief A pointer to a const definition of the class name. */
    std::string name;
  };

  /*! \brief A Base class which contains a writable pointer to a
   * SimData structure.
   *
   * This class must be able to change the SimData struct it points
   * to.
   */
  class SimBase: public Base
  {
  public:
    /*! \brief Constructor
     *
     * \param SD Pointer to the SimData class.
     * \param aName The name of the class deriving from this class.
     */
    SimBase(SimData* const SD,
	    const std::string aName):
      Base(aName),
      Sim(SD)
    {};
    
  protected:
    //! This constructor is only available for virtual
    //! inheritance. The concrete derived class must call the other
    //! constructor.
    SimBase() { M_throw() << "Calling the default constructor!"; }

    /*! \brief A writable pointer to a simulations data.*/
    SimData* Sim;
  };
  
  /*! \brief Similar to the SimBase class except it contains a const
   * pointer to a SimData class.
   *
   * This class must be able to change the SimData struct it points to.
   */
  class SimBase_const: public Base
  {
  public:
    /*! \brief Constructor
     *
     * \param SD Const pointer to the SimData struct
     * \param aName The name of the class deriving from this
     * \param aColor The colour of the output from this class.
     */
    SimBase_const(const SimData* const SD, const std::string aName):
      Base(aName),
      Sim(SD)
    {};

  protected:
    //! This constructor is only available for virtual
    //! inheritance. The concrete derived class must call the other
    //! constructor.
    SimBase_const() { M_throw() << "Calling the default constructor!"; }

    /*! \brief A const pointer to a SimData class.*/
    const SimData* Sim;
  };

}
