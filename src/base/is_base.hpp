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

#include "is_stream_op.hpp"
#include "constants.hpp"
#include <magnet/exception.hpp>


namespace DYNAMO
{
  class SimData;

  /*! \brief Associates a name and colour with a class.
   * 
   * This merely helps format screen output from a class.
   */
  class Base_Class
  {
  public:
    
    /*! \brief Initialises the Base_Class
     *
     * \param aName The name of the class.
     * \param aColor A terminal colour code to colour this classes name with.
     */
    Base_Class(const std::string aName, const std::string aColor):
      name(aName),color(aColor) {};
    
    /*! \brief A private stream to format the standard output stream. */
    Stream_Operator I_cout() const
    {
      return ((std::cout << Stream_Operator(name,color)) << "\n");
    }

    /*! \brief A private stream to format the standard error stream. */
    Stream_Operator I_cerr() const
    {
      return ((std::cerr << Stream_Operator(name,IC_red)) << "\n");
    }
    
  protected:
    //! This constructor is only available for virtual
    //! inheritance. The concrete derived class must call the other
    //! constructor.
    Base_Class() { M_throw() << "Calling the default constructor!"; }

    /*! \brief A pointer to a const definition of the class name. */
    std::string name;

    /*! \brief A pointer to a const definition of the class terminal
        colour. */
    std::string color;
  };

  /*! \brief A Base_Class which contains a writable pointer to a
   * SimData structure.
   *
   * This class must be able to change the SimData struct it points to.
   */
  class SimBase: public Base_Class
  {
  public:
    /*! \brief Constructor
     *
     * \param SD Pointer to the SimData struct
     * \param aName The name of the class deriving from this
     * \param aColor The colour of the output from this class.
     */
    SimBase(SimData* const SD,const std::string aName, 
	    const std::string aColor):
      Base_Class(aName, aColor),
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
   * pointer to a SimData struct.
   *
   * This class must be able to change the SimData struct it points to.
   */
  class SimBase_const: public Base_Class
  {
  public:
    /*! \brief Constructor
     *
     * \param SD Const pointer to the SimData struct
     * \param aName The name of the class deriving from this
     * \param aColor The colour of the output from this class.
     */
    SimBase_const(const SimData* const SD, const std::string aName, 
		  const std::string aColor):
      Base_Class(aName,aColor),
      Sim(SD)
    {};

  protected:
    //! This constructor is only available for virtual
    //! inheritance. The concrete derived class must call the other
    //! constructor.
    SimBase_const() { M_throw() << "Calling the default constructor!"; }

    /*! \brief A un-writable pointer to a simulations data.*/
    const SimData* Sim;
  };

}
