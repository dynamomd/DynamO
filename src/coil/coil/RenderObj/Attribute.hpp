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
#include <magnet/GL/buffer.hpp>
#include <vector>

namespace coil {
  /*! \brief This class encapsulates attributes (data) associated with
   * some topology.
   *
   * This class is the primary communication interface between a
   * simulation and the coil library. After the visualizer is
   * initialised, all data to be rendered should be passed through
   * this class.
   * 
   * The topology may be a collection of points or cells and the data
   * may be ordinates (positions of the points), extensive properties
   * (like the mass) or intensive properties (like the density). Some
   * data is scalar (like the temperature) and some data will have
   * several components per value (e.g. vector quantities like the
   * velocity).
   *
   */
  class Attribute
  {
  public:
    enum AttributeType { 
      INTENSIVE, //!< Intensive property (e.g., Temperature, density)
      EXTENSIVE, //!< Extensive property (e.g., mass, momentum)
      COORDINATE //!< A special attribute which specifies the location of the attribute.
    };

    Attribute(size_t N, AttributeType type = EXTENSIVE, size_t components = 1):
      _hostData(N * components)
      _components(components),
      _type(type),
      _references(0)
    {}
    
    /*! \brief Initialises the OpenGL resources of this object.
     */
    void init() { _glData.init(_hostData); }

    /*! \brief Releases the OpenGL resources of this object.
     */
    void deinit() { _glData.deinit(); }

    /*! \brief Returns the GL buffer associated with the Attribute
     * data.
     */
    magnet::GL::Buffer<GLfloat>& getBuffer() { return _glData; }

    /** @name The host code interface. */
    /**@{*/
    
    /*! \brief Returns a reference to the host-cache of the Attribute
     * data.
     *
     * The Attribute data may be directly updated by the host program,
     * but flagNewData() must be called for the update to take effect.
     */
    std::vector<GLfloat>& getData() { return _hostData; }
    
    /*! \brief Marks that the data in the buffer has been updated, and
     * should be uploaded to the GL system.
     *
     * This function just inserts a callback in the GL system to
     * reinitialise the Attribute.
     */
    void flagNewData()
    { _glData.getContext().queueTask(magnet::function::Task::makeTask(&Attribute::init, this)); }

    /*! \brief Test if the attribute is in use and should be
     * updated. 
     */
    bool active() const { return _references; }

    /*! \brief Returns the number of elements.
     */
    size_t size() const { return _hostData.size() / _components; }

    /**@}*/

  protected:
    /*! \brief The OpenGL representation of the attribute data.
     *
     * There are N * _components floats of attribute data.
     */
    magnet::GL::Buffer<GLfloat> _glData;
    
    /*! \brief A host side cache of \ref _glData.
     *
     * This is used as a communication buffer, both when the host
     * program is writing into coil, and when coil passes the data
     * into OpenGL.
     */
    std::vector<GLfloat> _hostData;

    //! \brief The number of components per value.
    size_t _components;

    //! \brief The type of data stored in this Attribute.
    AttributeType _type;

    /*! \brief The number of glyphs, filters and other render objects
     * currently using this type.
     */
    size_t _references;
  };
}
