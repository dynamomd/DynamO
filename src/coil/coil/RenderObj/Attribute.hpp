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
    some topology.
   
    This class is the primary communication interface between a
    simulation and the coil library. After the visualizer is
    initialised, all data to be rendered should be passed through this
    class.
    
    The topology may be a collection of points or cells and the data
    may be ordinates (positions of the points), extensive properties
    (like the mass) or intensive properties (like the density). Some
    data is scalar (like the temperature) and some data will have
    several components per value (e.g. vector quantities like the
    velocity).

    The inherited container is used as a communication buffer, both
    when the host program is writing into coil, and when coil passes
    the data into OpenGL.
  */
  class Attribute : public std::vector<GLfloat>
  {
  public:
    enum AttributeType { 
      INTENSIVE  = 1 << 0, //!< Intensive property (e.g., Temperature, density)
      EXTENSIVE  = 1 << 1, //!< Extensive property (e.g., mass, momentum)
      COORDINATE = 1 << 2, //!< A special attribute which specifies the location of the attribute.
      
      DEFAULT_GLYPH_POSITION = 1 << 3, //!< This flag marks that the attribute should be used as the initial position value for a glyph.
      DEFAULT_GLYPH_SCALING  = 1 << 4, //!< This flag marks that the attribute should be used as the initial scaling field for a glyph.
      DEFAULT_GLYPH_ORIENTATION  = 1 << 5 //!< This flag marks that the attribute should be used as the initial orientation for a glyph
    };

    inline Attribute(size_t N, int type, size_t components, 
		     magnet::GL::Context::ContextPtr context):
      std::vector<GLfloat>(N * components),
      _usedInLastRender(false),
      _usedInCurrentRender(false),
      _context(context),
      _dataUpdates(0),
      _components(components),
      _type(type),
      _references(0)
    {
      if (_components > 4)
	M_throw() << "We don't support greater than 4 component attributes due to the way "
		  << "data is sometimes directly passed to the shaders (e.g. positional data)";
    }
    
    /*! \brief Releases the OpenGL resources of this object.
     */
    inline void deinit() { _glData.deinit(); }

    /*! \brief Returns the GL buffer associated with the Attribute
     * data.
     */
    inline magnet::GL::Buffer<GLfloat>& getBuffer()
    { 
      //Initialise on demand
      if (!_glData.size()) initGLData();
      _usedInCurrentRender = true;
      return _glData; 
    }

    inline size_t getUpdateCount() const { return _dataUpdates; }

    /** @name The host code interface. */
    /**@{*/
    
    /*! \brief Marks that the data in the buffer has been updated, and
     * should be uploaded to the GL system.
     *
     * This function just inserts a callback in the GL system to
     * reinitialise the Attribute.
     */
    inline void flagNewData()
    { _context->queueTask(magnet::function::Task::makeTask(&Attribute::initGLData, this)); }

    /*! \brief Test if the attribute is in use and should be
     * updated. 
     */
    inline bool active() const { return _references; }

    /*! \brief Returns the number of elements */
    inline size_t num_elements() const { return size() / _components; }

    inline size_t components() const { return _components; }

    inline int getType() const { return _type; }

    /**@}*/

    inline void bindAttribute(size_t attrnum, bool normalise = false, size_t divisor = 1)
    { getBuffer().attachToAttribute(attrnum, _components, divisor, normalise); }

    inline const std::vector<GLfloat>& minVals() const { return _minVals; }
    inline const std::vector<GLfloat>& maxVals() const { return _maxVals; }

    inline volatile bool inUse() { return _usedInLastRender; }

    inline void renderComplete()
    {
      _usedInLastRender = _usedInCurrentRender;
      _usedInCurrentRender = false;
    }

  protected:
    volatile bool _usedInLastRender;
    volatile bool _usedInCurrentRender;
    
    magnet::GL::Context::ContextPtr _context;
    std::vector<GLfloat> _minVals;
    std::vector<GLfloat> _maxVals;
    
    /*! \brief A function which actually performs the copy of data to
        the OpenGL buffer.
	
	This function must be called in the OpenGL thread and is
	usually invoked as a callback from \ref flagNewData(). This
	function also emits the _glDataUpdated signal for any
	post-upload data processing to occur.
     */
    void initGLData()
    { 
      _glData.init(*this);
      //Increase the updates counter
      ++_dataUpdates;

      //Also generate any statistics we report on the OpenGL data
      size_t comps = components();
      _minVals = std::vector<GLfloat>(begin(), begin() + components());
      _maxVals = std::vector<GLfloat>(begin(), begin() + components());
      for (size_t i = 1; i < num_elements(); ++i)
	for (size_t j = 0; j < comps; ++j)
	  {
	    _minVals[j] = std::min(_minVals[j], (*this)[i * comps + j]);
	    _maxVals[j] = std::max(_maxVals[j], (*this)[i * comps + j]);
	  }
    }

    /*! \brief The OpenGL representation of the attribute data.
     *
     * There are N * _components floats of attribute data.
     */
    magnet::GL::Buffer<GLfloat> _glData;

    /*! \brief A counter of how many updates have been applied to the
      data.
      
      This is used to track when the data has been updated.
     */
    size_t _dataUpdates;

    //! \brief The number of components per value.
    size_t _components;

    //! \brief The type of data stored in this Attribute.
    int _type;

    /*! \brief The number of glyphs, filters and other render objects
     * currently using this type.
     */
    size_t _references;
  };
}
