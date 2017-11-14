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

#include <gtkmm.h>
#include <coil/RenderObj/RenderObj.hpp>
#include <magnet/GL/texture.hpp>
#include <magnet/GL/shader/sphere.hpp>
#include <magnet/GL/buffer.hpp>
#include <magnet/GL/camera.hpp>
#include <memory>
#include <array>

namespace coil {
  class RLight : public RenderObj, public magnet::GL::CameraHeadTracking
  {
  public:
    RLight(std::string name, magnet::math::Vector position, magnet::math::Vector lookAtPoint,
	   GLfloat zNearDist = 8.0f, GLfloat zFarDist = 10000.0f, magnet::math::Vector up = magnet::math::Vector{0,1,0},
	   GLfloat simLength = 25.0f, GLfloat size = 1.0):
      RenderObj(name),
      CameraHeadTracking(position, lookAtPoint, zNearDist, zFarDist, up, simLength, magnet::math::Vector{0,0,20}),
      _intensity(1.0 / simLength),
      _specularExponent(32),
      _specularFactor(1),
      _maxvariance(0.1),
      _bleedreduction(0.2),
      _size(size / simLength)
    {
      std::array<GLfloat, 3> tmp = {{1.0, 1.0, 1.0}};
      _color = tmp;
      _shadowCasting = false;
    }
  
    virtual void init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue);
    virtual void deinit();

    virtual void glRender(const magnet::GL::Camera& cam, RenderMode mode, const uint32_t offset);

    virtual void interfaceRender(const magnet::GL::Camera& camera, magnet::GL::objects::CairoSurface& cairo);

    virtual void showControls(Gtk::ScrolledWindow* win);

    void setIntensity(double);
    float getIntensity() const { return _intensity; }
    float getSpecularExponent() const { return _specularExponent; }
    float getSpecularFactor() const { return _specularFactor; }
    float getMaxVariance() const { return _maxvariance; }
    float getBleedReduction() const { return _bleedreduction; }    

    virtual bool deletable() { return true; }

    virtual void dragCallback(Vector cursorPos, uint32_t objID);

    const std::array<GLfloat, 3>& getColor() const { return _color; }
    std::array<GLfloat, 3> getLightColor() const 
    { 
      std::array<GLfloat, 3> retval = {{_color[0] * _intensity,
					     _color[1] * _intensity,
					     _color[2] * _intensity}};
      return retval;
    }

    void setPosition(magnet::math::Vector newposition);

    /*! \brief Load the specified OpenGL texture matrix with the
      projection required for shadow mapping.
      
      \note The current OpenGL model view matrix must be the matrix
      used for rendering.
      
      \param textureUnit The texture unit whose matrix is to be
      setup for shadowmapping.
    */
    inline magnet::GL::GLMatrix getShadowTextureMatrix()
    {
      return magnet::GL::translate(magnet::math::Vector{0.5, 0.5, 0.5})
	* magnet::GL::scale(magnet::math::Vector{0.5, 0.5, 0.5})
	* getProjectionMatrix()
	* getViewMatrix();
    }
    
    /*! \brief Returns a projected light position.
     */
    magnet::math::Vector getEyespacePosition(const magnet::GL::Camera& camera) const
    {
      magnet::math::Vector vec = getPosition();
      magnet::math::NVector<GLfloat, 4> lightPos = {{GLfloat(vec[0]), GLfloat(vec[1]), GLfloat(vec[2]), 1.0f}};
      magnet::math::NVector<GLfloat, 4> lightPos_eyespace
	= camera.getViewMatrix() * lightPos;
      return magnet::math::Vector{lightPos_eyespace[0], lightPos_eyespace[1], lightPos_eyespace[2]};
    }

    virtual uint32_t pickableObjectCount()
    { return visible(); }

    virtual std::string getCursorText(uint32_t objID)
    { return _name; }

    virtual magnet::math::NVector<GLfloat, 4> getCursorPosition(uint32_t objID)
    {
      magnet::math::Vector loc = getPosition();
      return magnet::math::NVector<GLfloat, 4>{GLfloat(loc[0]), GLfloat(loc[1]), GLfloat(loc[2]), 1.0f};
    }

    void setSize(double val);
    GLfloat getSize() const { return _size; }

  protected:
    void initGTK();
    void guiUpdate();

    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    //GTK gui stuff
    std::unique_ptr<Gtk::VBox> _optList;
    std::unique_ptr<Gtk::Entry> _intensityEntry;
    std::unique_ptr<Gtk::ColorButton> _lightColor;
    std::unique_ptr<Gtk::Entry> _specularExponentEntry;
    std::unique_ptr<Gtk::Entry> _specularFactorEntry;
    std::unique_ptr<Gtk::Entry> _positionXEntry;
    std::unique_ptr<Gtk::Entry> _positionYEntry;
    std::unique_ptr<Gtk::Entry> _positionZEntry;    
    std::unique_ptr<Gtk::Entry> _sizeEntry;

    std::unique_ptr<Gtk::Entry> _maxvarianceEntry;
    std::unique_ptr<Gtk::Entry> _bleedreductionEntry;

    float _intensity, _specularExponent, _specularFactor;
    float _maxvariance, _bleedreduction;
    GLfloat _size;
    std::array<GLfloat, 3> _color;

    magnet::GL::shader::SphereShader _sphereShader;
    magnet::GL::Buffer<GLfloat> _glposition;
    magnet::GL::Context::ContextPtr _context;
  };
}
