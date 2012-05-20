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
#include <magnet/gtk/transferFunction.hpp>
#include <magnet/GL/texture.hpp>
#include <magnet/GL/shader/sphere.hpp>
#include <magnet/GL/buffer.hpp>
#include <magnet/GL/camera.hpp>
#include <memory>
#include <tr1/array>

namespace coil {
  class RLight : public RenderObj, public magnet::GL::Camera
  {
  public:
    RLight(std::string name, magnet::math::Vector position, 
	   magnet::math::Vector lookAtPoint, GLfloat fovY = 45.0f,
	   GLfloat zNearDist = 0.05f, GLfloat zFarDist = 10.0f,
	   magnet::math::Vector up = magnet::math::Vector(0,1,0)): 
      RenderObj(name), 
      Camera(1,1,position, lookAtPoint, fovY, zNearDist, zFarDist, up),
      _intensity(1), 
      _specularExponent(96),
      _specularFactor(1),
      _size(0.015)
    {
      std::tr1::array<GLfloat, 3> tmp = {{1.0, 1.0, 1.0}};
      _color = tmp;
    }
  
    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue);
    virtual void deinit();

    virtual void glRender(const magnet::GL::Camera& cam, 
			  RenderMode mode);

    virtual void interfaceRender(const magnet::GL::Camera& camera, magnet::GL::objects::CairoSurface& cairo);

    virtual void showControls(Gtk::ScrolledWindow* win);

    float getIntensity() const { return _intensity; }
    float getSpecularExponent() const { return _specularExponent; }
    float getSpecularFactor() const { return _specularFactor; }

    virtual bool deletable() { return true; }

    virtual void dragCallback(Vector cursorPos, uint32_t objID);

    const std::tr1::array<GLfloat, 3>& getColor() const { return _color; }
    std::tr1::array<GLfloat, 3> getLightColor() const 
    { 
      std::tr1::array<GLfloat, 3> retval = {{_color[0] * _intensity,
					     _color[1] * _intensity,
					     _color[2] * _intensity}};
      return retval;
    }

    /*! \brief Load the specified OpenGL texture matrix with the
      projection required for shadow mapping.
      
      \note The current OpenGL model view matrix must be the matrix
      used for rendering.
      
      \param textureUnit The texture unit whose matrix is to be
      setup for shadowmapping.
    */
    inline magnet::GL::GLMatrix getShadowTextureMatrix()
    {
      return magnet::GL::GLMatrix::translate(magnet::math::Vector(0.5, 0.5, 0.5))
	* magnet::GL::GLMatrix::scale(magnet::math::Vector(0.5, 0.5, 0.5))
	* getProjectionMatrix()
	* getViewMatrix();
    }
    
    /*! \brief Returns a projected light position.
     */
    magnet::math::Vector getEyespacePosition(const magnet::GL::Camera& camera) const
    {
      magnet::math::Vector vec = getPosition();
      std::tr1::array<GLfloat, 4> lightPos = {{vec[0], vec[1], vec[2], 1.0}};
      std::tr1::array<GLfloat, 4> lightPos_eyespace
	= camera.getViewMatrix() * lightPos;
      return magnet::math::Vector(lightPos_eyespace[0], lightPos_eyespace[1], 
				  lightPos_eyespace[2]);
    }

    virtual uint32_t pickableObjectCount()
    { return visible(); }

    virtual void pickingRender(const magnet::GL::Camera& cam, 
			       const uint32_t offset);
    
    virtual std::string getCursorText(uint32_t objID)
    { return _name; }

    virtual std::tr1::array<GLfloat, 4> getCursorPosition(uint32_t objID)
    {
      magnet::math::Vector loc = getPosition();
      std::tr1::array<GLfloat, 4> vec = {{loc[0], loc[1], loc[2], 1.0}};
      return vec;
    }

  protected:
    void initGTK();
    void guiUpdate();

    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    //GTK gui stuff
    std::auto_ptr<Gtk::VBox> _optList;
    std::auto_ptr<Gtk::Entry> _intensityEntry;
    std::auto_ptr<Gtk::ColorButton> _lightColor;
    std::auto_ptr<Gtk::Entry> _specularExponentEntry;
    std::auto_ptr<Gtk::Entry> _specularFactorEntry;
    std::auto_ptr<Gtk::Entry> _positionXEntry;
    std::auto_ptr<Gtk::Entry> _positionYEntry;
    std::auto_ptr<Gtk::Entry> _positionZEntry;    
    std::auto_ptr<Gtk::Entry> _sizeEntry;

    float _intensity, _specularExponent, _specularFactor;
    GLfloat _size;
    std::tr1::array<GLfloat, 3> _color;

    magnet::GL::shader::SphereShader<true> _sphereShader;
    magnet::GL::Buffer<GLfloat> _glposition;
    magnet::GL::Context::ContextPtr _context;
  };
}
