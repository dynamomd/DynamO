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
      _attenuation(0.5), _specularExponent(96),
      _specularFactor(1) 
    {
      std::tr1::array<GLfloat, 3> tmp = {{1.0, 1.0, 1.0}};
      _color = tmp;
    }
  
    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue);
    virtual void deinit();

    virtual void glRender(magnet::GL::FBO& fbo, 
			  const magnet::GL::Camera& cam, 
			  RenderMode mode);

    virtual void clTick(const magnet::GL::Camera&) {}

    virtual void showControls(Gtk::ScrolledWindow* win);

    float getIntensity() const { return _intensity; }
    float getAttenuation() const { return _attenuation; }
    float getSpecularExponent() const { return _specularExponent; }
    float getSpecularFactor() const { return _specularFactor; }

    const std::tr1::array<GLfloat, 3>& getColor() { return _color; }

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
      magnet::math::Vector vec = getEyeLocationObjSpace();
      std::tr1::array<GLfloat, 4> lightPos = {{vec[0], vec[1], vec[2], 1.0}};
      std::tr1::array<GLfloat, 4> lightPos_eyespace
	= camera.getViewMatrix() * lightPos;
      return magnet::math::Vector(lightPos_eyespace[0], lightPos_eyespace[1], 
				  lightPos_eyespace[2]);
    }

    virtual uint32_t pickableObjectCount()
    { return visible(); }

    virtual void pickingRender(magnet::GL::FBO& fbo,
			       const magnet::GL::Camera& cam, 
			       const uint32_t offset);
    
    virtual std::string getCursorText(uint32_t objID)
    { return _name; }

    virtual std::tr1::array<GLfloat, 4> getCursorPosition(uint32_t objID)
    {
      magnet::math::Vector loc = getEyeLocationObjSpace();
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
    std::auto_ptr<Gtk::Entry> _attenuationEntry;
    std::auto_ptr<Gtk::Entry> _specularExponentEntry;
    std::auto_ptr<Gtk::Entry> _specularFactorEntry;
    
    float _intensity, _attenuation, _specularExponent, _specularFactor;
    std::tr1::array<GLfloat, 3> _color;

    magnet::GL::shader::SphereShader<true> _sphereShader;
    magnet::GL::Buffer<GLfloat> _glposition;
    magnet::GL::Context::ContextPtr _context;
  };
}
