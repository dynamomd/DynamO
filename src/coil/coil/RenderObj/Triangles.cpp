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
#include "Triangles.hpp"
#include <iostream>
#include <set>

namespace coil {
  RTriangles::RTriangles(magnet::GL::Context::ContextPtr context, std::string name):
    RenderObj(context, name),
    _RenderMode(TRIANGLES),
    _triangleComponents(3)
  {}

  RTriangles::~RTriangles()
  { deinit(); }

  void 
  RTriangles::init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  { 
    RenderObj::init(systemQueue); 
    initGTK(); 
    _renderShader.build(); 
    _renderVSMShader.build(); 
    _initialised = true;
  }

  void 
  RTriangles::glRender(const magnet::GL::Camera& cam, RenderMode mode, const uint32_t offset)
  {
    if (!_visible) return;

    magnet::GL::Buffer<GLubyte> _pickingColorBuff;
    if (mode & RenderObj::PICKING)
      {
	size_t N = (_posBuff.size() / _triangleComponents);
	std::vector<GLubyte> vertexColors;
	vertexColors.reserve(N * 4);
	
	for (size_t i(offset); i < N + offset; ++i)
	  {
	    vertexColors.push_back((i >>  0) & 0xFF);
	    vertexColors.push_back((i >>  8) & 0xFF);
	    vertexColors.push_back((i >> 16) & 0xFF);
	    vertexColors.push_back((i >> 24) & 0xFF);
	  }
	
	_pickingColorBuff.init(vertexColors, 4, magnet::GL::buffer_usage::STREAM_DRAW);
	_pickingColorBuff.attachToColor();
      }
    else
      _colBuff.attachToColor();

    using namespace magnet::GL::shader::detail;
    Shader& shader = (mode == RenderObj::SHADOW) ? static_cast<Shader&>(_renderVSMShader) : static_cast<Shader&>(_renderShader);
    shader.attach();
    shader["ProjectionMatrix"] = cam.getProjectionMatrix();
    shader["ViewMatrix"] = cam.getViewMatrix();

    _posBuff.getContext()->cleanupAttributeArrays();
    if (_normBuff.size()) _normBuff.attachToNormal();
    _posBuff.attachToVertex();
  
    switch (_RenderMode)
      {
      case TRIANGLES:
	_elementBuff.drawElements(magnet::GL::element_type::TRIANGLES);
	break;
      case LINES:
	_specialElementBuff.drawElements(magnet::GL::element_type::LINES);
	break;
      case POINTS:
	_specialElementBuff.drawElements(magnet::GL::element_type::POINTS);
	break;
      }
    shader.detach();
  }

  void 
  RTriangles::setGLColors(const std::vector<GLubyte>& VertexColor)
  {
    if (!VertexColor.size())
      throw std::runtime_error("VertexColor.size() == 0!");

    if (_posBuff.size())
      if ((VertexColor.size() / 4) != (_posBuff.size() / _triangleComponents))
	throw std::runtime_error("VertexColor.size() / 4 != posBuffSize/ comps");

    _colBuff.init(VertexColor, 4, magnet::GL::buffer_usage::STREAM_DRAW);
  }

  void 
  RTriangles::setGLPositions(const std::vector<float>& VertexPos)
  {
    if (!VertexPos.size())
      throw std::runtime_error("VertexPos.size() == 0!");
  
    if (VertexPos.size() % _triangleComponents)
      throw std::runtime_error("VertexPos.size() is not a multiple of 3!");

    const size_t N = VertexPos.size() / _triangleComponents;

    if (_colBuff.size())
      if (_colBuff.size() != N)
	throw std::runtime_error("VertexPos.size()/3 != colBuffSize/4 ");
  
    if (_normBuff.size())
      if (_normBuff.size() != VertexPos.size())
	throw std::runtime_error("VertexPos.size() != normBuffSize!");

    _posBuff.init(VertexPos, _triangleComponents, magnet::GL::buffer_usage::STREAM_DRAW);
  }

  void 
  RTriangles::setGLNormals(const std::vector<float>& VertexNormals)
  {
    if (!VertexNormals.size())
      throw std::runtime_error("VertexNormals.size() == 0!");

    if (VertexNormals.size() % 3)
      throw std::runtime_error("VertexNormals.size() is not a multiple of 3!");

    const size_t posSize = _posBuff.size() / _triangleComponents;
    if (_posBuff.size())
      if (VertexNormals.size() != 3 * posSize)
	throw std::runtime_error("VertexNormals.size() != posBuffsize!");

    _normBuff.init(VertexNormals, 3);
  }

  void 
  RTriangles::setGLElements(const std::vector<GLuint>& Elements)
  {
    if (!Elements.size())
      throw std::runtime_error("Elements.size() == 0!");

    if (Elements.size() % 3) 
      throw std::runtime_error("Elements.size() is not a multiple of 3!");

    _elementBuff.init(Elements, 3);
  }

  void 
  RTriangles::deinit()
  {
    _colBuff.deinit();
    _posBuff.deinit();
    _normBuff.deinit();
    _elementBuff.deinit();
    _specialElementBuff.deinit();
    _renderShader.deinit();
    _renderVSMShader.deinit();
  }

  void 
  RTriangles::initGTK()
  {
    _gtkOptList.reset(new Gtk::VBox);//The Vbox of options

    _gtkTriangleRender.reset(new Gtk::RadioButton("Solid"));
    _gtkLineRender.reset(new Gtk::RadioButton("Wireframe"));
    _gtkPointRender.reset(new Gtk::RadioButton("Vertex Points"));
  
    {
      Gtk::RadioButton::Group group = _gtkTriangleRender->get_group();
      _gtkLineRender->set_group(group);
      _gtkPointRender->set_group(group);

      _gtkTriangleRender->set_active();
    
      Gtk::HBox* box = manage(new Gtk::HBox);
      box->pack_start(*_gtkTriangleRender, true, true);
      box->pack_start(*_gtkLineRender, true, true);
      box->pack_start(*_gtkPointRender, true, true);
    
      _gtkTriangleRender->show();
      _gtkLineRender->show();
      _gtkPointRender->show();
      box->show();
      
      _gtkOptList->add(*box);
    }
    
    _gtkOptList->show();
    
    _gtkTriangleRender->signal_toggled()
      .connect(sigc::mem_fun(*this, &RTriangles::guiUpdate));
    _gtkLineRender->signal_toggled()
      .connect(sigc::mem_fun(*this, &RTriangles::guiUpdate));
    _gtkPointRender->signal_toggled()
      .connect(sigc::mem_fun(*this, &RTriangles::guiUpdate));

    guiUpdate();
  }

  void 
  RTriangles::showControls(Gtk::ScrolledWindow* win)
  {
    win->remove();
    _gtkOptList->unparent();
    win->add(*_gtkOptList);
    win->show();
  }

  void 
  RTriangles::guiUpdate()
  {
    RenderModeType rmode = TRIANGLES;
    if (_gtkLineRender->get_active()) rmode = LINES;
    if (_gtkPointRender->get_active()) rmode = POINTS;
  
    setRenderMode(rmode);
  }

  void 
  RTriangles::setRenderMode(RenderModeType rm)
  {
    if (rm != _RenderMode)
      {
	_specialElementBuff.deinit();

	//If we're in one of the special modes, build the special
	//element buffer
	switch (rm)
	  {
	  case LINES:
	    {
	      //Here, we find all the unique edges and build the
	      //corresponding element buffer
	      typedef std::pair<int,int> SetKey;
	      std::set<SetKey> edges;

	      size_t size = _elementBuff.size();

	      GLuint* elements =  _elementBuff.map();

	      for (size_t t(0); t < size; t += 3)
		{
		  edges.insert(SetKey(elements[t+0], elements[t+1]));
		  edges.insert(SetKey(elements[t+1], elements[t+2]));
		  edges.insert(SetKey(elements[t+0], elements[t+2]));
		}
	    
	      std::vector<GLuint> line_elements;
	      line_elements.reserve(edges.size() * 2);
	    
	      for (std::set<SetKey>::const_iterator iPtr = edges.begin();
		   iPtr != edges.end(); ++iPtr)
		{
		  line_elements.push_back(iPtr->first);
		  line_elements.push_back(iPtr->second);
		}

	      _specialElementBuff.init(line_elements, 2);

	      _elementBuff.unmap();
	      break;
	    }
	  case POINTS:
	    {
	      std::vector<GLuint> point_elements;
	      point_elements.resize(_posBuff.size() / 3);
	    
	      for (size_t i(0); i < point_elements.size(); ++i)
		point_elements[i] = i;

	      _specialElementBuff.init(point_elements, 1);
	      break;
	    }
	  default:
	    break;
	  }
      }
    _RenderMode = rm;
  }
}
