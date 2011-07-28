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
#include <magnet/GL/context.hpp>
#include <magnet/thread/taskQueue.hpp>
#include <magnet/GL/camera.hpp>
#include <magnet/GL/FBO.hpp>
#include <tr1/memory>

namespace Gtk { class ScrolledWindow; }

namespace coil {
  class RenderObjectsGtkTreeView;

  /*! \brief The base class for all renderable objects in the system.
    
    Any object which represents a 3D object or a 2D interface object
    must derive from this class. This class provides the abstract
    interface for the main \ref CLGLWindow window to interact with the
    object, requesting it to render its interface, the object and so
    on whenever it is required.
  */
  class RenderObj
  {
  public:
    /*! \brief Default constructor which just sets the name of the
      object.
    */
    RenderObj(std::string name): _name(name), _visible(true), _shadowCasting(true) {}
  
    /* \brief Initialises the object and any OpenCL, OpenGL or GTK
       resources it contains.
       
       \param systemQueue A reference to the task queue executed by
       the (simulation) thread which is providing the rendered
       data. This is to allow callbacks to the (simulation) thread
       when user-generated interface events occur.
     */
    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue) 
    { _systemQueue = systemQueue; }

    /*! \brief Release any OpenCL, OpenGL and GTK resources held by
        the object.
     */
    virtual void deinit() {}
    
    /*! \brief Called to perform any OpenCL calls before the OpenGL
        part of the render loop.
	
	\param cam The camera used for the next \ref glRender call.
     */
    virtual void clTick(const magnet::GL::Camera& cam) = 0;

    /*! \brief Called when the RenderObject must be drawn in the OpenGL scene.
	
	\param fbo The target framebuffer object which the scene is
	being rendered to.
	
	\param cam The active camera for the render.
     */
    virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam) = 0;

    
    /*! \brief Called when the RenderObject should draw its 2D interface controls.
	\param camera The active camera for the render.
     */
    virtual void interfaceRender(const magnet::GL::Camera& camera) {}


    /*! \brief The initialisation phase of the picking render.
      
      The picking render determines the current object underneath the
      cursor by drawing every object in a unique color and sampling
      the pixel underneath the mouse. This function is used to setup
      the unique coloring of the object ready for a picking render
      phase.
      
      An \ref offset value is passed into this function to allow the
      render object to determine unique colors for its own objects. If
      the RenderObj represents 12 "pickable" objects it must increase
      this offset value by 12 before returning.

      As the colors of objects are specified in coil using 4 8-bit
      values, we can directly convert a 32bit cl_uint type into a
      cl_uchar4 object to generate a unique color from an object
      ID. \ref offset represents the number of pickable objects that
      will be rendered before pickingRender is called on this
      object. Thus \ref offset is an offset to be applied to the
      unique colors generated for this object.

      \param offset This number is the running sum of "pickable"
      objects. This value should be increased by the number of unique
      objects drawn in the \ref pickingRender() function before
      initPicking returns.

      \sa pickingRender
     */
    virtual void initPicking(cl_uint& offset) {}
    
    /*! \brief The render phase of the picking render.

      Typically, this will just call \ref glRender() but with a unique
      color buffer generated in \ref initPicking() . However, if the
      object has special needs (like a custom shader), then extra
      logic will need to be implemented here.

      \param fbo The target framebuffer object which the scene is
      being rendered to.
      
      \param cam The active camera for the render.

      \sa initPicking()
     */
    virtual void pickingRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam) {}


    /*! \brief The cleanup and callback phase of the picking render.
      
      In this stage of the picking process, the objects can clean up
      any resources allocated in \ref initPicking() and also can
      determine if they were the object selected. 
      
      \param offset This value is identical to the value passed in
      \ref initPicking(), and must also be incremented by the number
      of pickable objects rendered by this RenderObj before this
      function returns..

      \param val This is the unique ID of the object selected. If this
      lies in the range of [offset, offset + "no. of pickable
      objects"), then an object was picked from this RenderObj, and
      some appropriate callback should be performed.
      
      \sa initPicking()
     */
    virtual void finishPicking(cl_uint& offset, const cl_uint val) {}

    
    /*! \brief Callback for when the RenderObj is to make its Gtk
        controls visible.
	
	\param win The ScrolledWindow to which this RenderObj must add
	its controls.
     */
    virtual void showControls(Gtk::ScrolledWindow* win) {}


    /*! \brief Callback for when the RenderObj is to create its
	entries in the RenderObjectsGtkTreeView.
	
	This function should add a line in the TreeView interface used
	to list the available RenderObj class instances. If this class
	has "children" entries, it is responsible to add these entries
	to the view.

	The iterator corresponding to the row of this RenderObj is
	returned to allow derived types to inherit the base classes
	code for generating entries in the TreeView, while still being
	able to add children to that generated row.

	\param view The TreeView used to list the available render
	objects.
	
	\returns An iterator pointing to the row in the TreeView
	generated by this RenderObj.
     */
    virtual inline Gtk::TreeModel::iterator addViewRows(RenderObjectsGtkTreeView& view);

    /*! \brief Sets the object's visibility.

      If this object is not visible, it will not have \ref glRender()
      called during the next render phase and it will not cast
      shadows.
     */
    inline void setVisible(bool val = true) { _visible = val; }
    inline bool visible() const { return _visible; }


    /*! \brief Sets the object's shadow casting.

      If this object is set to not cast shadows, it will not have \ref
      glRender() called during the next light render phase. 
     */
    inline void setShadowCasting(bool val = true) { _shadowCasting = val; }

    /*! \brief A test whether the RenderObj will cast a shadow.
      
      Any RenderObj which can never cast a shadow should override this
      function to always return false.
     */
    virtual bool shadowCasting() const { return _shadowCasting; }

    /*! \brief Returns the name of the RenderObj.
     */
    inline const std::string& getName() const { return _name; }

    /*! \brief Returns the system queue.
     */    
    std::tr1::shared_ptr<magnet::thread::TaskQueue> getQueue() { return _systemQueue; }

  protected:
    std::string _name;
    bool _visible;
    bool _shadowCasting;
    std::tr1::shared_ptr<magnet::thread::TaskQueue> _systemQueue;
  };

  class RenderObjectsGtkTreeView
  {
  public:
    void init(Gtk::TreeView* tree)
    { 
      _columns.reset(new ModelColumns); 
      _view = tree;
      _store = Gtk::TreeStore::create(*_columns);
      _view->set_model(_store);

      { //The cell renderer 
	Gtk::CellRendererToggle* renderer =
	  Gtk::manage(new Gtk::CellRendererToggle());
	int visible_col = _view->append_column("Visible", *renderer);
	Gtk::TreeViewColumn* column = _view->get_column(visible_col -1);
	if (column)
	  column->add_attribute(renderer->property_active(), _columns->m_visible);
	renderer->signal_toggled().connect(sigc::mem_fun(*this, &RenderObjectsGtkTreeView::visibleToggled));
      }

      { //The cell renderer 
	Gtk::CellRendererToggle* renderer =
	  Gtk::manage(new Gtk::CellRendererToggle());
	int shadow_col = _view->append_column("Shadow Casting", *renderer);
	Gtk::TreeViewColumn* column = _view->get_column(shadow_col -1);
	if (column)
	  column->add_attribute(renderer->property_active(), _columns->m_shadowcasting);
	renderer->signal_toggled().connect(sigc::mem_fun(*this, &RenderObjectsGtkTreeView::shadowingToggled));
      }

      _view->append_column("Object Name", _columns->m_name);

    }
    
    struct ModelColumns : Gtk::TreeModelColumnRecord
    {
      ModelColumns()
      { add(m_name); add(m_visible); add(m_shadowcasting); add(m_obj);}
      
      Gtk::TreeModelColumn<Glib::ustring> m_name;
      Gtk::TreeModelColumn<bool> m_visible;
      Gtk::TreeModelColumn<bool> m_shadowcasting;
      Gtk::TreeModelColumn<RenderObj*> m_obj;
    };
    
    std::auto_ptr<ModelColumns> _columns;
    Glib::RefPtr<Gtk::TreeStore> _store;
    Gtk::TreeView* _view;
    
  protected:
    
    void visibleToggled(const Glib::ustring& path_string)
    {
      Gtk::TreeModel::iterator iter 
	= _store->get_iter(path_string);

      bool visible = !(*iter)[_columns->m_visible];
      (*iter)[_columns->m_visible] = visible;

      RenderObj* obj = (*iter)[_columns->m_obj];
      obj->setVisible(visible);
    }

    void shadowingToggled(const Glib::ustring& path_string)
    {
      Gtk::TreeModel::iterator iter 
	= _store->get_iter(path_string);

      bool casting = !(*iter)[_columns->m_shadowcasting];
      (*iter)[_columns->m_shadowcasting] = casting;

      RenderObj* obj = (*iter)[_columns->m_obj];
      obj->setShadowCasting(casting);
    }

  };

  Gtk::TreeModel::iterator RenderObj::addViewRows(RenderObjectsGtkTreeView& view)
  {
    Gtk::TreeModel::iterator iter = view._store->append();
    (*iter)[view._columns->m_name] = getName();
    (*iter)[view._columns->m_visible] = visible();
    (*iter)[view._columns->m_shadowcasting] = shadowCasting();
    (*iter)[view._columns->m_obj] = this;
    return iter;
  }
}
