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

namespace magnet { namespace GL { namespace objects { class CairoSurface; }}}

namespace coil {
  class RenderObjectsGtkTreeView;

  class RLight;

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
    enum RenderMode {
      COLOR = 1, //!< The object is to render its color information.
      DEPTH = 2, //!< The object is to render its depth information.
      SHADOW = 4, //!< This is a shadow calculation pass (typically
		  //!depth only).
      PICKING = 8, //!< This is an object picking pass
      DEFAULT = COLOR | DEPTH, //!< By default, color and depth
			       //!information should be rendered.
      SHADOW_PASS = SHADOW | DEPTH, //!< Shadow passes only need depth
				    //!information.
      PICKING_PASS = COLOR | PICKING //!< In a picking pass, we only
				     //!render flat shaded images
				     //!using unique colors to
				     //!identify objects selected by
				     //!the user.
    };

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
    
    /*! \brief Called when the RenderObject must be drawn in the
	OpenGL scene using deferred shading.
	
	\param cam The active camera for the render.
     */
    virtual void glRender(const magnet::GL::Camera& cam, RenderMode mode)
    {}
    
    /*! \brief Called when the RenderObject must be drawn in the
	OpenGL scene using forward rendering.
	
	\param fbo The target framebuffer object which the scene is
	being rendered to.
	
	\param cam The active camera for the render.
     */
    virtual void forwardRender(magnet::GL::FBO& fbo, 
			       const magnet::GL::Camera& cam,
			       std::vector<std::tr1::shared_ptr<RLight> >& lights,
			       GLfloat ambientLight,
			       RenderMode mode) 
    {}

    /*! \brief Called when the RenderObject should draw its 2D interface controls.
	\param camera The active camera for the render.
     */
    virtual void interfaceRender(const magnet::GL::Camera& camera, 
				 magnet::GL::objects::CairoSurface& cairo) {}


    /*! \brief The render phase of the picking render.
      
      The picking render determines the current object underneath the
      cursor by drawing every object in a unique color and sampling
      the pixel underneath the mouse. 
      
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

      Typically, this will just call \ref glRender() but with a unique
      color buffer generated in \ref initPicking() . However, if the
      object has special needs (like a custom shader), then extra
      logic will need to be implemented here.

      \param offset This number is the running sum of "pickable"
      objects. This value should be increased by the number of unique
      objects drawn in the \ref pickingRender() function before
      initPicking returns.

      \param cam The active camera for the render.
     */
    virtual void pickingRender(const magnet::GL::Camera& cam, 
			       const uint32_t offset) 
    {}

    /*! \brief The number of objects available for picking rendering.
      
      This should return 0 if no objects would be rendered in a call
      to \ref pickingRender (i.e., the object is invisible).
     */
    virtual uint32_t pickableObjectCount() { return 0; }


    /*! \brief this Render object has been selected/picked. Return a
        pointer to the actual object selected.
	
	By default the 

	\param objID An ID number of the object selected in the
	picking.  
	
	\param my_ptr A shared pointer of *this object.
     */
    virtual std::tr1::shared_ptr<RenderObj> 
    getPickedObject(uint32_t objID, const std::tr1::shared_ptr<RenderObj>& my_ptr)
    { return my_ptr; }

    
    /*! \brief Get the text to be displayed about the picked object.

	\param objID The object ID which was picked, returned from a
	picking render.
	
	\sa pickingRender()
     */
    virtual std::string getCursorText(uint32_t objID)
    {
      M_throw() << "This object is not pickable";
    }

    /*! \brief Get the object space coordinates of the picked part of
        this RenderObject.

	\param objID The object ID which was picked, returned from a
	picking render.
	
	\sa pickingRender()
     */
    virtual std::tr1::array<GLfloat, 4> getCursorPosition(uint32_t objID)
    {
      M_throw() << "This object is not pickable";
    }

    /*! \brief Used to notify the RenderObject that it has been
        dragged by the user.

	\param cursorPos The cursors current position in Object space
	(3D untransformed coordinates).
       
	\param objID The object ID which was dragged, returned from a
	picking render.

	\sa pickingRender()
     */
    virtual void dragCallback(Vector cursorPos, uint32_t objID) {}

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
    virtual Gtk::TreeModel::iterator addViewRows(RenderObjectsGtkTreeView& view, Gtk::TreeModel::iterator& iter);

    /*! \brief Return the icon used for the object in the render view.
     */
    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

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

    inline void setName(const std::string& name) { _name = name; }

    /*! \brief Returns the system queue.
     */    
    std::tr1::shared_ptr<magnet::thread::TaskQueue> getQueue() { return _systemQueue; }

    
    /*! \brief Called when the object should be deleted.
     */
    virtual void request_delete() {}

    virtual bool deletable() { return false; }

  protected:
    std::string _name;
    bool _visible;
    bool _shadowCasting;
    std::tr1::shared_ptr<magnet::thread::TaskQueue> _systemQueue;
  };

  class RenderObjectsGtkTreeView
  {
  public:
    void init(Gtk::TreeView* tree);    

    void buildRenderView();

    struct ModelColumns : Gtk::TreeModelColumnRecord
    {
      ModelColumns()
      { add(m_name); add(m_visible); add(m_shadowcasting); add(m_delete); add(m_obj); add(m_icon); }
      
      Gtk::TreeModelColumn<Glib::ustring> m_name;
      Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > m_visible;
      Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > m_shadowcasting;
      Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > m_delete;
      Gtk::TreeModelColumn<RenderObj*> m_obj;
      Gtk::TreeModelColumn<Glib::RefPtr<Gdk::Pixbuf> > m_icon;
    };
    
    std::auto_ptr<ModelColumns> _columns;
    Glib::RefPtr<Gtk::TreeStore> _store;
    Gtk::TreeView* _view;

    std::vector<std::tr1::shared_ptr<RenderObj> > _renderObjects;

  protected:
    

    int _obj_col, _visible_col, _shadow_col;

    bool button_press(GdkEventButton*);
    void name_edited(const Glib::ustring& path,const Glib::ustring& string);
  };
}
