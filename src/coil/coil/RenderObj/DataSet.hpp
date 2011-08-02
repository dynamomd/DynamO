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
#include <coil/RenderObj/RenderObj.hpp>
#include <magnet/GL/buffer.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <boost/lexical_cast.hpp>
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
   * An important implementation note on Attributes is that they're
   * initialised on access. This is to facilitate the main thread
   * adding attributes after the GL threads initialisation phase.
   */
  class Attribute
  {
  public:
    enum AttributeType { 
      INTENSIVE = 1, //!< Intensive property (e.g., Temperature, density)
      EXTENSIVE = 2, //!< Extensive property (e.g., mass, momentum)
      COORDINATE = 3 //!< A special attribute which specifies the location of the attribute.
    };

    inline Attribute(size_t N, AttributeType type = EXTENSIVE, size_t components = 1):
      _hostData(N * components),
      _components(components),
      _type(type),
      _references(0)
    {}
    
    inline Attribute() {}

    /*! \brief Releases the OpenGL resources of this object.
     */
    inline void deinit() { _glData.deinit(); }

    /*! \brief Returns the GL buffer associated with the Attribute
     * data.
     */
    inline magnet::GL::Buffer<GLfloat>& getBuffer() { return _glData; }

    /** @name The host code interface. */
    /**@{*/
    
    /*! \brief Returns a reference to the host-cache of the Attribute
     * data.
     *
     * The Attribute data may be directly updated by the host program,
     * but flagNewData() must be called for the update to take effect.
     */
    inline std::vector<GLfloat>& getData() { return _hostData; }
    
    /*! \brief Marks that the data in the buffer has been updated, and
     * should be uploaded to the GL system.
     *
     * This function just inserts a callback in the GL system to
     * reinitialise the Attribute.
     */
    inline void flagNewData()
    { _glData.getContext().queueTask(magnet::function::Task::makeTask(&Attribute::initGLData, this)); }

    /*! \brief Test if the attribute is in use and should be
     * updated. 
     */
    inline bool active() const { return _references; }

    /*! \brief Returns the number of elements.
     */
    inline size_t size() const { return _hostData.size() / _components; }

    inline size_t getNComponents() const { return _components; }

    inline AttributeType getType() const { return _type; }

    /**@}*/

    inline void bindAttribute(size_t attrnum, bool normalise = false)
    {
      //Initialise on demand
      if (!_glData.size()) initGLData();
      _glData.attachToAttribute(attrnum, _components, 1, normalise); 
    }

  protected:
    void initGLData() 
    { _glData.init(_hostData); }

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

  class DataSet; 

  class DataSetChild: public RenderObj
  {
  public:
    inline DataSetChild(std::string name, DataSet& ds): RenderObj(name), _ds(ds) {}

    virtual Gtk::TreeModel::iterator addViewRows(RenderObjectsGtkTreeView& view, Gtk::TreeModel::iterator&) = 0;
    
  protected:

    DataSet& _ds;
  };
  
  /*! \brief A container class for a collection of \ref Attribute
   * instances forming a dataset, and any active filters/glyphs or any
   * other type derived from \ref DataSetChild.
   */
  class DataSet: public RenderObj, public std::map<std::string, Attribute>
  {
  public:
    DataSet(std::string name, size_t N): 
      RenderObj(name), 
      _context(NULL),
      _N(N) {}
    
    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue);

    virtual void deinit();

    /** @name The host code interface. */
    /**@{*/

    /*! \brief Add an Attribute to the DataSet.
     *
     * \param name The name of the Attribute.
     * \param type The type of the attribute.
     * \param components The number of components per value of the attribute.
     */
    void addAttribute(std::string name, Attribute::AttributeType type, size_t components);

    /*! \brief Looks up an attribute by its name.
     */
    inline Attribute& operator[](const std::string& name)
    {
      iterator iPtr = find(name);
      if (iPtr == end())
	M_throw() << "No attribute named " << name << " in Data set";
      
      return iPtr->second;
    }

    inline size_t size() const { return _N; }

    /**@}*/
        
    virtual void clTick(const magnet::GL::Camera& cam) 
    {
      for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	   iPtr != _children.end(); ++iPtr)
	(*iPtr)->clTick(cam);
    }

    virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam, RenderMode mode = DEFAULT)
    {
      for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	   iPtr != _children.end(); ++iPtr)
	if ((*iPtr)->visible() && (!(mode & SHADOW) || (*iPtr)->shadowCasting()))
	  (*iPtr)->glRender(fbo, cam, mode);
    }
    
    virtual Gtk::TreeModel::iterator addViewRows(RenderObjectsGtkTreeView& view)
    {
      _view = &view;
      _iter = RenderObj::addViewRows(view);
      
      for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	   iPtr != _children.end(); ++iPtr)
	(*iPtr)->addViewRows(view, _iter);

      return _iter;
    }

    virtual void showControls(Gtk::ScrolledWindow* win);

  protected:
    /*! \brief An iterator to this DataSet's row in the Render object
      treeview.
     */
    Gtk::TreeModel::iterator _iter;
    RenderObjectsGtkTreeView* _view;
    magnet::GL::Context* volatile _context;
    std::auto_ptr<Gtk::VBox> _gtkOptList;
    size_t _N;
    std::vector<std::tr1::shared_ptr<DataSetChild> > _children;

    void initGtk();
    void rebuildGui();

    void addGlyphs();

    struct ModelColumns : Gtk::TreeModelColumnRecord
    {
      ModelColumns()
      { add(name); add(components); add(type); }
      
      Gtk::TreeModelColumn<Glib::ustring> name;
      Gtk::TreeModelColumn<size_t> components;
      Gtk::TreeModelColumn<Attribute::AttributeType> type;
    };
    
    std::auto_ptr<ModelColumns> _attrcolumns;
    Glib::RefPtr<Gtk::TreeStore> _attrtreestore;
    std::auto_ptr<Gtk::TreeView> _attrview;
  };

  class AttributeSelector : public Gtk::HBox
  {
  public:
    AttributeSelector():
      _components(0)
    {
      //Label
      _label.show();
      pack_start(_label, false, false, 5);
      _context = &(magnet::GL::Context::getContext());
      //combo box
      _model = Gtk::ListStore::create(_modelColumns);
      _comboBox.set_model(_model);      
      _comboBox.pack_start(_modelColumns.m_name);
      _comboBox.show();
      pack_start(_comboBox, false, false, 5);
      
      _singleValueLabel.show();
      _singleValueLabel.set_text("Value:");
      _singleValueLabel.set_alignment(1.0, 0.5);
      pack_start(_singleValueLabel, true, true, 5);
      for (size_t i(0); i < 4; ++i)
	{
	  pack_start(_scalarvalues[i], false, false, 0);
	  _scalarvalues[i].signal_changed()
	    .connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, _scalarvalues[i]));
	  _scalarvalues[i].set_text("1.0");
	  _scalarvalues[i].set_max_length(0);
	  _scalarvalues[i].set_width_chars(5);	  
	}

      show();
    }
    
    void buildEntries(std::string name, DataSet& ds, Attribute::AttributeType typemask, size_t components)
    {
      _label.set_text(name);
      _model->clear();
      _components = components;

      for (size_t i(0); i < components; ++i)
	_scalarvalues[i].show();

      for (size_t i(components); i < 4; ++i)
	_scalarvalues[i].hide();

      for (DataSet::iterator iPtr = ds.begin();
	   iPtr != ds.end(); ++iPtr)
	if ((iPtr->second.getType() & typemask) && (iPtr->second.getNComponents() ==  _components))
	  {
	    Gtk::TreeModel::Row row = *(_model->append());
	    row[_modelColumns.m_name] = iPtr->first;
	    row[_modelColumns.m_ptr] = &(iPtr->second);
	  }

      {
	Gtk::TreeModel::Row row = *(_model->append());
	row[_modelColumns.m_name] = "Single Value";
	row[_modelColumns.m_ptr] = NULL;
      }
      
      _comboBox.set_active(0);
    }
    
    struct ModelColumns : Gtk::TreeModelColumnRecord
    {
      ModelColumns()
      { add(m_name); add(m_ptr); }
      
      Gtk::TreeModelColumn<Glib::ustring> m_name;
      Gtk::TreeModelColumn<Attribute*> m_ptr;
    };

    void bindAttribute(size_t attrnum, bool normalise = false)
    {
      Gtk::TreeModel::iterator iter = _comboBox.get_active();

      if (!iter) { setConstantAttribute(attrnum); return; }
      
      Attribute* ptr = (*iter)[_modelColumns.m_ptr];
      if (!ptr) { setConstantAttribute(attrnum); return; }

      ptr->bindAttribute(attrnum, normalise);
    }

    ModelColumns _modelColumns;
    Gtk::ComboBox _comboBox;
    Gtk::Label _label;
    Gtk::Label _singleValueLabel;
    Glib::RefPtr<Gtk::ListStore> _model;
    Gtk::Entry _scalarvalues[4];
    
  protected:
    magnet::GL::Context* _context;
    
    size_t _components;
    
    void setConstantAttribute(size_t attr)
    {
      _context->disableAttributeArray(attr);

      GLfloat val[4] = {1,1,1,1};

      for (size_t i(0); i < 4; ++i)
	try {
	  val[i] = boost::lexical_cast<double>(_scalarvalues[i].get_text());
	} catch (...) {}
      
      _context->setAttribute(attr, val[0], val[1], val[2], val[3]);
    }
  };

}
