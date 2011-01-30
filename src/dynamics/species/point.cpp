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

#ifdef DYNAMO_visualizer
# include <gtkmm.h>
#endif

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <cstring>
#include "include.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../ranges/1range.hpp"
#include "../ranges/1RAll.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include <boost/tokenizer.hpp>

#ifdef DYNAMO_visualizer
# include <magnet/thread/mutex.hpp>
# include "../liouvillean/CompressionL.hpp"
# include <magnet/HSV.hpp>
# include <coil/RenderObj/Spheres.hpp>
#endif

void 
SpPoint::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
    
    try {
      mass = boost::lexical_cast<double>(XML.getAttribute("Mass"))
	* Sim->dynamics.units().unitMass();
      spName = XML.getAttribute("Name");
      intName = XML.getAttribute("IntName");
    } 
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in SpPoint";
      }

}

void 
SpPoint::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Mass") 
      << mass / Sim->dynamics.units().unitMass()
      << xml::attr("Name") << spName
      << xml::attr("IntName") << intName
      << xml::attr("Type") << "Point"
      << range;
}

void
SpPoint::initialise()
{ 
  if (IntPtr == NULL)
    M_throw() << "SpPoint missing a matching interaction";
}


#ifdef DYNAMO_visualizer
# include <magnet/gtk/numericEntry.hpp>

namespace { 
  class DynamoSphereRenderer: public RTSpheres
  {
  public:
    DynamoSphereRenderer(size_t N, std::string name, 
			 magnet::function::Delegate1<magnet::CL::CLGLState&, void> updateColorFunc):
      RTSpheres(N, name),
      _mode(COLOR_BY_ID),
      _updateColorFunc(updateColorFunc)
    {
      for (size_t i(0); i < 4; ++i)
	_colorStatic.s[i] = _colorFixed.s[i] = 255;
    }
    
    virtual void initGTK() 
    {
      _optList.reset(new Gtk::VBox);//The Vbox of options
      _singleColorMode.reset(new Gtk::RadioButton("Single Color"));
      _colorByIDMode.reset(new Gtk::RadioButton("Color by ID"));
      _colorBySpeedMode.reset(new Gtk::RadioButton("Color by Speed"));

      _RFixed.reset(new Gtk::SpinButton);
      _GFixed.reset(new Gtk::SpinButton);
      _BFixed.reset(new Gtk::SpinButton);
      _AFixed.reset(new Gtk::SpinButton);
      _characteristicV.reset(new Gtk::Entry);

      _colorIfStatic.reset(new Gtk::CheckButton("Color if Static/Sleeping"));
      _RStatic.reset(new Gtk::SpinButton);
      _GStatic.reset(new Gtk::SpinButton);
      _BStatic.reset(new Gtk::SpinButton);
      _AStatic.reset(new Gtk::SpinButton);


      {//Static color and RGBA boxes
	Gtk::HBox* box = manage(new Gtk::HBox);
	
	box->pack_start(*_colorIfStatic, true, true);_colorIfStatic->show();
	_colorIfStatic->set_active();

	Gtk::Label* label = manage(new Gtk::Label("RGBA"));
	box->pack_start(*label, false, false); label->show();
	
	_RStatic->set_increments(1.0, 1.0);
	_RStatic->set_range(0.0, 255.0);
	_RStatic->set_value(_colorStatic.s[0]);
	_GStatic->set_increments(1.0, 1.0);
	_GStatic->set_range(0.0, 255.0);
	_GStatic->set_value(_colorStatic.s[1]);
	_BStatic->set_increments(1.0, 1.0);
	_BStatic->set_range(0.0, 255.0);
	_BStatic->set_value(_colorStatic.s[2]);
	_AStatic->set_increments(1.0, 1.0);
	_AStatic->set_range(0.0, 255.0);
	_AStatic->set_value(_colorStatic.s[3]);      

	box->pack_start(*_RStatic, false, false); _RStatic->show();
	box->pack_start(*_GStatic, false, false); _GStatic->show();
	box->pack_start(*_BStatic, false, false); _BStatic->show();
	box->pack_start(*_AStatic, false, false); _AStatic->show();
		
	_optList->add(*box);
	box->show();

	//Horizontal seperator
	Gtk::HSeparator* line = manage(new Gtk::HSeparator);
	line->show();
	_optList->add(*line);
      }

      {//Single color and RGBA boxes
	Gtk::HBox* box = manage(new Gtk::HBox);
	
	box->pack_start(*_singleColorMode, true, true);_singleColorMode->show();
	if (_mode == SINGLE_COLOR) _singleColorMode->set_active();

	Gtk::Label* label = manage(new Gtk::Label("RGBA"));
	box->pack_start(*label, false, false); label->show();
	
	_RFixed->set_increments(1.0, 1.0);
	_RFixed->set_range(0.0, 255.0);
	_RFixed->set_value(_colorFixed.s[0]);
	
	_GFixed->set_increments(1.0, 1.0);
	_GFixed->set_range(0.0, 255.0);
	_GFixed->set_value(_colorFixed.s[1]);
	
	_BFixed->set_increments(1.0, 1.0);
	_BFixed->set_range(0.0, 255.0);
	_BFixed->set_value(_colorFixed.s[2]);
	
	_AFixed->set_increments(1.0, 1.0);
	_AFixed->set_range(0.0, 255.0);
	_AFixed->set_value(_colorFixed.s[3]);      

	box->pack_start(*_RFixed, false, false); _RFixed->show();
	box->pack_start(*_GFixed, false, false); _GFixed->show();
	box->pack_start(*_BFixed, false, false); _BFixed->show();
	box->pack_start(*_AFixed, false, false); _AFixed->show();
		
	_optList->add(*box);
	box->show();

	//Horizontal seperator
	Gtk::HSeparator* line = manage(new Gtk::HSeparator);
	line->show();
	_optList->add(*line);
      }

      Gtk::RadioButton::Group group = _singleColorMode->get_group();
      
      {//Color by ID
	_colorByIDMode->set_group(group);
	_colorByIDMode->show();
	if (_mode == COLOR_BY_ID) _colorByIDMode->set_active();

	_optList->add(*_colorByIDMode);
	
	//Horizontal seperator
	Gtk::HSeparator* line = manage(new Gtk::HSeparator);
	line->show();
	_optList->add(*line);
      }

      {//Color by Speed
	Gtk::HBox* box = manage(new Gtk::HBox);	
	
	_colorBySpeedMode->set_group(group);
	_colorBySpeedMode->show();
	if (_mode == COLOR_BY_SPEED) _colorBySpeedMode->set_active();

	box->pack_start(*_colorBySpeedMode, true, true);
	
	Gtk::Label* label = manage(new Gtk::Label("Max Speed"));
	box->pack_start(*label, false, false); label->show();

	box->pack_start(*_characteristicV, false, false);
	_characteristicV->show();
	_characteristicV->set_text("1.0");

	box->show();
	_optList->add(*box);
      }
      
      _singleColorMode->signal_toggled()
	  .connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));
      _RFixed->signal_value_changed()
	.connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));
      _GFixed->signal_value_changed()
	.connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));
      _BFixed->signal_value_changed()
	.connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));
      _AFixed->signal_value_changed()
	.connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));

      _colorByIDMode->signal_toggled()
	  .connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));

      _colorBySpeedMode->signal_toggled()
	  .connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));
      _characteristicV->signal_changed()
	.connect(sigc::bind<Gtk::Entry&>(&magnet::Gtk::forceNumericEntry, *_characteristicV));
      _characteristicV->signal_activate()
	.connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));


      _colorIfStatic->signal_toggled()
	  .connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));
      _RStatic->signal_value_changed()
	.connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));
      _GStatic->signal_value_changed()
	.connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));
      _BStatic->signal_value_changed()
	.connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));
      _AStatic->signal_value_changed()
	.connect(sigc::mem_fun(*this, &DynamoSphereRenderer::guiUpdate));

      _optList->show();

      guiUpdate();
    }

    virtual void showControls(Gtk::ScrolledWindow* win) 
    {
      win->remove();
      _optList->unparent();
      win->add(*_optList);
      win->show();
    }
    
    typedef enum {
      SINGLE_COLOR = 1,
      COLOR_BY_ID = 2,
      COLOR_BY_SPEED = 3
    } DrawMode;

    inline volatile const DrawMode& getDrawMode() { return _mode; }
    inline volatile const double& getScaleV() { return _scaleV; }
    inline volatile const cl_uchar4& getColorFixed() { return _colorFixed; }
    inline volatile const cl_uchar4& getColorStatic() { return _colorStatic; }
    inline volatile const bool& getRecolorOnUpdate() { return _recolorOnUpdate; }
    inline volatile bool getColorIfStatic() { return _colorStaticParticles; }
  protected:

    void guiUpdate()
    {
      std::string val = _characteristicV->get_text();
      if (val.empty()) {val = "1.0"; _characteristicV->set_text("1.0"); }

      _scaleV = boost::lexical_cast<double>(val);
      if (_scaleV == 0) {_scaleV = 1; _characteristicV->set_text("1.0"); }
      
      if (_singleColorMode->get_active())  { _mode = SINGLE_COLOR; _recolorOnUpdate = false; }
      if (_colorByIDMode->get_active())    { _mode = COLOR_BY_ID; _recolorOnUpdate = false; }
      if (_colorBySpeedMode->get_active()) { _mode = COLOR_BY_SPEED; _recolorOnUpdate = true; }
      _colorFixed.s[0] = _RFixed->get_value();
      _colorFixed.s[1] = _GFixed->get_value();
      _colorFixed.s[2] = _BFixed->get_value();
      _colorFixed.s[3] = _AFixed->get_value();

      _colorStatic.s[0] = _RStatic->get_value();
      _colorStatic.s[1] = _GStatic->get_value();
      _colorStatic.s[2] = _BStatic->get_value();
      _colorStatic.s[3] = _AStatic->get_value();
      
      _colorStaticParticles = _colorIfStatic->get_active();

      _systemQueue->queueTask(new magnet::function::Task1<void, magnet::CL::CLGLState&>
			      (_updateColorFunc, *_CLState));
    }

    std::auto_ptr<Gtk::VBox> _optList;
    std::auto_ptr<Gtk::CheckButton> _colorIfStatic;
    std::auto_ptr<Gtk::SpinButton> _RStatic;
    std::auto_ptr<Gtk::SpinButton> _GStatic;
    std::auto_ptr<Gtk::SpinButton> _BStatic;
    std::auto_ptr<Gtk::SpinButton> _AStatic;
    std::auto_ptr<Gtk::RadioButton> _singleColorMode;
    std::auto_ptr<Gtk::RadioButton> _colorByIDMode;
    std::auto_ptr<Gtk::RadioButton> _colorBySpeedMode;
    std::auto_ptr<Gtk::SpinButton> _RFixed;
    std::auto_ptr<Gtk::SpinButton> _GFixed;
    std::auto_ptr<Gtk::SpinButton> _BFixed;
    std::auto_ptr<Gtk::SpinButton> _AFixed;
    std::auto_ptr<Gtk::Entry> _characteristicV;
    
    volatile cl_uchar4 _colorFixed;
    volatile cl_uchar4 _colorStatic;
    volatile DrawMode _mode;
    volatile double _scaleV;
    volatile bool _recolorOnUpdate;
    volatile bool _colorStaticParticles;

    magnet::function::Delegate1<magnet::CL::CLGLState&, void> _updateColorFunc;
  };
}


magnet::thread::RefPtr<RenderObj>& 
SpPoint::getCoilRenderObj() const
{
  if (!_renderObj.isValid())
    {
      _renderObj = new DynamoSphereRenderer(range->size(), "Species: " + spName,
					    magnet::function::MakeDelegate(this, &SpPoint::updateColorObj));
      particleData.resize(range->size());
      particleColorData.resize(range->size());
    }

  return _renderObj;
}

void
SpPoint::updateRenderObj(magnet::CL::CLGLState& CLState) const
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  ///////////////////////POSITION DATA UPDATE
  //Check if the system is compressing and adjust the radius scaling factor
  float factor = 1;
  if (Sim->dynamics.liouvilleanTypeTest<LCompression>())
    factor = (1 + static_cast<const LCompression&>(Sim->dynamics.getLiouvillean()).getGrowthRate() * Sim->dSysTime);
 
  double diam = getIntPtr()->hardCoreDiam() * factor;
  
  size_t sphID(0);
  BOOST_FOREACH(unsigned long ID, *range)
    {
      Vector pos = Sim->particleList[ID].getPosition();
      
      Sim->dynamics.BCs().applyBC(pos);
      
      for (size_t i(0); i < NDIM; ++i)
	particleData[sphID].s[i] = pos[i];
      
      particleData[sphID++].w = diam * 0.5;
    }

  {
    CLState.getCommandQueue().enqueueWriteBuffer
      (static_cast<RTSpheres&>(*_renderObj).getSphereDataBuffer(),
       false, 0, range->size() * sizeof(cl_float4), &particleData[0]);
  }
  
  if (_renderObj.as<DynamoSphereRenderer>().getRecolorOnUpdate())
    updateColorObj(CLState);
}

void 
SpPoint::updateColorObj(magnet::CL::CLGLState& CLState) const
{
  bool colorIfSleeping = _renderObj.as<DynamoSphereRenderer>().getColorIfStatic();
  cl_uchar4 sleepcolor;
  for (size_t cc(0); cc < 4; ++cc)
    sleepcolor.s[cc] = _renderObj.as<DynamoSphereRenderer>().getColorStatic().s[cc];

  switch (_renderObj.as<DynamoSphereRenderer>().getDrawMode())
    {
    case DynamoSphereRenderer::SINGLE_COLOR:
      {
	cl_uchar4 color;
	for (size_t cc(0); cc < 4; ++cc)
	  color.s[cc] = _renderObj.as<DynamoSphereRenderer>().getColorFixed().s[cc];

	size_t sphID(0);
	BOOST_FOREACH(unsigned long ID, *range)
	  {
	    if (colorIfSleeping && !Sim->particleList[ID].testState(Particle::DYNAMIC))
	      for (size_t cc(0); cc < 4; ++cc)
		particleColorData[sphID].s[cc] = sleepcolor.s[cc];
	    else
	      for (size_t cc(0); cc < 4; ++cc)
		particleColorData[sphID].s[cc] = color.s[cc];

	    ++sphID;
	  }
	break;
      }
    case DynamoSphereRenderer::COLOR_BY_ID:
      {
	size_t np = range->size();

	size_t sphID(0);
	BOOST_FOREACH(unsigned long ID, *range)
	  {
	    if (colorIfSleeping && !Sim->particleList[ID].testState(Particle::DYNAMIC))
	      for (size_t cc(0); cc < 4; ++cc)
		particleColorData[sphID].s[cc] = sleepcolor.s[cc];
	    else
	      magnet::color::HSVtoRGB(particleColorData[sphID], ((float)(sphID)) / np);

	    ++sphID;
	  }
	break;
      }
    case DynamoSphereRenderer::COLOR_BY_SPEED:
      {
	double scaleV = _renderObj.as<DynamoSphereRenderer>().getScaleV();
	size_t sphID(0);
	BOOST_FOREACH(unsigned long ID, *range)
	  {
	    if (colorIfSleeping && !Sim->particleList[ID].testState(Particle::DYNAMIC))
	      for (size_t cc(0); cc < 4; ++cc)
		particleColorData[sphID].s[cc] = sleepcolor.s[cc];
	    else
	      {
		Vector vel = Sim->particleList[ID].getVelocity();
		magnet::color::HSVtoRGB(particleColorData[sphID], magnet::clamp(vel.nrm() / scaleV, 0.0, 1.0));
	      }
	    ++sphID;
	  }
	break;
      }
    default:
      M_throw() << "Not Implemented";
    }

  {
    CLState.getCommandQueue().enqueueWriteBuffer
      (static_cast<RTSpheres&>(*_renderObj).getColorDataBuffer(),
       false, 0, range->size() * sizeof(cl_uchar4), &particleColorData[0]);
  }
}
#endif
