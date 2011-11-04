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

#include <magnet/color/transferFunction.hpp>
#include <magnet/color/HSV.hpp>
#include <magnet/function/delegate.hpp>
#include <gtkmm/drawingarea.h>
#include <gdk/gdkkeysyms.h>
#include <gtkmm/colorselection.h>
#include <tr1/array>

namespace magnet {
  namespace gtk {
    class TransferFunction : public ::Gtk::DrawingArea
    {
      typedef magnet::color::TransferFunction::Knot Knot;

      Gdk::Color ConvertKnotToGdk(const Knot& knot)
      {
	Gdk::Color retval;
	float color[4];
	magnet::color::HSVtoRGB<float>(color,knot._h, knot._s, knot._v, knot._a);
	retval.set_red(color[0] * G_MAXUSHORT);
	retval.set_green(color[1] * G_MAXUSHORT);
	retval.set_blue(color[2] * G_MAXUSHORT);
	return retval;
      }

      Knot ConvertGdkToKnot(const Gdk::Color& col, const guint16 alpha, double x)
      {
	float color[4];
	
	magnet::color::RGBtoHSV(color,
				float(col.get_red()) / G_MAXUSHORT,
				float(col.get_green()) / G_MAXUSHORT,
				float(col.get_blue()) / G_MAXUSHORT,
				float(alpha) / G_MAXUSHORT);

	Knot knot(x, color[0], color[1], color[2], color[3]);
	return knot;
      }

    public:
      TransferFunction(function::Delegate0<void> updated):
	_updatedCallback(updated),
	_grid_line_width(1), _selectedNode(-1), _dragMode(false)
      {
	set_events(::Gdk::KEY_PRESS_MASK | ::Gdk::EXPOSURE_MASK | ::Gdk::BUTTON_PRESS_MASK
		   | ::Gdk::BUTTON_RELEASE_MASK | ::Gdk::POINTER_MOTION_MASK);
	set_flags(::Gtk::CAN_FOCUS);

	_transferFunction.addKnot(0.0, 1, 0, 1, 0.0);
	_transferFunction.addKnot(0.5, 1, 0, 1, 0.5);
	_transferFunction.addKnot(1.0, 1, 0, 1, 1.0);
      }

      std::vector<uint8_t> getColorMap()
      { return _transferFunction.getIndexRGBMap(); }

      std::vector<float>& getHistogram() { return _histogram; }

      virtual ~TransferFunction() {}

    protected:
      void forceRedraw()
      {
	Glib::RefPtr<Gdk::Window> win = get_window();
	if (win)
	  {
	    Gdk::Rectangle r(0, 0, get_allocation().get_width(),
			     get_allocation().get_height());
	    win->invalidate_rect(r, false);
	  }
      }

      std::pair<float,float> getGraphPosition(GdkEventButton* event)
      {
	const ::Gtk::Allocation& allocation = get_allocation();
	//Calculate the mouse coordinates relative to the top left of the screen
	float x = event->x - allocation.get_x();
	float y = event->y - allocation.get_y();

	//Flip the y direction so its at the bottom
	return std::pair<float, float>(x, y);
      }

      int getClickedKnot(GdkEventButton* event)
      {
	std::pair<float,float> pos(getGraphPosition(event));
      
	typedef magnet::color::TransferFunction::const_iterator it;

	for (it iPtr = _transferFunction.begin(); iPtr != _transferFunction.end(); ++iPtr)
	  {
	    std::pair<float,float> pointpos = to_graph_transform(iPtr->_x, iPtr->_a);
	    pointpos.first -= pos.first;
	    pointpos.first *= pointpos.first;
	    pointpos.second -= pos.second;
	    pointpos.second *= pointpos.second;

	    if (pointpos.first + pointpos.second  <= 0.25f * getPointSize() * getPointSize())
	      return iPtr - _transferFunction.begin();
	  }

	return -1;
      }

      virtual bool on_key_press_event(GdkEventKey* event) 
      {
#ifndef GDK_KEY_Delete
# define GDK_KEY_Delete GDK_Delete 
#endif
	if ((event->type == GDK_KEY_PRESS)
	    && (event->keyval == GDK_KEY_Delete)
	    && (_selectedNode >= 0)
	    && (_transferFunction.size() > 2))
	  {
	    _transferFunction.eraseKnot(_transferFunction.begin() + _selectedNode);
	    _selectedNode = -1;
	    forceRedraw();
	    _updatedCallback();
	  }
	return ::Gtk::DrawingArea::on_key_press_event(event);
      }

      virtual bool on_button_press_event(GdkEventButton* event) 
      {
	grab_focus();

	if (event->button == 1)//Left mouse click
	  switch (event->type)
	    {
	    case GDK_BUTTON_PRESS: //Single click
	      {
		_selectedNode = getClickedKnot(event);
		if (_selectedNode >= 0) _dragMode = true;
		forceRedraw();
	      }
	      break;
	    case GDK_2BUTTON_PRESS: //Double click
	      if ((_selectedNode = getClickedKnot(event)) == -1)
		{//Add a new node!
		  const ::Gtk::Allocation& allocation = get_allocation();
		  std::pair<double, double> newPlace 
		    = from_graph_transform(event->x - allocation.get_x(), 
					   event->y - allocation.get_y());
		  _transferFunction.addInterpolatedKnot(newPlace.first);
		  _updatedCallback();
		  forceRedraw();
		}
	      else
		{//Color an existing node!
		  _dragMode = false;
		  ::Gtk::ColorSelectionDialog select("Node Color Selection");

		  typedef magnet::color::TransferFunction::Knot Knot;

		  magnet::color::TransferFunction::const_iterator 
		    iPtr = _transferFunction.begin() + _selectedNode;
		  
		  Knot knot = *iPtr;
		  
		  select.get_color_selection()->set_current_color(ConvertKnotToGdk(knot));
		  select.get_color_selection()->set_current_alpha(knot._a * G_MAXUSHORT);
		  select.get_color_selection()->set_has_opacity_control(true);
		    
		  switch(select.run())
		    {
		    case ::Gtk::RESPONSE_OK:
		      {
			_transferFunction
			  .setKnot(iPtr, ConvertGdkToKnot(select.get_color_selection()->get_current_color(), 
							  select.get_color_selection()->get_current_alpha(), 
							  knot._x));
			_updatedCallback();
			forceRedraw();
		      }
		      break;
		    case ::Gtk::RESPONSE_CANCEL:
		      break;
		    default:
		      M_throw() << "Unexpected return value!";
		    }
		}
	    default:
	      break;
	    }

	return ::Gtk::DrawingArea::on_button_press_event(event);
      }
  
      virtual bool on_button_release_event(GdkEventButton* event) 
      { 
	if ((event->button == 1)//Left mouse click
	    && (event->type == GDK_BUTTON_RELEASE))
	  _dragMode = false;

	return ::Gtk::DrawingArea::on_button_release_event(event); 
      }
  
      virtual bool on_motion_notify_event(GdkEventMotion* event)
      { 
	const ::Gtk::Allocation& allocation = get_allocation();
	if (_dragMode && (_selectedNode >= 0))
	  {
	    std::pair<double, double> newPlace 
	      = from_graph_transform(event->x - allocation.get_x(), 
				     event->y - allocation.get_y());

	    magnet::color::TransferFunction::const_iterator 
	      iPtr = _transferFunction.begin() + _selectedNode;

	    _transferFunction.setKnot
	      (iPtr, magnet::color::TransferFunction::Knot
	       (newPlace.first, iPtr->_h, iPtr->_s, iPtr->_v, newPlace.second));
	
	    forceRedraw();
	    _updatedCallback();
	  }

	return ::Gtk::DrawingArea::on_motion_notify_event(event); 
      }

      inline double getPointSize() const
      { return 15 * _grid_line_width; }

      std::pair<double, double> to_graph_transform(double x, double y)
      {
	const ::Gtk::Allocation& allocation = get_allocation();
	return std::pair<double, double>(x * (allocation.get_width() - getPointSize()) 
					 + 0.5 * getPointSize(), 
					 (1 - y) * (allocation.get_height() - getPointSize())
					 + 0.5 * getPointSize());
      }

      std::pair<double, double> from_graph_transform(double x, double y)
      {
	const ::Gtk::Allocation& allocation = get_allocation();

	double xdraw((x - 0.5 * getPointSize()) / (allocation.get_width() - getPointSize()));
	double ydraw(1 - (y - 0.5 * getPointSize()) / (allocation.get_height() - getPointSize()));
	return std::pair<double, double>(xdraw, ydraw);
      }

      void graph_move_to(Cairo::RefPtr<Cairo::Context>& cr, double x, double y)
      {
	const std::pair<double, double>& pos(to_graph_transform(x, y));
	cr->move_to(pos.first, pos.second);
      }

      void graph_line_to(Cairo::RefPtr<Cairo::Context>& cr, double x, double y)
      {
	const std::pair<double, double>& pos(to_graph_transform(x, y));
	cr->line_to(pos.first, pos.second);
      }

      void graph_rectangle(Cairo::RefPtr<Cairo::Context>& cr, double x, double y,
			   double width, double height)
      {
	const ::Gtk::Allocation& allocation = get_allocation();
	const std::pair<double, double>& pos(to_graph_transform(x, y));
	cr->rectangle(pos.first, pos.second, 
		      width * (allocation.get_width() - getPointSize()), 
		      -height * (allocation.get_height() - getPointSize()));
      }

      bool on_expose_event(GdkEventExpose* event)
      {
	// This is where we draw on the window
	Glib::RefPtr<Gdk::Window> window = get_window();
	if(window)
	  {
	    Cairo::RefPtr<Cairo::Context> cr = window->create_cairo_context();

	    const ::Gtk::Allocation& allocation = get_allocation();
	
	    //Scale so that x is in the range [0,1] and that y goes from 0 (bottom) to 1 (top)	
	    if(event)
	      {
		// clip to the area indicated by the expose event so that we only
		// redraw the portion of the window that needs to be redrawn
		cr->rectangle(event->area.x, event->area.y,
			      event->area.width, event->area.height);
		cr->clip();
	      }

	    //Draw background colors
	    cr->save();
	    {
	      //First build the gradient 
	      Cairo::RefPtr<Cairo::LinearGradient> grad
		= Cairo::LinearGradient::create(0.5 * getPointSize(), 0,
						allocation.get_width() 
						- 0.5 * getPointSize(), 0);
	      for (size_t i(1); i < 256; ++i)
		{
		  std::tr1::array<float,4> val = _transferFunction.getValue(i / 255.0);
		  
		  float color[4];
		  magnet::color::HSVtoRGB(color, val[0], val[1], val[2], val[3]);

		  grad->add_color_stop_rgba(i / 255.0, color[0], color[1], color[2], color[3]);
		}
	      
	      cr->set_source(grad);
	      graph_rectangle(cr, 0, 0, 1, 1);
	      cr->fill();
	    }
	    cr->restore();

	    //Draw Grid
	    cr->set_source_rgba(0, 0, 0, 1);
	    cr->set_line_width(_grid_line_width);
	    //horizontal lines
	    for (size_t i(0); i < 5; ++i)
	      {
		graph_move_to(cr, 0.0f, i / 4.0f);
		graph_line_to(cr, 1.0f, i / 4.0f);
		cr->stroke();
	      }

	    //Vertical lines
	    for (size_t i(0); i < 5; ++i)
	      {
		graph_move_to(cr, i / 4.0f, 0.0f);
		graph_line_to(cr, i / 4.0f, 1.0f);
		cr->stroke();
	      }

	    { //draw the curve of the graph
	      cr->set_source_rgba(0, 0, 0, 1);
	      cr->set_line_width(5 * _grid_line_width);

	      //Initial point
	      graph_move_to(cr, 0, _transferFunction.getValue(0)[3]);
	      
	      for (size_t i(1); i < 256; ++i)
		graph_line_to(cr, i / 255.0, magnet::clamp(_transferFunction.getValue(i / 255.0)[3], 0.0f, 1.0f));
	  
	      cr->stroke();
	    }

	    if (_histogram.size() == 256)
	      { //Draw the histogram of the data set
		cr->set_source_rgba(0.2, 0.2, 0.2, 1);
		cr->set_line_width(2 * _grid_line_width);
		//Initial point
		graph_move_to(cr, 0, _histogram[0]);
		for (size_t i(1); i < 256; ++i)
		  graph_line_to(cr, i / 255.0, _histogram[i]);
		cr->stroke();
	      }

	    
	    { // draw the nodes of the graph
	      typedef magnet::color::TransferFunction::const_iterator it;

	      for (it iPtr = _transferFunction.begin(); iPtr != _transferFunction.end(); ++iPtr)
		{
		  cr->set_source_rgba(0, 0, 0, 1);
		  cr->set_line_width(5 * _grid_line_width);
		  const std::pair<double, double>& pos(to_graph_transform(iPtr->_x, iPtr->_a));
		  cr->arc(pos.first, pos.second,  getPointSize() / 2, 0, 2 * M_PI);

		  float color[4];
		  magnet::color::HSVtoRGB<float>(color, iPtr->_h, iPtr->_s, iPtr->_v, 1);


		  cr->set_source_rgba(color[0], color[1], color[2], 1.0);
		  cr->fill_preserve();

		  if (iPtr - _transferFunction.begin() == _selectedNode)
		    {
		      cr->set_source_rgba(1, 1, 1, 1); cr->stroke();
		      cr->arc(pos.first, pos.second,  
			      (getPointSize() + 4 * _grid_line_width) /  2, 0, 2 * M_PI);
		      cr->set_source_rgba(0, 0, 0, 1); cr->stroke();
		    }
		  else
		    { cr->set_source_rgba(0, 0, 0, 1); cr->stroke(); }
		}
	    }
	  }

	return true;
      }

      function::Delegate0<void> _updatedCallback;
      color::TransferFunction _transferFunction;
      std::vector<float> _histogram;
      double _grid_line_width;
      int _selectedNode;
      bool _dragMode;
    };
  }
}
