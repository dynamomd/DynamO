/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef GTKMM_EXAMPLE_TREEVIEW_WITHPOPUP_H
#define GTKMM_EXAMPLE_TREEVIEW_WITHPOPUP_H

#include <gtkmm.h>
#include "../extcode/xmlParser.h"
#include <iostream>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/lexical_cast.hpp>
namespace io = boost::iostreams;
#include "../base/constants.hpp"
#include "../base/is_exception.hpp"


class TagView : public Gtk::TreeView
{
public:
  TagView();
  virtual ~TagView();

  virtual void loadfile(std::string);

protected:
  // Override Signal handler:
  // Alternatively, use signal_button_press_event().connect_notify()
  virtual bool on_button_press_event(GdkEventButton *ev);

  //Signal handler for popup menu items:
  virtual void on_menu_file_popup_generic();

  virtual void save_correlator();

  virtual void grace_bulk_avg();

  virtual void grace_shear_avg();

  virtual void thermal_cond_avg();

  virtual void contactMap();

  virtual void vacf_avg();

  virtual void Histogram();

  template<int i>
  void VHistogram()
  {
    Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
    
    if (refSelection)
      {
	Gtk::TreeModel::iterator iter = refSelection->get_selected();
	if (iter)
	  {
	    int id = (*iter)[m_Columns.m_col_id];
	    XMLNode xBrowseNode = xMainNode.getChildNode(id), xhistogram;
	    if (xBrowseNode.hasChild("WeightHistogram"))
	      for (int j = 1; j <= i; j++)
		{
		  xhistogram = xBrowseNode.getChildNode("WeightHistogram",j);
		  io::filtering_ostream coutputFile;
		  coutputFile.push(io::file_sink(std::string("tmpfile.dat.f") + boost::lexical_cast<std::string>(j) ));
		  coutputFile << xhistogram.getText();
		}
	    else
	      if (xBrowseNode.hasChild("Histogram"))
	      for (int j = 1; j <= i; j++)
		{
		  xhistogram = xBrowseNode.getChildNode("Histogram",j);
		  io::filtering_ostream coutputFile;
		  coutputFile.push(io::file_sink(std::string("tmpfile.dat.f") + boost::lexical_cast<std::string>(j) ));
		  coutputFile << xhistogram.getText();
		}
	    std::string tmp = "xmgrace ";
	    for (int j = 1; j <= i; j++)	    
	      tmp += " tmpfile.dat.f" + j;
	    system((tmp + " &").c_str());
	  }
      }
  }

  template<const int i, const int j>
  void VHistogram()
  {
    bool multlevel = false;
    std::string Name;
    switch (j)
      {
      case (0):
	Name = "NematicOrderParameter";
	multlevel = true;
	break;
      case (1):
	Name = "GyrationRadii";
	multlevel = true;
      case (2):
	Name = "CubaticOrderParameter";
	multlevel = true;
	break;
      case (3):
	Name = "CBTSize";
	multlevel = false;
	break;
      };

    Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
    
    if (refSelection)
      {
	Gtk::TreeModel::iterator iter = refSelection->get_selected();
	if (iter)
	  {
	    int id = (*iter)[m_Columns.m_col_id];
	    XMLNode xBrowseNode, xhistogram;

	    if (multlevel)
	      xBrowseNode = xMainNode.getChildNode(id).getChildNode(0).getChildNode(Name.c_str());
	    else
	      xBrowseNode = xMainNode.getChildNode(id).getChildNode(Name.c_str());

	    std::cout << "\n\n\n" << xBrowseNode.getName() << "\n";

	    if (xBrowseNode.nChildNode("WeightHistogram"))
	      for (int j = 0; j < i; j++)
		{
		  std::cout << "\nFound WeightHistogram\n";
		  xhistogram = xBrowseNode.getChildNode("WeightHistogram",j);
		  io::filtering_ostream coutputFile;
		  coutputFile.push(io::file_sink(std::string("tmpfile.dat.f") + boost::lexical_cast<std::string>(j) ));
		  coutputFile << xhistogram.getText();
		}
	    else
	      if (xBrowseNode.nChildNode("Histogram"))
		for (int j = 0; j < i; j++)
		  {
		    std::cout << "\nFound Histogram\n";
		    xhistogram = xBrowseNode.getChildNode("Histogram",j);
		    io::filtering_ostream coutputFile;
		    coutputFile.push(io::file_sink(std::string("tmpfile.dat.f") + boost::lexical_cast<std::string>(j) ));
		    coutputFile << xhistogram.getText();
		  }

	    std::string tmp = "xmgrace ";
	    for (int j = 0; j < i; j++)	    
	      tmp = tmp + std::string(" tmpfile.dat.f") + boost::lexical_cast<std::string>(j);
	    system((tmp + " &").c_str());
	  }
      }
  }

  virtual void molGamma();
  virtual void sysGamma();

  //Tree model columns:
  class ModelColumns : public Gtk::TreeModel::ColumnRecord
  {
  public:

    ModelColumns()
    { add(m_col_id); add(m_col_name); add(m_tag_name); add(m_col_val); }

    Gtk::TreeModelColumn<unsigned int> m_col_id;
    Gtk::TreeModelColumn<Glib::ustring> m_col_name;
    Gtk::TreeModelColumn<Glib::ustring> m_tag_name;
    Gtk::TreeModelColumn<Iflt> m_col_val;
  };
  
  ModelColumns m_Columns;

  //The Tree model:
  Glib::RefPtr<Gtk::ListStore> m_refTreeModel;
  Gtk::Menu m_Menu_Popup;
  XMLNode xMainNode;

};

#endif //GTKMM_EXAMPLE_TREEVIEW_WITHPOPUP_H
