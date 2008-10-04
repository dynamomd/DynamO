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

#include "treeview.hpp"

#include <stdio.h>

TagView::TagView()
{

}

void
TagView::loadfile(std::string fileName)
{
  if (std::string(fileName.end()-4, fileName.end()) == ".xml")
    {
      if (!boost::filesystem::exists(fileName))
	D_throw() << "Could not open XML file";
      
      std::cout << "Uncompressed XML input file " << fileName << " loading";
      xMainNode=XMLNode::openFileHelper(fileName.c_str(), "DYNAMOconfig");
    }
  else if (std::string(fileName.end()-8, fileName.end()) == ".xml.bz2")
    {
      if (!boost::filesystem::exists(fileName))
	D_throw() << "Could not open XML file";
      
      io::filtering_istream inputFile;
      inputFile.push(io::bzip2_decompressor());
      inputFile.push(io::file_source(fileName));
      //Copy file to a string
      std::string line, fileString;
          
      std::cout << "Bzip compressed XML input file found\nDecompressing file "
		<< fileName;
      
      while(getline(inputFile,line)) 
	{
	  fileString.append(line);
	  fileString.append("\n");
	}
      XMLNode tmpNode = XMLNode::parseString(fileString.c_str());
      xMainNode = tmpNode.getChildNode("OutputData");
    }
  else
    D_throw() << "Unrecognised extension for input file";
  
  std::cout << "Parsing XML file";
  //XMLNode xSubNode= xMainNode.getChildNode("OutputData");
  XMLNode xBrowseNode;

  //Create the Tree model:
  m_refTreeModel = Gtk::ListStore::create(m_Columns);
  set_model(m_refTreeModel);

  //Fill the TreeView's model
  for (int i = 0; i < xMainNode.nChildNode(); i++)
    {
      xBrowseNode = xMainNode.getChildNode(i);
      Gtk::TreeModel::Row row = *(m_refTreeModel->append());
      row[m_Columns.m_col_id] = i;
      row[m_Columns.m_tag_name] = xBrowseNode.getName();
      if (xBrowseNode.isAttributeSet("Name"))
	row[m_Columns.m_col_name] = xBrowseNode.getAttribute("Name");
      if (xBrowseNode.isAttributeSet("Val"))
	row[m_Columns.m_col_val] = atof(xBrowseNode.getAttribute("Val"));
      if (xBrowseNode.isAttributeSet("val"))
	row[m_Columns.m_col_val] = atof(xBrowseNode.getAttribute("val"));
    }

  //Add the TreeView's view columns:
  append_column("ID", m_Columns.m_col_id);
  append_column("Type", m_Columns.m_tag_name);
  append_column("Name", m_Columns.m_col_name);
  append_column("Value", m_Columns.m_col_val);
}

TagView::~TagView()
{
}

bool TagView::on_button_press_event(GdkEventButton* event)
{
  //Call base class, to allow normal handling,
  //such as allowing the row to be selected by the right-click:
  bool return_value = TreeView::on_button_press_event(event);

  //Then do our custom stuff:
  if( (event->type == GDK_BUTTON_PRESS) && (event->button == 3) )
  {
    Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();

    if(refSelection)
      {	
	//Fill popup menu:
	{
	  Gtk::Menu::MenuList& menulist = m_Menu_Popup.items();
	  menulist.clear();
	  if ((*(refSelection->get_selected()))[m_Columns.m_tag_name] == "EinsteinCorrelator")
	    {
	      menulist.push_back( Gtk::Menu_Helpers::MenuElem("_Save to file",
							      sigc::mem_fun(*this, &TagView::save_correlator) ) );

	      if ((*(refSelection->get_selected()))[m_Columns.m_col_name] == "Viscosity")
		{
		  menulist.push_back( Gtk::Menu_Helpers::MenuElem("xmGrace Bulk Average",
								  sigc::mem_fun(*this, &TagView::grace_bulk_avg) ) );

		  menulist.push_back( Gtk::Menu_Helpers::MenuElem("xmGrace Shear Average",
								  sigc::mem_fun(*this, &TagView::grace_shear_avg) ) );
		}

	      if ((*(refSelection->get_selected()))[m_Columns.m_col_name] == "ThermalConductivity")
		{
		  menulist.push_back( Gtk::Menu_Helpers::MenuElem("xmGrace Average",
								  sigc::mem_fun(*this, &TagView::thermal_cond_avg) ) );
		}

	    }
	  else if ((*(refSelection->get_selected()))[m_Columns.m_tag_name] == "Correlator")
	    {
	      if ((*(refSelection->get_selected()))[m_Columns.m_col_name] == "VACF")
		{
		  menulist.push_back( Gtk::Menu_Helpers::MenuElem("xmGrace Average",
								  sigc::mem_fun(*this, &TagView::vacf_avg) ) );
		}
	    }
	  else if ((*(refSelection->get_selected()))[m_Columns.m_tag_name] == "ChainTorsion")
	    {
	      menulist.push_back( Gtk::Menu_Helpers::MenuElem("xmGrace molGamma",
								  sigc::mem_fun(*this, &TagView::molGamma) ) );	      
	      menulist.push_back( Gtk::Menu_Helpers::MenuElem("xmGrace sysGamma",
								  sigc::mem_fun(*this, &TagView::sysGamma) ) );	      
	    }
	  else if ((*(refSelection->get_selected()))[m_Columns.m_tag_name] == "EnergyHist")
	    {
	      menulist.push_back( Gtk::Menu_Helpers::MenuElem("xmGrace histogram",
							      sigc::mem_fun(*this, &TagView::Histogram) ) );	      
	    }
	  else if ((*(refSelection->get_selected()))[m_Columns.m_tag_name] == "ChainGyration")
	    {
	      menulist.push_back( Gtk::Menu_Helpers::MenuElem("Radius squared",
							      sigc::mem_fun(*this, &TagView::VHistogram<3, 1>) ) );	      
	      menulist.push_back( Gtk::Menu_Helpers::MenuElem("Nematic Order Parameter",
							      sigc::mem_fun(*this, &TagView::VHistogram<3, 0>) ) );	      
	      menulist.push_back( Gtk::Menu_Helpers::MenuElem("Cubatic Order Parameters",
							      sigc::mem_fun(*this, &TagView::VHistogram<5, 2>) ) );	      
	    }
	  else if ((*(refSelection->get_selected()))[m_Columns.m_tag_name] == "boundedQstats")
	    {
	      menulist.push_back( Gtk::Menu_Helpers::MenuElem("CB-Tree Size Histogram",
							      sigc::mem_fun(*this, &TagView::VHistogram<1, 3>) ) );	      
	    }
	  else if ((*(refSelection->get_selected()))[m_Columns.m_tag_name] == "ContactMap")
	    {
	      menulist.push_back( Gtk::Menu_Helpers::MenuElem("Contact Map",
							      sigc::mem_fun(*this, &TagView::contactMap) ) );	      
	    }

	  if (menulist.empty())
	    menulist.push_back(Gtk::Menu_Helpers::MenuElem("Nothing"));

	    
	}
	
	m_Menu_Popup.popup(event->button, event->time);
      }
  
  }

  return return_value;
}

void 
TagView::Histogram()
{
  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
  {
    Gtk::TreeModel::iterator iter = refSelection->get_selected();
    if(iter)
      {
	    int id = (*iter)[m_Columns.m_col_id];
	    XMLNode xBrowseNode = xMainNode.getChildNode(id), xhistogram;
	    if (xBrowseNode.hasChild("WeightHistogram"))
	      xhistogram = xBrowseNode.getChildNode("WeightHistogram");
	    else
	      if (xBrowseNode.hasChild("Histogram"))
		xhistogram = xBrowseNode.getChildNode("Histogram");
		  

	{
	  io::filtering_ostream coutputFile;
	  coutputFile.push(io::file_sink("tmpfile.dat"));
	  coutputFile << xhistogram.getText();
	}
	
	system("cat tmpfile.dat | xmgrace -pipe &");
      }
  }
  
}

void 
TagView::molGamma()
{
  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
  {
    Gtk::TreeModel::iterator iter = refSelection->get_selected();
    if(iter)
      {
	int id = (*iter)[m_Columns.m_col_id];
	XMLNode xBrowseNode = xMainNode.getChildNode(id).getChildNode("HelixPolymer").getChildNode("MolecularHistogram").getChildNode("Histogram");
	
	{
	  io::filtering_ostream coutputFile;
	  coutputFile.push(io::file_sink("tmpfile.dat"));
	  coutputFile << xBrowseNode.getText();
	}
	
	system("cat tmpfile.dat | xmgrace -pipe &");
      }
  }
  
}

void 
TagView::sysGamma()
{
  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
  {
    Gtk::TreeModel::iterator iter = refSelection->get_selected();
    if(iter)
      {
	int id = (*iter)[m_Columns.m_col_id];
	XMLNode xBrowseNode = xMainNode.getChildNode(id).getChildNode("HelixPolymer").getChildNode("SystemHistogram").getChildNode("Histogram");
	
	{
	  io::filtering_ostream coutputFile;
	  coutputFile.push(io::file_sink("tmpfile.dat"));
	  coutputFile << xBrowseNode.getText();
	}
	
	system("cat tmpfile.dat | xmgrace -pipe &");
      }
  }
}

void 
TagView::contactMap()
{
  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
  {
    Gtk::TreeModel::iterator iter = refSelection->get_selected();
    if(iter)
      {
	int id = (*iter)[m_Columns.m_col_id];
	XMLNode xBrowseNode = xMainNode.getChildNode(id).getChildNode(0);
	
	{
	  io::filtering_ostream coutputFile;
	  coutputFile.push(io::file_sink("tmpfile.dat"));
	  coutputFile << xBrowseNode.getText();
	}
	
	system("echo -e \"set pm3d map\\nset xrange [-0.5 to 19.5]\\nset yrange [-0.5 to 19.5]\\nsplot 'tmpfile.dat'\\npause 500000\" | gnuplot &");
      }
  }
}

void 
TagView::grace_bulk_avg()
{
  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
  {
    Gtk::TreeModel::iterator iter = refSelection->get_selected();
    if(iter)
    {
      int id = (*iter)[m_Columns.m_col_id];
      XMLNode xBrowseNode = xMainNode.getChildNode(id);
      
      {
	io::filtering_ostream coutputFile;
	coutputFile.push(io::file_sink("tmpfile.dat"));
	coutputFile << xBrowseNode.getText();
      }
      
      system("cat tmpfile.dat | gawk '{print $1,($2+$6+$10)/3.0}' | xmgrace -pipe &");
    }
  }

}

void 
TagView::vacf_avg()
{
  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
  {
    Gtk::TreeModel::iterator iter = refSelection->get_selected();
    if(iter)
    {
      int id = (*iter)[m_Columns.m_col_id];
      XMLNode xBrowseNode = xMainNode.getChildNode(id);
      
      {
	io::filtering_ostream coutputFile;
	coutputFile.push(io::file_sink("tmpfile.dat"));
	coutputFile << xBrowseNode.getText();
      }
      
      system("cat tmpfile.dat | gawk '{print $1,($2+$3+$4)/3.0}' | xmgrace -pipe &");
    }
  }

}

void 
TagView::thermal_cond_avg()
{
  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
    {
      Gtk::TreeModel::iterator iter = refSelection->get_selected();
      if(iter)
	{
	  int id = (*iter)[m_Columns.m_col_id];
	  XMLNode xBrowseNode = xMainNode.getChildNode(id);
	  
	  {
	    io::filtering_ostream coutputFile;
	    coutputFile.push(io::file_sink("tmpfile.dat"));
	    coutputFile << xBrowseNode.getText();
	  }
	  
	  system("cat tmpfile.dat | gawk '{print $1,($2+$3+$4)/3.0}' | xmgrace -pipe &");
	}
    } 
}

void 
TagView::grace_shear_avg()
{
  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
  {
    Gtk::TreeModel::iterator iter = refSelection->get_selected();
    if(iter)
    {
      int id = (*iter)[m_Columns.m_col_id];
      XMLNode xBrowseNode = xMainNode.getChildNode(id);
      
      {
	io::filtering_ostream coutputFile;
	coutputFile.push(io::file_sink("tmpfile.dat"));
	coutputFile << xBrowseNode.getText();
      }
      
      system("cat tmpfile.dat | gawk '{print $1,($3+$4+$5+$7+$8+$9)/6.0}' | xmgrace -pipe &");
    }
  }
}

void 
TagView::save_correlator()
{
  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
  {
    Gtk::TreeModel::iterator iter = refSelection->get_selected();
    if(iter)
    {
      Gtk::FileChooserDialog dialog("Please choose an output file", Gtk::FILE_CHOOSER_ACTION_SAVE);
      //dialog.set_transient_for(*this);
      //Add response buttons the the dialog:
      dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
      dialog.add_button(Gtk::Stock::SAVE, Gtk::RESPONSE_OK);
      
      Gtk::FileFilter filter_compressed;
      filter_compressed.set_name("Data");
      filter_compressed.add_pattern(".dat");
      dialog.add_filter(filter_compressed);
      
      //Show the dialog and wait for a user response:
      int result = dialog.run();
      
      //Handle the response:
      switch(result)
	{
	case(Gtk::RESPONSE_OK):
	  {
	    io::filtering_ostream coutputFile;
	    coutputFile.push(io::file_sink(dialog.get_filename()));
	    int id = (*iter)[m_Columns.m_col_id];
	    XMLNode xBrowseNode = xMainNode.getChildNode(id);
      
	    coutputFile << xBrowseNode.getText();
	    break;
	  }
	case(Gtk::RESPONSE_CANCEL):
	  {
	    break;
	  }
	default:
	  {
	    std::cout << "Unexpected button clicked." << std::endl;
	    break;
	  }
	}
    }
  }  
}

void TagView::on_menu_file_popup_generic()
{
  std::cout << "A popup menu item was selected." << std::endl;

  Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
  if(refSelection)
  {
    Gtk::TreeModel::iterator iter = refSelection->get_selected();
    if(iter)
    {
      int id = (*iter)[m_Columns.m_col_id];
      std::cout << "  Selected ID=" << id << std::endl;
    }
  }
}
