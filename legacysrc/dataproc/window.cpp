/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <iostream>
#include "window.hpp"

MainWindow::MainWindow(std::string filename)
{
  m_TreeView.loadfile(filename);

  set_title("DYNAMO Tag Browser");
  set_border_width(5);
  set_default_size(400, 200);

  add(m_VBox);

  //Add the TreeView, inside a ScrolledWindow, with the button underneath:
  m_ScrolledWindow.add(m_TreeView);

  //Only show the scrollbars when they are necessary:
  m_ScrolledWindow.set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

  m_VBox.pack_start(m_ScrolledWindow);
  m_VBox.pack_start(m_ButtonBox, Gtk::PACK_SHRINK);

  m_ButtonBox.pack_start(m_Button_Quit, Gtk::PACK_SHRINK);
  m_ButtonBox.set_border_width(5);
  m_ButtonBox.set_layout(Gtk::BUTTONBOX_END);
  m_Button_Quit.signal_clicked().connect( sigc::mem_fun(*this,
              &MainWindow::on_button_quit) );

  show_all_children();

}
MainWindow::MainWindow()
: m_Button_Quit("Quit")
{
  //Find the file to load
  Gtk::FileChooserDialog dialog("Please choose an output file", Gtk::FILE_CHOOSER_ACTION_OPEN);
  dialog.set_transient_for(*this);
  set_default_size(800, 400);

  //Add response buttons the the dialog:
  dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
  dialog.add_button(Gtk::Stock::OPEN, Gtk::RESPONSE_OK);

  Gtk::FileFilter filter_compressed;
  filter_compressed.set_name("Compressed Output");
  filter_compressed.add_pattern("*.xml.bz2");
  dialog.add_filter(filter_compressed);
  Gtk::FileFilter filter_xml;
  filter_xml.set_name("unCompressed Output");
  filter_xml.add_pattern("*.xml");
  dialog.add_filter(filter_xml);

  //Show the dialog and wait for a user response:
  int result = dialog.run();

  //Handle the response:
  switch(result)
  {
    case(Gtk::RESPONSE_OK):
    {
      m_TreeView.loadfile(dialog.get_filename());
      break;
    }
    case(Gtk::RESPONSE_CANCEL):
    {
      exit(0);
      break;
    }
    default:
    {
      std::cout << "Unexpected button clicked." << std::endl;
      exit(0);
      break;
    }
  }

  set_title("DYNAMO Tag Browser");
  set_border_width(5);
  set_default_size(400, 200);

  add(m_VBox);

  //Add the TreeView, inside a ScrolledWindow, with the button underneath:
  m_ScrolledWindow.add(m_TreeView);

  //Only show the scrollbars when they are necessary:
  m_ScrolledWindow.set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

  m_VBox.pack_start(m_ScrolledWindow);
  m_VBox.pack_start(m_ButtonBox, Gtk::PACK_SHRINK);

  m_ButtonBox.pack_start(m_Button_Quit, Gtk::PACK_SHRINK);
  m_ButtonBox.set_border_width(5);
  m_ButtonBox.set_layout(Gtk::BUTTONBOX_END);
  m_Button_Quit.signal_clicked().connect( sigc::mem_fun(*this,
              &MainWindow::on_button_quit) );

  show_all_children();
}

MainWindow::~MainWindow()
{
}

void MainWindow::on_button_quit()
{
  hide();
}
