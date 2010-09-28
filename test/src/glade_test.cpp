#include <gtkmm.h>
#include <iostream>

int main(int argc, char **argv)
{
  Gtk::Main kit(argc, argv);

  Glib::RefPtr<Gtk::Builder> refXml 
    = Gtk::Builder::create_from_file("visualizer/visualizer.glade");

  kit.run();

  return 0;
}
