#!/bin/bash

echo "/*Generated Image header file for coil*/" > images.hpp
echo "#include <gtkmm.h>" >> images.hpp
echo "namespace coil { namespace images {" >> images.hpp

echo "/*Generated Image function file for coil*/" > images.cpp
echo "#include <coil/images/images.hpp>" >> images.cpp
echo "namespace coil { namespace images {" >> images.cpp
for image in $(ls *.svg)
do
    basename=$(basename $image)
    basename=${basename%.*}
    gdk-pixbuf-csource --static --name $basename"_data" $image >> images.cpp
    echo "Glib::RefPtr<Gdk::Pixbuf> $basename();" >> images.hpp
    echo "Glib::RefPtr<Gdk::Pixbuf> $basename()" >> images.cpp
    echo "{ return Gdk::Pixbuf::create_from_inline(sizeof("$basename"_data), "$basename"_data); }" >> images.cpp
done

echo "}}" >> images.hpp
echo "}}" >> images.cpp
