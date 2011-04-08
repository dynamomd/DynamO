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

#pragma once

#include <magnet/detail/rapidXML/rapidxml.hpp>
#include <magnet/exception.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/copy.hpp>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

namespace magnet {
  //!Namespace enclosing the XML tools included in magnet.
  namespace xml {
    namespace detail {
      //! Generates an error string including the full XML path of the passed node.
      //! \param error A custom error message to incorporate into the whole message.
      //! \param _node A valid node whose path is written in the error message.
      //! \return A error message string.
      inline std::string xmlerror(const std::string error, rapidxml::xml_node<> *_node)
      {
	std::vector<std::pair<std::string, size_t> > pathTree;
	while (_node != 0)
	  {
	    size_t index(0);
	    { rapidxml::xml_node<>* sibling = _node->previous_sibling(_node->name());
	      while (sibling) { sibling = sibling->previous_sibling(); ++index; }
	    }

	    pathTree.push_back(std::pair<std::string, size_t>(_node->name(), index));
	    _node = _node->parent();
	  }
	std::ostringstream os;
	os << "XML error: " << error << "\n"
	   << "XML Path: ";

	for (std::vector<std::pair<std::string, size_t> >::const_reverse_iterator 
	       iPtr = pathTree.rbegin(); iPtr != pathTree.rend(); ++iPtr)
	  os << "/" << iPtr->first << "[" << iPtr->second << "]";

	os << "/INVALID";
	return os.str();
      }
    }

    //! Represents an Attribute of an XML Document.
    class Attribute {
    public:
      //! Converts the attributes value to a type.
      template<class T> inline T as() const 
      { 
	if (!valid()) M_throw() << detail::xmlerror("Missing attribute being converted", _parent);
	return boost::lexical_cast<T>(_attr->value()); 
      }

      //! Converts the attributes value to a type.  This version of
      //! as() will return the passed value if the Attribute is
      //! invalid. 
      //! \param defaultval The default value to be returned if the Attribute is invalid.
      //! \sa as()
      template<class T> inline T as(T defaultval) const 
      { 
	if (!valid()) return defaultval;
	return boost::lexical_cast<T>(_attr->value()); 
      }

      //! Conversion operator returning the value of the attribute.
      inline operator const char*() const 
      { 
	if (!valid()) M_throw() << detail::xmlerror("Missing attribute being read", _parent);
	return _attr->value(); 
      }

      //! Test if this Attribute is valid and has a value.
      inline bool valid() const { return _attr != NULL; }

    private:
      friend class Node;
      Attribute();
      inline Attribute(rapidxml::xml_attribute<>* attr,
		       rapidxml::xml_node<>* parent):_attr(attr), _parent(parent) {}
      
      rapidxml::xml_attribute<> *_attr;
      rapidxml::xml_node<> *_parent;
    };

    //! Represents a Node of an XML Document.
    class Node {
    public:
      //! Fetches the first Attribute with a given name. Test if the
      //! Attribute is Attribute::valid() before using it.
      //! \param name The name of the attribute to return.
      template<class T> inline Attribute getAttribute(T name) const 
      { 
	if (!valid()) M_throw() << detail::xmlerror("Missing node's attribute being accessed", _parent);
	return Attribute(_node->first_attribute(name), _node); 
      }

      //! Fetches the first Node with a given name. Test if the
      //! Node is Node::valid() before using it.
      //! \param name The name of the Node to return.
      template<class T> inline Node getNode(T name) const 
      { 
	if (!valid()) M_throw() << detail::xmlerror("Cannot fetch sub node of invalid node", _parent);
	return Node(_node->first_node(name), _node); 
      }

      //! Returns the value of the Node.
      inline operator const char*() const 
      { 
	if (!valid()) M_throw() << detail::xmlerror("Cannot read invalid node", _parent);
	return _node->value(); 
      }

      //! Test if the Node is valid.
      inline bool valid() const { return _node != NULL; }

      //! Replace this Node with the next Node in the parent Node with the same name.
      inline void operator++()
      { 
	if (!valid()) M_throw() << detail::xmlerror("Cannot increment invalid node", _parent);
	_node = _node->next_sibling(_node->name()); 
      }
      //! Replace this Node with the previous Node in the parent Node with the same name.
      inline void operator--() 
      { 
	if (!valid()) M_throw() << detail::xmlerror("Cannot decrement invalid node", _parent);
	_node = _node->previous_sibling(_node->name()); 
      }

      //! Constructor to build a Node from rapidxml::xml_node data.
      inline Node(rapidxml::xml_node<>* node,
		  rapidxml::xml_node<>* parent):_node(node), _parent(parent) {}

   private:
      rapidxml::xml_node<> *_node;
      rapidxml::xml_node<> *_parent;
    };

    template<>
    inline Attribute Node::getAttribute<const std::string&>(const std::string& name) const
    { 
      if (!valid()) M_throw() << detail::xmlerror("Missing node's attribute being accessed", _parent);
      return Attribute(_node->first_attribute(name.c_str()), _node); 
    }

    template<>
    inline Node Node::getNode<const std::string&>(const std::string& name) const
    { 
      if (!valid()) M_throw() << detail::xmlerror("Cannot fetch sub node of invalid node", _parent);
      return Node(_node->first_node(name.c_str()), _node); 
    }


    //! A class which represents a whole XML Document, including
    //! storage. 
    //!
    //! This class thinly wraps the rapidXML parser, providing a
    //! simple type-casting interface and automatic error handling,
    //! similar to Boost's property_tree. Unlike property_tree this
    //! class is space-efficient and does not copy the XML data.  This
    //! class must outlive any Node or Attribute generated from it.
    class Document {
    public:
      //! Construct an XML Document from a file.
      //! \param fileName Path to the file to load.
      inline Document(std::string fileName)
      {
	namespace io = boost::iostreams;

	if (!boost::filesystem::exists(fileName))
	  M_throw() << "Could not open XML file";
	{ //This scopes out the file objects
	  
	  //We use the boost iostreams library to load the file into a
	  //string which may be compressed.
	  
	  //We make our filtering iostream
	  io::filtering_istream inputFile;
	  
	  //Now check if we should add a decompressor filter
	  if (std::string(fileName.end()-8, fileName.end()) == ".xml.bz2")
	      inputFile.push(io::bzip2_decompressor());
	  else if (!(std::string(fileName.end()-4, fileName.end()) == ".xml"))
	    M_throw() << "Unrecognized extension for xml file";

	  //Finally, add the file as a source
	  inputFile.push(io::file_source(fileName));
	  
	  //Force the copy to occur
	  io::copy(inputFile, io::back_inserter(_data));
	}

	_doc.parse<0>(&_data[0]);
      }
      
      //! Return the first Node with a certain name in the
      //! Document. Test if the Node is Node::valid before using it.
      //!\param name Name of the node to return.
      template<class T>
      inline Node getNode(T name) { return Node(_doc.first_node(name), &_doc); }

    protected:
      std::string _data;
      rapidxml::xml_document<> _doc;
    };

    template<>
    inline Node Document::getNode<const std::string&>(const std::string& name)
    { return Node(_doc.first_node(name.c_str()), &_doc); }

  }
}
