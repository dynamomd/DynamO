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
      //! \brief Generates a string representation of the passed XML rapidxml node.
      inline std::string getPath(rapidxml::xml_node<> *_node)
      {
	std::vector<std::pair<std::string, size_t> > pathTree;
	while (_node != 0)
	  {
	    size_t index(0);
	    if (_node->parent())
	      { rapidxml::xml_node<>* sibling = _node->previous_sibling(_node->name());
		while (sibling) { sibling = sibling->previous_sibling(); ++index; }
	      }
	    
	    pathTree.push_back(std::pair<std::string, size_t>(_node->name(), index));
	    _node = _node->parent();
	  }
	
	std::ostringstream os;
	
#ifdef __PATHCC__ 
	os << "NO PATH FOR PATHSCALE COMPILER";
#else
	for (std::vector<std::pair<std::string, size_t> >::const_reverse_iterator 
	       iPtr = pathTree.rbegin(); iPtr != pathTree.rend(); ++iPtr)
	  os << "/" << iPtr->first << "[" << iPtr->second << "]";
#endif	
	return os.str();
      }
    }

    //! \brief Represents an Attribute of an XML Document.
    class Attribute {
    public:
      //! \brief Converts the attributes value to a type.
      template<class T> inline T as() const 
      { 
	if (!valid()) 
	  M_throw() << (std::string("XML error: Missing attribute being converted\nXML Path: ")
			+ detail::getPath(_parent) + "/INVALID");
	try {
	  return boost::lexical_cast<T>(_attr->value()); 
	} catch (boost::bad_lexical_cast&)
	  {
	    M_throw() << "Failed to cast XML attribute to type, XMLPath: " << getPath();
	  }
      }

      //! \brief Conversion operator returning the value of the attribute.
      inline operator const char*() const 
      { 
	if (!valid()) 
	  M_throw() << (std::string("XML error: Missing attribute being converted\nXML Path: ")
			+ detail::getPath(_parent) + "/INVALID");
	return _attr->value(); 
      }

      //! \brief Returns the value of the attribute.
      const char* getValue() const { return _attr->value(); }

      //! \brief Test if this Attribute is valid and has a value.
      inline bool valid() const { return _attr != NULL; }

      /*! \brief Returns a string representation of this attributes
       * XML location.
       */
      inline std::string getPath() const
      {	
	if (!valid()) 
	  M_throw() << (std::string("XML error: Cannot get path of invalid attribute\nXML Path: ")
			+ detail::getPath(_parent) + "/INVALID");
	return detail::getPath(_parent) + std::string("/@") + std::string(_attr->name()); 
      }

    private:
      friend class Node;

      //! \brief Hidden default constructor to stop its use.
      Attribute();

      //! \brief True Attribute constructor, used by the Node class.
      inline Attribute(rapidxml::xml_attribute<>* attr,
		       rapidxml::xml_node<>* parent):_attr(attr), _parent(parent) {}
      
      rapidxml::xml_attribute<> *_attr;
      rapidxml::xml_node<> *_parent;
    };

    //! \brief Represents a Node of an XML Document.
    class Node {
    public:
      /*! \brief Fetches the first Attribute with a given name. 
       *
       * Test if the Attribute is Attribute::valid() before using it.
       * \param name The name of the attribute to return.
       */
      template<class T> inline Attribute getAttribute(T name) const 
      { 
	if (!valid()) 
	  M_throw() << (std::string("XML error: Invalid node's attribute being accessed\nXML Path:")
			+ detail::getPath(_parent) + "/INVALID");

	Attribute attr(_node->first_attribute(name), _node);

	if (!attr.valid())
	  M_throw() << "XML error: Attribute does not exist."
		    << name << (std::string("\nXML Path: ") + detail::getPath(_node) + "/@" + std::string(name));

	return attr; 
      }

      /*! \brief Fetches the first Node with a given name.
       *
       * \param name The name of the Node to return.
       */
      template<class T> inline Node getNode(T name) const 
      { 
	Node child(fastGetNode(name));
	
	if (!child.valid())
	  M_throw() << "XML error: Node does not exist."
		    << (std::string("\nXML Path: ") + detail::getPath(_node) + "/" + std::string(name));
	return child;
      }

      /*! \brief Fetches the first Node with a given name without
       * testing if its a valid node.
       *
       * \param name The name of the Node to return.
       */
      template<class T> inline Node fastGetNode(T name) const 
      { 
	if (!valid())
	  M_throw() << (std::string("XML error: Cannot fetch sub node of invalid node\nXML Path: ")
			+ detail::getPath(_parent) + "/INVALID");
	return Node(_node->first_node(name), _node);
      }

      /*! \brief Tests if the Node has a named child Node.
       *
       * \param name The name of the child Node to test for.
       */
      template<class T> inline bool hasNode(T name) const 
      { 
	if (!valid()) 
	  M_throw() << (std::string("XML error: Cannot fetch sub node of invalid node\nXML Path: ")
			+ detail::getPath(_parent) + "/INVALID");
	
	return _node->first_node(name);
      }

      /*! \brief Tests if the Node has a named child Node.
       *
       * \param name The name of the child Node to test for.
       */
      template<class T> inline bool hasAttribute(T name) const 
      { 
	if (!valid()) 
	  M_throw() << (std::string("XML error: Cannot fetch attribute of invalid node\nXML Path: ")
			+ detail::getPath(_parent) + "/@" + std::string(name));	
	return _node->first_attribute(name);
      }

      //! \brief Returns the value of the Node.
      inline operator const char*() const { return getValue(); }

      //! \brief Returns the value of the Node.
      const char* getValue() const
      { 
	if (!valid()) 
	  M_throw() << (std::string("XML error: Cannot get the value of an invalid node\nXML Path: ")
			+ detail::getPath(_parent) + "/INVALID");
	return _node->value(); 
      }

      //! \brief Test if the Node is valid.
      inline bool valid() const { return _node != NULL; }

      /*! \brief Replace this Node with the next Node in the parent
       * Node with the same name.
       */
      inline void operator++()
      { 
	if (!valid()) 
	  M_throw() << (std::string("XML error: Cannot increment invalid node\nXML Path: ")
			+ detail::getPath(_parent) + "/INVALID");
	_node = _node->next_sibling(_node->name()); 
      }

      /*! \brief Replace this Node with the previous Node in the
       * parent Node with the same name.
       */
      inline void operator--() 
      { 
	if (!valid()) 
	  M_throw() << (std::string("XML error: Cannot decrement invalid node\nXML Path: ")
			+ detail::getPath(_parent) + "/INVALID");
	_node = _node->previous_sibling(_node->name()); 
      }

      //! \brief Fetch the parent Node of this Node.
      inline Node getParent() const
      {
	if (_parent)
	  return Node(_parent, _parent->parent());
	else
	  M_throw() << "No parent node for node " << detail::getPath(_node);
      }

      //! \brief Constructor to build a Node from rapidxml::xml_node data.
      inline Node(rapidxml::xml_node<>* node,
		  rapidxml::xml_node<>* parent):_node(node), _parent(parent) {}

      /*! \brief Returns a string representation of the Node's
       * location in the XML document.
       */
      inline std::string getPath() const 
      { 
	if (!valid()) 
	  M_throw() << (std::string("XML error: Cannot get path of invalid node\nXML Path: ")
			+ detail::getPath(_parent) + "/INVALID");

	return detail::getPath(_node); 
      }
   private:
      rapidxml::xml_node<> *_node;
      rapidxml::xml_node<> *_parent;
    };

    template<>
    inline Attribute Node::getAttribute<std::string>(std::string name) const
    { 
      if (!valid()) 
	M_throw() << (std::string("XML error: Missing node's attribute being accessed\nXML Path: ")
		      + detail::getPath(_parent) + "/INVALID");
      return Attribute(_node->first_attribute(name.c_str()), _node); 
    }

    template<>
    inline Node Node::getNode<std::string>(std::string name) const
    { 
      if (!valid()) 
	M_throw() << (std::string("XML error: Cannot fetch sub node of invalid node\nXML Path: ")
		      + detail::getPath(_parent) + "/INVALID");
      return Node(_node->first_node(name.c_str()), _node); 
    }


    /*! \brief A class which represents a whole XML Document, including
     * storage. 
     *
     * This class thinly wraps the rapidXML parser, providing a simple
     * type-casting interface and automatic error handling, similar to
     * Boost's property_tree. Unlike property_tree this class is
     * space-efficient and does not copy the XML data.  This class
     * must outlive any Node or Attribute generated from it.
     */
    class Document {
    public:
      /*! \brief Construct an XML Document from a file.
       *
       * \param fileName Path to the file to load.
       */
      inline Document(std::string fileName)
      {
	namespace io = boost::iostreams;

	if (!boost::filesystem::exists(fileName))
	  M_throw() << "Could not find the XML file named " << fileName
		    << "\nPlease check the file exists.";
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

	_doc.parse<rapidxml::parse_trim_whitespace>(&_data[0]);
      }
      
      /*! \brief Return the first Node with a certain name in the
       * Document.
       * 
       * Test if the Node is Node::valid before using it.
       * \param name Name of the node to return.
       */
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
