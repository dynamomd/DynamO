/*  dynamo:- Event driven molecular dynamics simulator
    http://www.dynamomd.org
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

#include <boost/lexical_cast.hpp>
#include <magnet/exception.hpp>
#include <rapidXML/rapidxml.hpp>
#ifdef DYNAMO_bzip2_support
#include <bzlib.h>
#endif
#include <fstream>
#include <iostream>
#include <vector>

namespace magnet {
//! Namespace enclosing the XML tools included in magnet.
namespace xml {
namespace detail {
/*! \brief Generates a string representation of the passed XML
    rapidxml node.
 */
inline std::string getPath(rapidxml::xml_node<> *_node) {
  std::vector<std::pair<std::string, int>> pathTree;
  while (_node != 0) {
    // Determine what is the index of the current node.
    int index(0);
    if (_node->parent()) {
      rapidxml::xml_node<> *sibling =
          _node->previous_sibling(_node->name(), _node->name_size());
      while (sibling) {
        sibling = sibling->previous_sibling();
        ++index;
      }
      // Set the index to -1 if this is the only node with its name (it has no
      // siblings)
      if ((index == 0) &&
          !(_node->next_sibling(_node->name(), _node->name_size())))
        index = -1;
    }

    pathTree.push_back(std::make_pair(
        std::string(_node->name(), _node->name() + _node->name_size()), index));
    _node = _node->parent();
  }

  std::ostringstream os;

  for (std::vector<std::pair<std::string, int>>::const_reverse_iterator iPtr =
           pathTree.rbegin();
       iPtr != pathTree.rend(); ++iPtr) {
    os << "/" << iPtr->first;
    // Only output the index of the node if it has no siblings
    if (iPtr->second >= 0)
      os << "[" << iPtr->second << "]";
  }
  return os.str();
}
} // namespace detail

/*! \brief Represents an Attribute of an XML Document.
 */
class Attribute {
public:
  //! \brief Converts the attributes value to a type.
  template <class T> inline T as() const {
    try {
      return boost::lexical_cast<T>(getValue());
    } catch (boost::bad_lexical_cast &) {
      M_throw() << "The value \"" << getValue()
                << "\" will not cast to the correct type. Please check the "
                   "attribute at the following XMLPath: "
                << getPath();
    }
  }

  /*! \brief Returns the value of the attribute.
   */
  std::string getValue() const {
    if (!valid())
      M_throw()
          << (std::string(
                  "XML error: Missing attribute being converted\nXML Path: ") +
              detail::getPath(_parent) + "/INVALID");

    return std::string(_attr->value(), _attr->value() + _attr->value_size());
  }

  /*! \brief A convenient method to access the attributes value.
   */
  operator std::string() const { return getValue(); }

  /*! \brief Test if this Attribute is valid and has a value.
   */
  inline bool valid() const { return _attr != NULL; }

  /*! \brief Returns a string representation of this attributes
    XML location.
   */
  inline std::string getPath() const {
    if (!valid())
      M_throw() << (std::string("XML error: Cannot get path of invalid "
                                "attribute\nXML Path: ") +
                    detail::getPath(_parent) + "/INVALID");
    return detail::getPath(_parent) + std::string("/@") +
           std::string(_attr->name(), _attr->name() + _attr->name_size());
  }

  inline std::string getName() const {
    return std::string(_attr->name(), _attr->name() + _attr->name_size());
  }

private:
  friend class Node;

  /*! \brief Hidden default constructor to stop its use.
   */
  Attribute();

  /*! \brief True Attribute constructor, used by the Node class.
   */
  inline Attribute(rapidxml::xml_attribute<> *attr,
                   rapidxml::xml_node<> *parent)
      : _attr(attr), _parent(parent) {}

  rapidxml::xml_attribute<> *_attr;
  rapidxml::xml_node<> *_parent;
};

/*! \brief Represents a Node of an XML Document.
 */
class Node {
public:
  /*! \brief Fetches the first Attribute with a given name.

    Test if the Attribute is Attribute::valid() before using it.
    \param name The name of the attribute to return.
   */
  template <class T> inline Attribute getAttribute(T name) const {
    if (!valid())
      M_throw() << (std::string("XML error: Invalid node's attribute being "
                                "accessed\nXML Path:") +
                    detail::getPath(_parent) + "/INVALID");

    Attribute attr(_node->first_attribute(name), _node);

    if (!attr.valid())
      M_throw() << "XML error: Attribute \"" << name << "\" does not exist."
                << (std::string("\nXML Path: ") + detail::getPath(_node) +
                    "/@" + std::string(name));

    return attr;
  }

  /*! \brief Fetches the first Node with a given name.

    \param name The name of the Node to return.
   */
  template <class T> inline Node getNode(T name) const {
    Node child(findNode(name));

    if (!child.valid())
      M_throw() << "XML error: Node \"" << name << "\" does not exist."
                << (std::string("\nXML Path: ") + detail::getPath(_node) + "/" +
                    std::string(name));
    return child;
  }

  /*! \brief Fetches the first Node with a given name without
    testing if its a valid node.

    \param name The name of the Node to return.
   */
  template <class T> inline Node findNode(T name) const {
    if (!valid())
      M_throw() << (std::string("XML error: Cannot fetch sub node of invalid "
                                "node\nXML Path: ") +
                    detail::getPath(_parent) + "/INVALID");
    return Node(_node->first_node(name), _node);
  }

  /*! \brief Tests if the Node has a named child Node.

    \param name The name of the child Node to test for.
   */
  template <class T> inline bool hasNode(T name) const {
    if (!valid())
      M_throw() << (std::string("XML error: Cannot fetch sub node of invalid "
                                "node\nXML Path: ") +
                    detail::getPath(_parent) + "/INVALID");

    return _node->first_node(name);
  }

  /*! \brief Tests if the Node has a named child Node.

    \param name The name of the child Node to test for.
   */
  template <class T> inline bool hasAttribute(T name) const {
    if (!valid())
      M_throw() << (std::string("XML error: Cannot fetch attribute of invalid "
                                "node\nXML Path: ") +
                    detail::getPath(_parent) + "/@" + std::string(name));
    return _node->first_attribute(name);
  }

  //! \brief Returns the value of the Node.
  std::string getValue() const {
    if (!valid())
      M_throw() << (std::string("XML error: Cannot get the value of an invalid "
                                "node\nXML Path: ") +
                    detail::getPath(_parent) + "/INVALID");

    return std::string(_node->value(), _node->value() + _node->value_size());
  }

  //! \brief Test if the Node is valid.
  inline bool valid() const { return _node != NULL; }

  /*! \brief Replace this Node with the next Node in the parent
    Node with the same name.
   */
  inline void operator++() {
    if (!valid())
      M_throw()
          << (std::string(
                  "XML error: Cannot increment invalid node\nXML Path: ") +
              detail::getPath(_parent) + "/INVALID");
    _node = _node->next_sibling(_node->name(), _node->name_size());
  }

  /*! \brief A convenient method to access the value of the node.
   */
  operator std::string() const { return getValue(); }

  /*! \brief Replace this Node with the previous Node in the
    parent Node with the same name.
   */
  inline void operator--() {
    if (!valid())
      M_throw()
          << (std::string(
                  "XML error: Cannot decrement invalid node\nXML Path: ") +
              detail::getPath(_parent) + "/INVALID");
    _node = _node->previous_sibling(_node->name(), _node->name_size());
  }

  //! \brief Fetch the parent Node of this Node.
  inline Node getParent() const {
    if (_parent)
      return Node(_parent, _parent->parent());
    else
      M_throw() << "No parent node for node " << detail::getPath(_node);
  }

  //! \brief Constructor to build a Node from rapidxml::xml_node data.
  inline Node(rapidxml::xml_node<> *node, rapidxml::xml_node<> *parent)
      : _node(node), _parent(parent) {}

  /*! \brief Returns a string representation of the Node's
    location in the XML document.
   */
  inline std::string getPath() const {
    if (!valid())
      M_throw()
          << (std::string(
                  "XML error: Cannot get path of invalid node\nXML Path: ") +
              detail::getPath(_parent) + "/INVALID");

    return detail::getPath(_node);
  }

  inline std::string getName() const {
    return std::string(_node->name(), _node->name() + _node->name_size());
  }

private:
  rapidxml::xml_node<> *_node;
  rapidxml::xml_node<> *_parent;
};

template <>
inline Attribute Node::getAttribute<std::string>(std::string name) const {
  if (!valid())
    M_throw() << (std::string("XML error: Missing node's attribute being "
                              "accessed\nXML Path: ") +
                  detail::getPath(_parent) + "/INVALID");
  return Attribute(_node->first_attribute(name.c_str()), _node);
}

template <> inline Node Node::getNode<std::string>(std::string name) const {
  if (!valid())
    M_throw() << (std::string("XML error: Cannot fetch sub node of invalid "
                              "node\nXML Path: ") +
                  detail::getPath(_parent) + "/INVALID");
  return Node(_node->first_node(name.c_str()), _node);
}

/*! \brief A class which represents a whole XML Document, including
  storage.

  This class thinly wraps the rapidXML parser, providing a simple
  type-casting interface and automatic error handling, similar to
  Boost's property_tree. Unlike property_tree this class is
  space-efficient and does not copy the XML data.  This class
  must outlive any Node or Attribute generated from it.

  To parse a document, you must first load the document text into
  this class's _data parameter using \ref
  getStoredXMLData(). Then you call \ref parseData().
*/
class Document {
public:
  /*! \brief Decompress (if needed) and parse an XML file. */
  Document(std::string filename) {
    _data.clear();

    if (std::string(filename.end() - 4, filename.end()) == ".bz2") {
#ifdef DYNAMO_bzip2_support
      FILE *f = fopen(filename.c_str(), "r");
      if (!f) {
        M_throw() << "Failed to open " << filename << " for reading.";
      }
      int bzerror;
      BZFILE *b = BZ2_bzReadOpen(&bzerror, f, 0, 0, NULL, 0);
      if (bzerror != BZ_OK) {
        BZ2_bzReadClose(&bzerror, b);
        fclose(f);
        M_throw() << "Failed beginning decompression of " << filename
                  << " for reading.";
      }
      char buf[1024 * 10];
      bzerror = BZ_OK;
      while (bzerror == BZ_OK) {
        size_t nBuf = BZ2_bzRead(&bzerror, b, buf, sizeof(buf));
        if ((bzerror == BZ_OK) || (bzerror == BZ_STREAM_END))
          _data.append(buf, nBuf);
      }

      if (bzerror != BZ_STREAM_END) {
        BZ2_bzReadClose(&bzerror, b);
        fclose(f);
        M_throw() << "Failed while decompressing " << filename
                  << " for reading. (bzerror=" << bzerror << ")";
      } else {
        BZ2_bzReadClose(&bzerror, b);
        fclose(f);
      }
#else
      M_throw() << "bz2 compressed file support was not built in! (only "
                   "available on linux)";
#endif
    } else {
      std::ifstream t(filename);
      if (!t.is_open())
        M_throw() << "Failed to open " << filename << " for reading.";
      t.seekg(0, std::ios::end);
      _data.reserve(t.tellg());
      t.seekg(0, std::ios::beg);
      _data.assign((std::istreambuf_iterator<char>(t)),
                   std::istreambuf_iterator<char>());
    }
    parseData();
  }

  /*! \brief Return the first root node with a certain name in the
    Document.

    \param name Name of the node to return.
  */
  inline Node getNode(const std::string &name) {
    Node node(_doc.first_node(name.c_str()), &_doc);
    if (!node.valid())
      M_throw() << "XML error: Root node \"" << name << "\" does not exist.";
    return node;
  }

protected:
  /*! \brief Parse the stored XML data.
   */
  inline void parseData() {
    try {
      // Parse in non-destructive mode to allow verbose error reporting
      _doc.parse<rapidxml::parse_non_destructive |
                 rapidxml::parse_validate_closing_tags |
                 rapidxml::parse_trim_whitespace>(&_data[0]);
    } catch (rapidxml::parse_error &err) {
      const char *error_loc_ptr = err.where<char>();

      // Find the line of the error
      size_t line_num = 1;
      for (const char *ptr = &_data[0]; ptr < error_loc_ptr; ++ptr)
        if (*ptr == '\n')
          ++line_num;

      // Determine the start of the error line
      const char *error_line_start = error_loc_ptr;
      while ((*error_line_start != '\n') && (error_line_start != &_data[0]))
        --error_line_start;
      ++error_line_start;

      // Determine the end of the error line
      const char *error_line_end = error_loc_ptr;
      while ((*error_line_end != '\n') && (*error_line_end != '\0'))
        ++error_line_end;

      M_throw() << "Parser error at line " << line_num << ": " << err.what()
                << "\n"
                << std::string(error_line_start, error_line_end) << "\n"
                << std::string(error_loc_ptr - error_line_start, ' ') << "^";
    }
  }

  std::string _data;
  rapidxml::xml_document<> _doc;
};

template <> inline bool Node::hasAttribute(std::string name) const {
  return hasAttribute(name.c_str());
}

} // namespace xml
} // namespace magnet
