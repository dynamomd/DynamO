#include <iostream>
#include "hamigaki/iostreams/filter/base64.hpp"
#include <string>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include "../../xmlwriter.hpp"
#include <boost/iostreams/device/stream_sink.hpp>
#include <boost/iostreams/device/stream_source.hpp>
#include "../../../base/is_exception.hpp"

int main()
{
  namespace io = boost::iostreams;
  
  std::string filename("filename.out.bz2");

  {
    io::filtering_ostream coutputFile;
    coutputFile.push(io::bzip2_compressor());
    coutputFile.push(io::file_sink(filename));

    coutputFile << "Heres a bunch of text stuff that has no meaning or length"
		<< "\n<EOXML />\n";
    
    float test(3.14159265);

    io::filtering_ostream base64Convertor;
    base64Convertor.push(hamigaki::iostreams::base64_encoder());
    base64Convertor.push(io::stream_sink<io::filtering_ostream>(coutputFile));
    
    base64Convertor << "A test sentance for base 64 conversion";
    //base64Convertor.write(reinterpret_cast<const char*>(&test), sizeof(test));
  }

  {
    io::filtering_istream cinputFile;
    cinputFile.push(io::bzip2_decompressor());
    cinputFile.push(io::file_source(filename));

    bool found(false);
    while (!cinputFile.eof())
      {
	std::string line;
	std::getline(cinputFile, line);
	if (line == "<EOXML />") 
	  {
	    found = true;
	    break;
	  }
      }

    if (!found) D_throw() << "Could not find the End of XML marker";

    //Switch to base64 

    io::filtering_istream base64Convertor;
    base64Convertor.push(hamigaki::iostreams::base64_decoder());
    base64Convertor.push(io::stream_source<io::filtering_istream>(cinputFile));

    
    while (!base64Convertor.eof())
      {
	std::string line;
	std::getline(base64Convertor, line);
	std::cout << line << "\n";
      }
    /*float test;
    base64Convertor.read(reinterpret_cast<char *>(&test), sizeof(float));    
    std::cout << "The value is " << test << "\n";*/
  }
}
