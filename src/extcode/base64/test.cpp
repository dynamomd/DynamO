#include <iostream>
#include <string>
#include "CyoDecode.h"
#include "CyoEncode.h"


int main()
{
  std::string sinput("HelloThere1");
  
  size_t enclength=CyoEncode::Base64EncodeGetLength(sinput.size());
  char* output 
    = new char[enclength+1];

  size_t length = CyoEncode::Base64Encode
    (output, sinput.c_str(), sinput.size());

  //Encode doesn't set an end char, but why would it?
  output[enclength] = '\0';

  std::cout << "\nResult is := " << output << "\nLast character is ";
  
  std::cout << output[CyoEncode::Base64EncodeGetLength(sinput.size())];

  delete[] output;

//  sinput = output;
//
//  char* output2
//    = new char[CyoDecode::Base64DecodeGetLength(sinput.size())];
//
//  length = CyoDecode::Base64Decode
//    (output2, sinput.c_str(), sinput.size());
//
//  std::cout << "\nResult is := " << output2 << "\n"; 
}
