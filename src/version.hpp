#pragma once
#include <string>

namespace DYNAMO { 
  struct buildInfo {
    static const std::string gitCheckoutHash;
    static const std::string buildDate;
    static const std::string toolChain;
    static const std::string buildHost;
  };
}
