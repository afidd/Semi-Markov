#ifndef _LOGGING_H_
#define _LOGGING_H_ 1

#include <iostream>
#include <map>
#include "stochnet.hpp"
#include "boost/log/core.hpp"
#include "boost/log/expressions.hpp"
#include "boost/algorithm/string/join.hpp"


namespace afidd
{

void LogInit(std::string level)
{
  namespace logging=boost::log;

  std::map<std::string,boost::log::trivial::severity_level> severities=
    {
      { "trace",   logging::trivial::trace },
      { "debug",   logging::trivial::debug },
      { "info",    logging::trivial::info },
      { "warning", logging::trivial::warning },
      { "error",   logging::trivial::error },
      { "fatal",   logging::trivial::fatal }
    };

  logging::trivial::severity_level assign_level=logging::trivial::info;
  auto logiter=severities.find(level);
  if (logiter!=severities.end()) {
    assign_level=logiter->second;
  } else {
    std::cout << "Could not set the logging level from " << level
        << ". Choices are: ";
    std::vector<std::string> names;
    for (auto& kv : severities) {
      names.push_back(kv.first);
    }
    std::cout << boost::algorithm::join(names, ", ")
        << ". Choosing info." << std::endl;
  }

  logging::core::get()->set_filter(logging::trivial::severity >= assign_level);
}


}

#endif /* _LOGGING_H_ */
