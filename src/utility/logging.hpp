// ===========================================================================
//
//                            PUBLIC DOMAIN NOTICE
//                       Agricultural Research Service
//                  United States Department of Agriculture
//
//   This software/database is a "United States Government Work" under the
//   terms of the United States Copyright Act.  It was written as part of
//   the author's official duties as a United States Government employee
//   and thus cannot be copyrighted.  This software/database is freely
//   available to the public for use. The Department of Agriculture (USDA)
//   and the U.S. Government have not placed any restriction on its use or
//   reproduction.
//
//   Although all reasonable efforts have been taken to ensure the accuracy
//   and reliability of the software and data, the USDA and the U.S.
//   Government do not and cannot warrant the performance or results that
//   may be obtained by using this software or data. The USDA and the U.S.
//   Government disclaim all warranties, express or implied, including
//   warranties of performance, merchantability or fitness for any
//   particular purpose.
//
//   Please cite the author in any work or product based on this material.
//
//   =========================================================================
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

template<typename LevelType>
void LogInit(LevelType level)
{
  namespace logging=boost::log;

  std::map<std::string,boost::log::trivial::severity_level> severities= {
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

#ifndef SMVHIDELOG
  #define SMVLOG( x ) x
#else
  #define SMVLOG( x )
#endif

}

#endif /* _LOGGING_H_ */
