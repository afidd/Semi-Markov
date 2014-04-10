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
//   =========================================================================== 



#include "weiss.hpp"

static const counter_type major_version = 0;
static const counter_type minor_version = 1;
static const string_type version =  
       boost::lexical_cast<string_type>(major_version) + 
       "." + 
       boost::lexical_cast<string_type>(minor_version);


int main(int argc, char* argv[])
{
  po::options_description general_options("General options");

  po::options_description method_options("Method-specific options");

  po::options_description combined_options("Usage: weiss");

  po::variables_map vm;

  general_options.add_options()
    (
     "help,h", 
     "Produce help message"
     );

  general_options.add_options()
    (
     "version,v", 
     "Print version string to standard output"
     );

  general_options.add_options()
    (
    "loglevel,l",
    po::value<string_type>()->default_value("error"),
    "Logging level: trace, debug, info, warning, error, or fatal"
    );

  general_options.add_options()
    (
     "output,o", 
     po::value<string_type>()->default_value("weiss.out"), 
     "Name of output file"
     );

  general_options.add_options()
    (
     "transitions,t", 
     po::value<counter_type>()->default_value(100),
     "Number of transitions"
     );

  general_options.add_options()
    (
     "replicates,r", 
     po::value<counter_type>()->default_value(100),
     "Number of replicates"
     );

  method_options.add_options()
  (
    "seed,s", 
    po::value<boost::uint32_t>()->default_value(2342387),
    "Seed for Boost random number generator (boost::uint32_t)"
  );

  method_options.add_options()
  (
   "alpha,a",
   po::value<real_type>()->default_value(2.0),
   "Shape parameter for Weibull distribution of residence times."
  );

  method_options.add_options()
  (
   "beta,b",
   po::value<real_type>()->default_value(3.0),
   "Scale parameter for Weibull distribution of residence times."
  );

  combined_options.add(general_options);
  combined_options.add(method_options);

  po::command_line_parser clp(argc, argv);
  clp.options(combined_options);
  po::store(clp.run(), vm);
  po::notify(vm);

    
  if (vm.count("help")) {
    std::cerr << combined_options << std::endl;
    return boost::exit_success;
  } 

  afidd::LogInit(vm["loglevel"].as<string_type>());

  if (vm.count("version")) {
    std::cout << argv[0] << " version " << version << std::endl;
    return boost::exit_success;
  };

  std::ofstream os;
  string_type ofname = vm["output"].as<string_type>();
  BOOST_LOG_TRIVIAL(info) << "Opening output file " << ofname;

  os.open(ofname.c_str(), std::ofstream::out);

  BOOST_LOG_TRIVIAL(info) << "Using " << vm["seed"].as<boost::uint32_t>() << " as seed";
  RandGen rng(vm["seed"].as<boost::uint32_t>());


  using BrownionState = afidd::smv::GSPNState<Mark,UserState>;
  using SemiMarkovKernel = afidd::smv::PartialCoreMatrix<BrownionGSPN, BrownionState, RandGen>;

  auto initialize_walkers=[](BrownionState& s)->void {
    afidd::smv::Add<0>(s.marking, PlaceKey{0,0,Mobility::immobile}, Walker());
  };

  auto reporter=[](BrownionState& s) -> void 
  {
    // Intentionally empty
  };

  for (counter_type repl=0; repl < vm["replicates"].as<counter_type>(); 
       repl++) {
    BOOST_LOG_TRIVIAL(info) << "Working on replicate " << repl;

    BrownionGSPN gspn;
    BrownionState state;
    SemiMarkovKernel Q(gspn, state);
    auto next = afidd::smv::PropagateCompetingProcesses(Q,
        initialize_walkers, rng);

    real_type elapsed_time = 0.0;

    for (counter_type tcount=0; tcount<vm["transitions"].as<counter_type>();
        tcount++) {
      auto tk = std::get<0>(next);
      auto residence_time = std::get<1>(next);
      elapsed_time += residence_time;
      
      BOOST_LOG_TRIVIAL(info) << "replicate: " << repl
                              << ", transition: " << tcount
                              << ", " << tk
                              << ", residence time: " << residence_time
                              << ", elapsed time: " << elapsed_time;

      os << repl << ", " << tcount << ", " 
         << residence_time << ", " << elapsed_time << ", "
         << tk.from.x << ", " << tk.from.y << ", \"" << tk.from.mobility  << "\"" 
         <<  std::endl;
    
      next = afidd::smv::PropagateCompetingProcesses(Q, reporter, rng);
    }
  }

  BOOST_LOG_TRIVIAL(info) << "Closing output file " << ofname;
  os.close();

  return boost::exit_success;

}

