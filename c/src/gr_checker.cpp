/*
 * gr_checker.cpp
 *
 *  Created on: 04.01.2015
 *      Author: schlund
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "datastructs/matrix.h"
#include "datastructs/equations.h"

#include "polynomials/commutative_polynomial.h"
#include "polynomials/non_commutative_polynomial.h"

#include "semirings/pseudo_linear_set.h"
#include "semirings/semilinear_set.h"

#include "semirings/semilinSetNdd.h"

#include "parser.h"


#include "solvers/newton_generic.h"
#include "solvers/solver_utils.h"

#include "utils/string_util.h"

/*
 * Tests whether two grammars generate the same language modulo commutativity.
 * We use semilinear sets in constant-period representation to represent Parikh images
 * and check their equivalence via NDDs.
 */
int main(int argc, char* argv[]) {
  namespace po = boost::program_options;

  po::options_description generic("Generic options");
  generic.add_options()
    ( "help,h", "print this help message" )
    ( "input", po::value<std::vector<std::string> >(), "input grammars (at least two): g1 g2 [g3] [...]" )
    ;

  po::positional_options_description pos;
  pos.add("input", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(generic).positional(pos).run(), vm);
  po::notify(vm);

  SemilinSetNdd::genepi_init();

  if(vm.count("help")) {
    std::cout << generic << std::endl;
    return EXIT_SUCCESS;
  }

  std::vector<std::string> input_files = vm["input"].as< std::vector<std::string> >();
  int num_grammars = input_files.size();
  std::vector<std::string> inputs;

  if (vm.count("input") && num_grammars > 1) {
    // we are reading the input from the given files
    std::ifstream file;
    
    for(auto& filename : input_files) {
      file.open(filename, std::ifstream::in);
      if (file.fail()) {
        std::cerr << "Could not open input file: " << filename << std::endl;
      }
      std::string line;
      std::vector<std::string> input;
      while (std::getline(file, line)) {
        input.push_back(line);
      }
      // join the input into one string
      inputs.push_back( std::accumulate(input.begin(), input.end(), std::string("")) );
      file.close();
    }
  } else {
    std::cout << "Please provide at least two input files!" << std::endl;
    return EXIT_FAILURE;
  }
  
  Parser p;

  // Use SL-sets with simplification
  auto equations_fst = MakeCommEquationsAndMap(p.free_parser(inputs[0]), [](const FreeSemiring &c) -> SemilinearSetL {
    auto srconv = SRConverter<SemilinearSetL>();
    return c.Eval(srconv);
  });

  ValuationMap<SemilinearSetL> sol_fst = apply_solver<NewtonCL, CommutativePolynomial>(equations_fst, true, false, 0, false);

  bool all_equal = true;
  for(int i=1; i<num_grammars; i++) {
    auto equations = MakeCommEquationsAndMap(p.free_parser(inputs[i]), [](const FreeSemiring &c) -> SemilinearSetL {
      auto srconv = SRConverter<SemilinearSetL>();
      return c.Eval(srconv);
    });

    ValuationMap<SemilinearSetL> sol = apply_solver<NewtonCL, CommutativePolynomial>(equations, true, false, 0, false);
    if(sol != sol_fst) {
      std::cout << "Difference found" << std::endl << "0:" << result_string(sol_fst)
                                   << std::endl << i << ":" << result_string(sol) << std::endl;
      all_equal = false;
      break;
    }
  }

  if(all_equal) {
    std::cout << "All grammars equivalent modulo commutativity" << std::endl;
  }

  //TODO: cmdline-switch for overapproximation


  SemilinSetNdd::genepi_dealloc();

  return EXIT_SUCCESS;
}




