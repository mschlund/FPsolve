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
#include "utils/timer.h"



// check whether a set of grammars generates the same language up to commuativity
// (and modulo additional overapproximations given by the semiring)
template <typename SR>
void check_all_equal_commutative(const std::string& startsymbol, const std::vector<std::string>& inputs) {

  Parser p;
  int num_grammars = inputs.size();

  auto nc_equations = p.free_parser(inputs[0]);

  std::cout << "Eq (non-comm) : " << std::endl;
  PrintEquations(nc_equations);

  // Use appropriate semiring (has to be commutative!)
  auto equations_fst = MakeCommEquationsAndMap(nc_equations, [](const FreeSemiring &c) -> SR {
    auto srconv = SRConverter<SR>();
    return c.Eval(srconv);
  });

  std::cout << "Eq (comm) : "  << std::endl;
  PrintEquations(equations_fst);

  Timer timer;
  timer.Start();

  ValuationMap<SR> sol_fst = apply_solver<NewtonCL, CommutativePolynomial>(equations_fst, true, false, 0, false);

  bool all_equal = true;
  for(int i=1; i<num_grammars; i++) {
    auto equations = MakeCommEquationsAndMap(p.free_parser(inputs[i]), [](const FreeSemiring &c) -> SR {
      auto srconv = SRConverter<SR>();
      return c.Eval(srconv);
    });

    ValuationMap<SR> sol = apply_solver<NewtonCL, CommutativePolynomial>(equations, true, false, 0, false);


    if(startsymbol.compare("") == 0) {
      if(sol[equations[0].first] != sol_fst[equations_fst[0].first]) {
        std::cout << "[DIFF] Difference found for startsymbols (" << equations_fst[0].first << "," << equations[0].first << ")" << std::endl;
        std::cout << "0:" << result_string(sol_fst) << std::endl << i << ":" << result_string(sol) << std::endl;
        all_equal = false;
        break;
      }
    }
    else {

      if(sol.find(Var::GetVarId(startsymbol)) == sol.end() || sol_fst.find(Var::GetVarId(startsymbol)) == sol_fst.end()) {
        std::cout << "[ERROR] startsymbol (" << startsymbol << ") does not occur!"<< std::endl;
        return;
      }
      else if(sol[Var::GetVarId(startsymbol)] != sol_fst[Var::GetVarId(startsymbol)]) {
        std::cout << "[DIFF] Difference found for startsymbol (" << startsymbol << ")" << std::endl << "0:" << result_string(sol_fst)
                                               << std::endl << i << ":" << result_string(sol) << std::endl;
        all_equal = false;
        break;
      }
    }

  }

  if(all_equal) {
    std::cout << "[EQUIV] All grammars equivalent modulo commutativity" << std::endl;
  }

  timer.Stop();
  std::cout
  << "Total checking time:\t" << timer.GetMilliseconds().count()
  << " ms" << " ("
  << timer.GetMicroseconds().count()
  << "us)" << std::endl;

}


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
        ( "startsymbol,s", po::value<std::string>(), "start symbol of the grammars")
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

  std::string startsymbol = "";

  if(vm.count("startsymbol")) {
    startsymbol = vm["startsymbol"].as<std::string>();
    std::cout << "Comparing startsymbols (" << startsymbol << ")" << std::endl;
  }
  else {
    std::cout << "Comparing *all* nonterminals." << std::endl;
  }

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


   std::cout << "Plain Semilinear Sets" << std::endl;
  // no overapproximation -- just plain semilinear sets with simplification
  // (the approximations are not sound for inequivalence-testing!)
  check_all_equal_commutative<SemilinearSetL>(startsymbol, inputs);


  /* TODO: incorporate lossy-check
   *
   *      auto equations_2 = p.lossy_fa_parser(input_2_all);
          VarId S_2 = equations_2[0].first;

          int refinementDepth = 0;
          if(vm.count("refine")) {
              refinementDepth = vm["refine"].as<int>();
          }

          auto witness = LossyFiniteAutomaton::refineCourcelle(equations, S_1, equations_2, S_2, refinementDepth);

          if(witness != LossyFiniteAutomaton::null()) {
              std::cout << "Found witness: " << witness.string() << std::endl;
          }
   *
   *
   *
   */

  SemilinSetNdd::genepi_dealloc();

  return EXIT_SUCCESS;
}




