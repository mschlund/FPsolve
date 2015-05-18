#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>

#include <boost/program_options.hpp>

#include "datastructs/matrix.h"
#include "datastructs/equations.h"


#include "polynomials/commutative_polynomial.h"
#include "polynomials/non_commutative_polynomial.h"
#include "polynomials/lossy_non_commutative_polynomial.h"


#include "semirings/commutativeRExp.h"
#include "semirings/float-semiring.h"
#include "semirings/tropical-semiring.h"
#include "semirings/pseudo_linear_set.h"
#include "semirings/semilinear_set.h"
#include "semirings/bool-semiring.h"
#include "semirings/why-set.h"
#include "semirings/viterbi-semiring.h"
#include "semirings/maxmin-semiring.h"
#include "semirings/lossy-finite-automaton.h"


#ifdef USE_GENEPI
#include "semirings/semilinSetNdd.h"
#endif

#include "parser.h"


#include "solvers/newton_generic.h"
#include "solvers/kleene_seminaive.h"
#include "solvers/solver_utils.h"

#include "utils/string_util.h"


template <typename SR, template <typename> class Poly>
ValuationMap<SR> call_solver(const std::string solver_name,  const GenericEquations<Poly, SR> &equations,
   const bool scc, const bool iteration_flag, const std::size_t iterations, const bool graphviz_output){
  if(0 == solver_name.compare("newtonSymb")) {
    std::cout << "Solver: Newton Symbolic" << std::endl;
    return apply_solver<Newton, Poly>(equations, scc, iteration_flag, iterations, graphviz_output);
  }
  else if(0 == solver_name.compare("newtonConc")) {
    std::cout << "Solver: Newton Concrete"<< std::endl;
    return apply_solver<NewtonCL, Poly>(equations, scc, iteration_flag, iterations, graphviz_output);
  }
  else if(0 == solver_name.compare("kleene")) {
    std::cout << "Solver: Kleene solver"<< std::endl;
    return apply_solver<KleeneComm, Poly>(equations, scc, iteration_flag, iterations, graphviz_output);
  }
  else {
    // default-case
    std::cout << "Solver: Newton Concrete"<< std::endl;
    return apply_solver<NewtonCL, Poly>(equations, scc, iteration_flag, iterations, graphviz_output);
  }
}

int main(int argc, char* argv[]) {

  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  desc.add_options()
    ( "scc", "apply newton method iteratively to strongly connected components of the equation graph" )
    ( "help,h", "print this help message" )
    ( "iterations,i", po::value<int>(), "specify the number of newton iterations. Default is number of equations + 1." )
    //( "verbose", "enable verbose output" )
    //( "debug", "enable debug output" )
    ( "test", "just for testing purposes ... explicit test defined in main()" )
    ( "file,f", po::value<std::string>(), "input file" )
    ( "float", "float semiring" )
    ( "bool", "boolean semiring" )
    ( "why", "why semiring" )
    ( "maxmin", "MaxMin semiring" )
    ( "viterbi", "Viterbi semiring" )
    ( "tropical", "tropical semiring over the integers" )
    ( "rexp", "commutative regular expression semiring" )
    ( "slset", "explicit semilinear sets semiring, no simplification " )
#ifdef USE_GENEPI
    ( "slsetndd", po::value<std::string>(), "ndd semilinear sets semiring, give plugin name (lash-msdf, mona)" )
    ( "n", po::value<int>(), "number of variables, used for slsetndd" )
#endif
    ( "mlset", "abstraction over semilinear sets" )
    ( "vec-simpl", "vector simplification only (only semilinear and multilinear sets)" )
    ( "lin-simpl", "linear set simplification (only semilinear sets)" )
    ( "free", "free semiring" )
    ( "lossy", "lossy semiring" )
    ( "prefix", po::value<int>(), "prefix semiring with given length")
    ( "graphviz", "create the file graph.dot with the equation graph (NOTE: currently only with option --scc) " )
    ( "solver,s", po::value<std::string>(), "solver type (currently: \"newtonSymb\", \"newtonConc\" or \"kleene\")" )
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("test")) {
    Newton<SemilinSetExp> newton;
    std::vector<VarId> variables;
    variables.push_back(Var::GetVarId("x"));
    std::cout << "- newton (cnt-SR):" << std::endl;

    std::vector<CommutativePolynomial<SemilinSetExp> > polynomials;
    CommutativePolynomial<SemilinSetExp> f1 = CommutativePolynomial<SemilinSetExp>({
        { SemilinSetExp(Var::GetVarId("a")), CommutativeMonomial{ {Var::GetVarId("x"),Var::GetVarId("x")} } },
        { SemilinSetExp(Var::GetVarId("c")), CommutativeMonomial{} } });

    polynomials.push_back(f1);

    Matrix<SemilinSetExp> result = newton.solve_fixpoint(polynomials, variables, 2);
    std::cout << result << std::endl;

    /*              auto s1 = CommutativeRExp(Var::GetVarId("a"));
                    auto s2 = CommutativeRExp(Var::GetVarId("b"));
                    auto m1 = Matrix<CommutativeRExp>(1,1,{s1});
                    auto m2 = Matrix<CommutativeRExp>(1,1,{s2});
                    */

    // this actually led to a strange bug with the ublas-matrix implementation!!
    /*              auto s1 = SemilinSetExp(Var::GetVarId("a"));
                    auto s2 = SemilinSetExp(Var::GetVarId("b"));
                    auto m1 = Matrix<SemilinSetExp>(1,1,{s1});
                    auto m2 = Matrix<SemilinSetExp>(1,1,{s2});

                    auto m3 = m1*m2;
                    std::cout << m3;
                    */
    return 0;
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

#ifdef USE_GENEPI
  SemilinSetNdd::genepi_init();
#endif

  int iterations=0;
  if (vm.count("iterations"))
          iterations = vm["iterations"].as<int>();

  // check if we can do something useful
  if (!vm.count("float") &&
      !vm.count("rexp") &&
      !vm.count("slset") &&
#ifdef USE_GENEPI
      !vm.count("slsetndd") &&
#endif
      !vm.count("free") &&
      !vm.count("bool") &&
      !vm.count("why") &&
      !vm.count("tropical") &&
      !vm.count("mlset") &&
      !vm.count("prefix") &&
      !vm.count("lossy") &&
      !vm.count("maxmin") &&
      !vm.count("viterbi"))
    {
    std::cout << "Please supply a supported semiring :)" << std::endl;
    return 0;
  }

  std::vector<std::string> input;
  std::string line;
  if (vm.count("file")) {
    // we are reading the input from the given file
    std::ifstream file;
    file.open(vm["file"].as<std::string>(), std::ifstream::in);
    if (file.fail()) {
      std::cerr << "Could not open input file: " << vm["file"].as<std::string>() << std::endl;
    }
    while (std::getline(file, line)) {
      input.push_back(line);
    }
  } else {
    // we are reading from stdin
    while (std::getline(std::cin, line)) {
      input.push_back(line);
    }
  }

  // join the input into one string
  std::string input_all =
    std::accumulate(input.begin(), input.end(), std::string(""));

  Parser p;
  std::string solver_name;

  if(vm.count("solver")) {
    solver_name = vm["solver"].as<std::string>();
  }
  else {
    solver_name = "newtonConc"; //seems to be the fastest
  }


  const auto iter_flag = vm.count("iterations");
  const auto graph_flag = vm.count("graphviz");
  const auto scc_flag = vm.count("scc");

  if (vm.count("slset")) {
    auto equations_raw = p.free_parser(input_all);
    auto equations = MakeCommEquationsAndMap(equations_raw, [](const FreeSemiring &c) -> SemilinSetExp {
      auto srconv = SRConverter<SemilinSetExp>();
      return c.Eval(srconv);
    });

    if (equations.empty()) return EXIT_FAILURE;
    //PrintEquations(equations);
    if (!vm.count("vec-simpl") && !vm.count("lin-simpl")) {
      DMSG("A");
      std::cout << result_string(
          call_solver(solver_name, equations, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;

    } else if (vm.count("vec-simpl") && !vm.count("lin-simpl")) {
      DMSG("B");
      auto equations2 = MapEquations(equations, [](const SemilinearSet<> &s) {
        return SemilinearSetV{s};
      });
      std::cout << result_string(
          call_solver(solver_name, equations2, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    } else {
      DMSG("C");
      auto equations2 = MapEquations(equations, [](const SemilinearSet<> &s) {
        return SemilinearSetL{s};
      });
      std::cout << result_string(
          call_solver(solver_name, equations2, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    }
#ifdef USE_GENEPI
  } else if (vm.count("slsetndd")) {
      SemilinSetNdd::solver_init(vm["n"].as<int>());
      auto equations = p.slsetndd_parser(input_all);
      if (equations.empty()) return EXIT_FAILURE;
      std::cout << result_string(
          call_solver(solver_name, equations, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
      SemilinSetNdd::solver_dealloc();
#endif
  } else if (vm.count("mlset")) {

    auto equations_raw = p.free_parser(input_all);
    auto equations = MakeCommEquationsAndMap(equations_raw, [](const FreeSemiring &c) -> SemilinSetExp {
      auto srconv = SRConverter<SemilinSetExp>();
      return c.Eval(srconv);
    });

    if (equations.empty()) return EXIT_FAILURE;
    if (vm.count("vec-simpl")) {
      auto m_equations = SemilinearToPseudoLinearEquations<
        DummyDivider, SparseVecSimplifier>(equations);
      //PrintEquations(m_equations);
      std::cout << result_string(
          call_solver(solver_name, equations, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    } else {
      auto m_equations = SemilinearToPseudoLinearEquations<
        DummyDivider, DummyVecSimplifier>(equations);
      //PrintEquations(m_equations);
      std::cout << result_string(
          call_solver(solver_name, equations, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    }

  } else if (vm.count("rexp")) {

    // parse the input into a list of (Var → Polynomial[SR])
    auto equations = p.rexp_parser(input_all);
    if (equations.empty()) return EXIT_FAILURE;

    //PrintEquations(equations);

    // apply solver to the equations
    std::cout << result_string(
        call_solver(solver_name, equations, scc_flag, iter_flag, iterations, graph_flag)
        ) << std::endl;

  } else if (vm.count("free")) {

    // parse the input into a list of (Var → Polynomial[SR])
	auto equations = p.free_parser(input_all);
  auto equations2 = MakeCommEquations(equations);

  if (equations2.empty()) return EXIT_FAILURE;

//    std::cout << result_string(
//        call_solver(solver_name, equations2, scc_flag, iter_flag, iterations, graph_flag)
//        ) << std::endl;

  } else if (vm.count("lossy")) {

    auto equations = p.lossy_fa_parser(input_all);
    if (equations.empty()) return EXIT_FAILURE;

    VarId S_1 = equations[0].first;

    NCEquationsBase<LossyFiniteAutomaton> eq2 = NCEquationsBase<LossyFiniteAutomaton>(equations.begin(), equations.end());

    auto approximation = NonCommutativePolynomial<LossyFiniteAutomaton>::downwardClosureCourcelle(eq2, S_1);
    std::cout << "dwc:\t" << approximation.string() << std::endl;
    std::cout << "size NFA for DWC:\t" << approximation.size() << std::endl;
    std::cout << "size minimal DFA for DWC:\t" << approximation.minimize().size() << std::endl;

  } else if (vm.count("prefix")) {

    // parse the input into a list of (Var → Polynomial[SR])
    auto equations = p.prefix_parser(input_all, vm["prefix"].as<int>());
    if (equations.empty()) return EXIT_FAILURE;

    //PrintEquations(equations);

    // apply the newton method to the equations
    //auto result = apply_solver<Kleene, PrefixSemiring>(equations,
    //                                           vm.count("scc"),
    //                                           vm.count("iterations"),
    //                                           iterations,
    //                                           vm.count("graphviz"));
    //std::cout << result_string(result) << std::endl;

    //std::cout << result_string(
    //    call_solver(solver_name, equations, scc_flag, iter_flag, iterations, graph_flag)
    //    ) << std::endl;

  } else if (vm.count("float")) {
    auto equations = p.free_parser(input_all);
    //PrintEquations(equations);
    auto equations2 = MakeCommEquationsAndMap(equations, [](const FreeSemiring &c) -> FloatSemiring {
      auto srconv = SRConverter<FloatSemiring>();
      return c.Eval(srconv);
    });

    if (equations2.empty()) return EXIT_FAILURE;

    //PrintEquations(equations);
    //PrintEquations(equations2);

    if(vm.count("solver")){
      //std::cout << "Solver: " << solver_name << std::endl;
      std::cout << result_string(
          call_solver(solver_name, equations2, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    }
    else {
#ifdef USE_NUMERICNEWTON
      //default if enabled: "Fast" numeric Newton
      std::cout << "Solver: Newton Numeric"<< std::endl;
      std::cout << result_string(
          apply_solver<NewtonNumeric>(equations2, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
#else
      //use default
      std::cout << result_string(
            call_solver("", equations2, scc_flag, iter_flag, iterations, graph_flag)
            ) << std::endl;
#endif
    }


  } else if (vm.count("bool")) {
    auto equations = p.free_parser(input_all);
    auto equations2 = MakeCommEquationsAndMap(equations, [](const FreeSemiring &c) -> BoolSemiring {
      auto srconv = SRConverter<BoolSemiring>();
      return c.Eval(srconv);
    });

    if (equations2.empty()) return EXIT_FAILURE;

    PrintEquations(equations);
    PrintEquations(equations2);
      std::cout << result_string(
          call_solver(solver_name, equations2, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;

  }
  else if (vm.count("why")) {
      auto equations = p.free_parser(input_all);
      auto equations2 = MakeCommEquationsAndMap(equations, [](const FreeSemiring &c) -> WhySemiring {
        auto srconv = SRConverter<WhySemiring>();
        return c.Eval(srconv);
      });

      if (equations2.empty()) return EXIT_FAILURE;

        std::cout << result_string(
            call_solver(solver_name, equations2, scc_flag, iter_flag, iterations, graph_flag)
            ) << std::endl;

    }
    else if (vm.count("tropical")) {
      auto equations = p.free_parser(input_all);
      auto equations2 = MakeCommEquationsAndMap(equations, [](const FreeSemiring &c) -> TropicalSemiring {
        auto srconv = SRConverter<TropicalSemiring>();
        return c.Eval(srconv);
      });

      if (equations2.empty()) return EXIT_FAILURE;

        std::cout << result_string(
            call_solver(solver_name, equations2, scc_flag, iter_flag, iterations, graph_flag)
            ) << std::endl;

    }
    else if (vm.count("viterbi")) {
          auto equations = p.free_parser(input_all);
          auto equations2 = MakeCommEquationsAndMap(equations, [](const FreeSemiring &c) -> ViterbiSemiring {
            auto srconv = SRConverter<ViterbiSemiring>();
            return c.Eval(srconv);
          });

          if (equations2.empty()) return EXIT_FAILURE;

            std::cout << result_string(
                call_solver(solver_name, equations2, scc_flag, iter_flag, iterations, graph_flag)
                ) << std::endl;
    }
    else if (vm.count("maxmin")) {
          auto equations = p.free_parser(input_all);
          auto equations2 = MakeCommEquationsAndMap(equations, [](const FreeSemiring &c) -> MaxMinSemiring {
            auto srconv = SRConverter<MaxMinSemiring>();
            return c.Eval(srconv);
          });

          if (equations2.empty()) return EXIT_FAILURE;

            std::cout << result_string(
                call_solver(solver_name, equations2, scc_flag, iter_flag, iterations, graph_flag)
                ) << std::endl;

    }

#ifdef USE_GENEPI
  SemilinSetNdd::genepi_dealloc();
#endif

  return EXIT_SUCCESS;
}
