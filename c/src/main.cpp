#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/program_options.hpp>

#include "datastructs/matrix.h"

#include "polynomials/polynomial.h"
#include "polynomials/non_commutative_polynomial.h"

#include "semirings/commutativeRExp.h"
#include "semirings/float-semiring.h"
#include "semirings/pseudo_linear_set.h"
#include "datastructs/finite_automaton.h"
#include "semirings/lossy-semiring.h"
#include "semirings/lossy-regular-expression.h"
#include "semirings/lossy-finite-automaton.h"

#ifdef OLD_SEMILINEAR_SET
#include "semirings/semilinSetExp.h"
#else
#include "semirings/semilinear_set.h"
#endif

#include "utils/timer.h"


#include "parser.h"

//#include "newton.h"
#include "newton_generic.h"

template <typename SR>
struct VertexProp {
  std::string name;   // used for graphviz output
  VarId var;         // var and rex combines the equations in the vertex
  Polynomial<SR> rex;
};

// group the equations to SCCs
template <typename SR>
std::vector< std::vector< std::pair< VarId, Polynomial<SR> > > >
group_by_scc(const std::vector< std::pair< VarId, Polynomial<SR> > > &equations,
             bool graphviz_output) {
  // create map of variables to [0..n]. this is used to enumerate important variables in a clean way from 0 to n during graph construction
  std::unordered_map<VarId, int> var_key;

  // build the graph
  boost::adjacency_list<boost::vecS,
                        boost::vecS,
                        boost::bidirectionalS,
                        VertexProp<SR>> graph(equations.size());
  // for (auto e_it = equations.begin(); e_it != equations.end(); ++e_it)
  for (const auto &eq : equations) {
    // if the variable is not yet in the map insert it together with the size of
    // the map. this way we get a unique identifier for the variable counting from 0 to n
    if (var_key.find(eq.first) == var_key.end()) {
      var_key.insert(std::make_pair(eq.first, var_key.size()));
    }
    int a = var_key.find(eq.first)->second; // variable key
    graph[a].var = eq.first; // store VarId to the vertex
    graph[a].name = Var::GetVar(eq.first).string(); // store the name to the vertex
    graph[a].rex = eq.second; // store the regular expression to the vertex

    auto variables = eq.second.get_variables(); // all variables of this rule;
    // for (auto v_it = v.begin(); v_it != v.end(); ++v_it)
    for (const auto &var : variables) {
      if (var_key.find(var) == var_key.end()) // variable is not yet in the map
        var_key.insert(var_key.begin(), std::pair<VarId,int>(var,var_key.size()));
      int b = var_key.find(var)->second; // variable key
      boost::add_edge(a, b, graph);
    }
  }

  if (graphviz_output) {
    // output the created graph to graph.dot
    boost::dynamic_properties dp;
    dp.property("label", boost::get(&VertexProp<SR>::name, graph)); // vertex name is the name of the equation
    dp.property("node_id", get(boost::vertex_index, graph)); // this is needed
    std::ofstream outf("graph.dot");
    boost::write_graphviz_dp(outf, graph, dp);
  }

  // calculate strongly connected components and store them in 'component'
  std::vector<int> component(boost::num_vertices(graph));
  boost::strong_components(graph,&component[0]);

  // group neccessary equations together
  int num_comp = *std::max_element(component.begin(), component.end()) + 1; // find the number of components
  std::vector< std::vector< std::pair< VarId,Polynomial<SR> > > >
    grouped_equations(num_comp);
  // grouped_equations.resize(num_comp);

  // iterate over all vertices (0 to n)
  // collect the necessary variables + equations for every component
  for (std::size_t j = 0; j != component.size(); ++j) {
    //std::cout << j << ", " << graph[j].var << " is in component " << component[j] << std::endl;
    grouped_equations[component[j]].push_back(std::pair<VarId,Polynomial<SR>>(graph[j].var, graph[j].rex));
  }

  return grouped_equations;
}

// apply Newton's method to the given input
template <template <typename> class NewtonType = Newton, typename SR>
ValuationMap<SR> apply_newton(
    const std::vector< std::pair< VarId, Polynomial<SR> > > &equations,
    bool scc, bool iteration_flag, std::size_t iterations, bool graphviz_output) {

  // TODO: sanity checks on the input!

  // generate an instance of the newton solver
  NewtonType<SR> newton;

  // if we use the scc method, group the equations
  // the outer vector contains SCCs starting with a bottom SCC at 0
  std::vector< std::vector< std::pair< VarId,Polynomial<SR> > > > equations2;
  if (scc) {
    auto tmp = group_by_scc(equations, graphviz_output);
    equations2.insert(equations2.begin(), tmp.begin(), tmp.end());
  }
  else if (!scc) {
    equations2.push_back(equations);
  }

  // this holds the solution
  ValuationMap<SR> solution;

  Timer timer;
  timer.Start();

  // the same loop is used for both the scc and the non-scc variant
  // in the non-scc variant, we just run once through the loop
  //for (auto it1 = equations2.begin(); it != equations2.end(); ++it)
  for (std::size_t j = 0; j != equations2.size(); ++j) {
    // use the solutions to get rid of variables in the remaining equations
    // does nothing in the first round
    std::vector<std::pair<VarId, Polynomial<SR>>> tmp1;
    for (auto it = equations2[j].begin(); it != equations2[j].end(); ++it)
    { // it = (VarId, Polynomial[SR])
      auto tmp2 = it->second.partial_eval(solution);
      tmp1.push_back(std::pair<VarId, Polynomial<SR>>(it->first, tmp2));
    }
    // replace old equations with simplified ones
    equations2[j] = tmp1;

    // dynamic iterations
    if (!iteration_flag) {
      // for commutative and idempotent SRs Newton has converged after n+1 iterations, so use this number as default
      iterations = equations2[j].size() + 1;
    }

    // do some real work here
    ValuationMap<SR> result = newton.solve_fixpoint(equations2[j], iterations);

    // copy the results into the solution map
    solution.insert(result.begin(), result.end());
  }

  timer.Stop();
  std::cout << "Solving time:\t" << timer.GetMicroseconds().count()
    << " us" << std::endl
    << "Solving time:\t" << timer.GetMilliseconds().count()
    << " ms" << std::endl;


  return solution;
}

template <typename SR>
std::string result_string(const ValuationMap<SR> &result) {
  std::stringstream ss;
  for (auto &x : result) {
    ss << x.first << " == " << x.second << std::endl;
  }
  return ss.str();
}

template <typename Container>
void PrintEquations(const Container &equations) {
  std::cout << "Equations:" << std::endl;
  for (auto &eq : equations) {
    std::cout << "* " << eq.first << " → " << eq.second << std::endl;
  }
}

int main(int argc, char* argv[]) {

  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  desc.add_options()
    ( "scc", "apply newton method iteratively to strongly connected components of the equation graph" )
    ( "help,h", "print this help message" )
    ( "iterations,i", po::value<int>(), "specify the number of newton iterations. default is optimal number" )
    //( "verbose", "enable verbose output" )
    //( "debug", "enable debug output" )
    ( "test", "just for testing purposes ... explicit test defined in main()" )
    ( "file,f", po::value<std::string>(), "input file" )
    ( "float", "float semiring" )
    ( "rexp", "commutative regular expression semiring" )
    ( "slset", "explicit semilinear sets semiring, no simplification " )
    ( "mlset", "abstraction over semilinear sets" )
    ( "vec-simpl", "vector simplification only (only semilinear and multilinear sets)" )
    ( "lin-simpl", "linear set simplification (only semilinear sets)" )
    ( "free", "free semiring" )
    ( "lossy", "lossy semiring" )
    ( "prefix", po::value<int>(), "prefix semiring with given length")
    ( "graphviz", "create the file graph.dot with the equation graph" )
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("test")) {
    Newton<SemilinSetExp> newton;
    std::vector<VarId> variables;
    variables.push_back(Var::GetVarId("x"));
    std::cout << "- newton (cnt-SR):" << std::endl;

    std::vector<Polynomial<SemilinSetExp> > polynomials;
    Polynomial<SemilinSetExp> f1 = Polynomial<SemilinSetExp>({
        { SemilinSetExp(Var::GetVarId("a")), Monomial{ {Var::GetVarId("x"),Var::GetVarId("x")} } },
        { SemilinSetExp(Var::GetVarId("c")), Monomial{} } });

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

  int iterations=0;
  if (vm.count("iterations"))
          iterations = vm["iterations"].as<int>();

  // check if we can do something useful
  if (!vm.count("float") &&
      !vm.count("rexp") &&
      !vm.count("slset") &&
      !vm.count("free") &&
      !vm.count("mlset") &&
      !vm.count("prefix") &&
      !vm.count("lossy")) {
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

  const auto iter_flag = vm.count("iterations");
  const auto graph_flag = vm.count("graphviz");
  const auto scc_flag = vm.count("scc");

  if (vm.count("slset")) {

    auto equations = p.slset_parser(input_all);
    if (equations.empty()) return EXIT_FAILURE;
    //PrintEquations(equations);
    if (!vm.count("vec-simpl") && !vm.count("lin-simpl")) {
      DMSG("A");
      std::cout << result_string(
          apply_newton(equations, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    } else if (vm.count("vec-simpl") && !vm.count("lin-simpl")) {
      DMSG("B");
      auto equations2 = MapEquations(equations, [](const SemilinearSet<> &s) {
        return SemilinearSetV{s};
      });
      std::cout << result_string(
          apply_newton(equations2, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    } else {
      DMSG("C");
      auto equations2 = MapEquations(equations, [](const SemilinearSet<> &s) {
        return SemilinearSetL{s};
      });
      std::cout << result_string(
          apply_newton(equations2, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    }

  } else if (vm.count("mlset")) {

    auto equations = p.slset_parser(input_all);
    if (equations.empty()) return EXIT_FAILURE;
    if (vm.count("vec-simpl")) {
      auto m_equations = SemilinearToPseudoLinearEquations<
        DummyDivider, SparseVecSimplifier>(equations);
      //PrintEquations(m_equations);
      std::cout << result_string(
          apply_newton(m_equations, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    } else {
      auto m_equations = SemilinearToPseudoLinearEquations<
        DummyDivider, DummyVecSimplifier>(equations);
      //PrintEquations(m_equations);
      std::cout << result_string(
          apply_newton(m_equations, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;
    }


  } else if (vm.count("rexp")) {

    // parse the input into a list of (Var → Polynomial[SR])
    auto equations = p.rexp_parser(input_all);
    if (equations.empty()) return EXIT_FAILURE;

    //PrintEquations(equations);

    // apply Newton's method to the equations
    auto result = apply_newton(equations, vm.count("scc"),
                               vm.count("iterations"), iterations,
                               vm.count("graphviz"));
    std::cout << result_string(result) << std::endl;

  } else if (vm.count("free")) {

    // parse the input into a list of (Var → Polynomial[SR])
	auto equations = p.free_parser(input_all);
    if (equations.empty()) return EXIT_FAILURE;

    PrintEquations(equations);

    // apply the newton method to the equations
    //auto result = apply_newton<FreeSemiring>(equations,
    //                                         vm.count("scc"),
    //                                         vm.count("iterations"),
    //                                         iterations,
    //                                         vm.count("graphviz"));
    //std::cout << result_string(result) << std::endl;

  } else if (vm.count("lossy")) {

    // parse to lossy semiring polynomial; an element "a" of the semiring
    // will be parsed to "1+a" while variables do not get the "1+" bit
    auto equations = p.lossy_fa_parser(input_all);
    std::cout << "hooray for parsing" << std::endl;
    if (equations.empty()) return EXIT_FAILURE;

    ValuationMap<LossyFiniteAutomaton> valuation = LossyFiniteAutomaton::solvePolynomialSystem(equations);

    PrintEquations(equations);
    std::cout << "Fixpoints:" << std::endl;
    std::cout << result_string(valuation) << std::endl;

//    for (auto it = valuation.begin(); it != valuation.end(); ++it ) {
//        std::cout << "fixpoint of " << it->first << ": " << it->second;
//        std::cout << std::endl;
//    }

  } else if (vm.count("prefix")) {

    // parse the input into a list of (Var → Polynomial[SR])
    auto equations = p.prefix_parser(input_all, vm["prefix"].as<int>());
    if (equations.empty()) return EXIT_FAILURE;

    //PrintEquations(equations);

    // apply the newton method to the equations
    //auto result = apply_newton<PrefixSemiring>(equations,
    //                                           vm.count("scc"),
    //                                           vm.count("iterations"),
    //                                           iterations,
    //                                           vm.count("graphviz"));
    //std::cout << result_string(result) << std::endl;

  } else if (vm.count("float")) {

    auto equations = p.float_parser(input_all);
    if (equations.empty()) return EXIT_FAILURE;

    //PrintEquations(equations);
      std::cout << result_string(
          apply_newton<NewtonCL>(equations, scc_flag, iter_flag, iterations, graph_flag)
          ) << std::endl;

  }



  return EXIT_SUCCESS;
}
