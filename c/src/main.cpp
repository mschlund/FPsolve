#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

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

int64_t timespecDiff(struct timespec *timeA_p, struct timespec *timeB_p)
{
  return ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
           ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec);
}

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

std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
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
    ( "file2,f2", po::value<std::string>(), "input file 2" )
//    ( "cfgequal", "cfgs equal via downward closure?" )
    ( "float", "float semiring" )
    ( "rexp", "commutative regular expression semiring" )
    ( "slset", "explicit semilinear sets semiring, no simplification " )
    ( "mlset", "abstraction over semilinear sets" )
    ( "vec-simpl", "vector simplification only (only semilinear and multilinear sets)" )
    ( "lin-simpl", "linear set simplification (only semilinear sets)" )
    ( "free", "free semiring" )
    //( "lossy", "lossy semiring" )
    ( "lossyC", "if used with --file: downward closure via Courcelle\nif used with --file and --file2: compare two input grammars for inequality via downward closure\nuse --refine n to fix refinement depth")
    ( "refine", po::value<int>(), "refinement depth for inequality check via downward closure; only available for --lossyC --file --file2")
   // ( "lossyComp", "")
//    ( "lossyArbTest", "")
//    ( "lossyIntersectTest", "lossy semiring intersection test" )
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
//      !vm.count("lossy") &&
//      !vm.count("lossyIntersectTest") &&
      !vm.count("lossyC") //&&
//      !vm.count("lossyComp") &&
//      !vm.count("lossyArbTest")) {
     ) {
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

  // see to the input of the second file
  std::string input_2_all = "";
  std::vector<std::string> input_2;
  if (vm.count("file2")) {
    // we are reading the input from the given second file
    std::ifstream file;
    file.open(vm["file2"].as<std::string>(), std::ifstream::in);
//    std::cout << "--file2 string:\t" << vm["file2"].as<std::string>() << std::endl;
    if (file.fail()) {
      std::cerr << "Could not open input file: " << vm["file2"].as<std::string>() << std::endl;
    } else {
        while (std::getline(file, line)) {
            input_2.push_back(line);
        }

        input_2_all = std::accumulate(input_2.begin(), input_2.end(), std::string(""));
    }
  }

  Parser p;



  int who = RUSAGE_SELF;
  struct rusage usage;
  int ret;
  struct timeval startT, endT;
  long mtime, seconds, useconds;



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

  } /*else if (vm.count("lossy")) {

//    struct timespec start, end;
    auto equations = p.lossy_fa_parser(input_all);
    if (equations.empty()) return EXIT_FAILURE;

    VarId S_1 = equations[0].first;

//              std::cout << "approximating language..." << std::endl;

    if (vm.count("file2")) {
        auto equations_2 = p.lossy_fa_parser(input_2_all);
        VarId S_2 = equations_2[0].first;

        auto approx_1 = LossyFiniteAutomaton::downwardClosureDerivationTrees(equations, S_1);
        auto approx_2 = LossyFiniteAutomaton::downwardClosureDerivationTrees(equations_2, S_2);

        bool A1_subset_A2 = approx_2.contains(approx_1);
        bool A2_subset_A1 = approx_1.contains(approx_2);

        if(A1_subset_A2 && A2_subset_A1) {
            std::cout << "Both languages have the same downward closure:\t" << approx_1.lossify().string() << std::endl;
        } else {
            LossyFiniteAutomaton L1_intersect_A2c = LossyFiniteAutomaton::null();
            bool L1_intersect_A2c_changed = false;
            LossyFiniteAutomaton L2_intersect_A1c = LossyFiniteAutomaton::null();
            bool L2_intersect_A1c_changed = false;

            if(!A1_subset_A2) {
                VarId startSymbol_1_2;
                auto intersectionGrammar = approx_2.complement().intersectionWithCFG(startSymbol_1_2, S_1, equations);
                std::queue<VarId> worklist;
                worklist.push(startSymbol_1_2);
                intersectionGrammar = NonCommutativePolynomial<LossyFiniteAutomaton>::cleanSystem(intersectionGrammar, worklist);

                L1_intersect_A2c = NonCommutativePolynomial<LossyFiniteAutomaton>::shortestWord
                        (intersectionGrammar, startSymbol_1_2);
                L1_intersect_A2c_changed = true;
            }

            if(!A2_subset_A1) {
                VarId startSymbol_2_1;
                auto intersectionGrammar_2 = approx_1.complement().intersectionWithCFG(startSymbol_2_1, S_2, equations_2);
                std::queue<VarId> worklist;
                worklist.push(startSymbol_2_1);
                intersectionGrammar_2 = NonCommutativePolynomial<LossyFiniteAutomaton>::cleanSystem(intersectionGrammar_2, worklist);

                L2_intersect_A1c = NonCommutativePolynomial<LossyFiniteAutomaton>::shortestWord
                        (intersectionGrammar_2, startSymbol_2_1);
                L2_intersect_A1c_changed = true;
            }

            if(L1_intersect_A2c_changed && L2_intersect_A1c_changed) {
                auto first = L1_intersect_A2c.string();
                auto second = L2_intersect_A1c.string();

                if(first.size() <= second.size()) {
                    std::cout << "There is a word in L1 that is not in L2." << std::endl;
                    std::cout << "Shortest such word:\t" << first << std::endl;
                } else {
                    std::cout << "There is a word in L2 that is not in L1." << std::endl;
                    std::cout << "Shortest such word:\t" << second << std::endl;
                }
            } else {
                if(L1_intersect_A2c_changed) {
                    std::cout << "There is a word in L1 that is not in L2." << std::endl;
                    std::cout << "Shortest such word:\t" << L1_intersect_A2c.string() << std::endl;
                } else {
                    std::cout << "There is a word in L2 that is not in L1." << std::endl;
                    std::cout << "Shortest such word:\t" << L2_intersect_A1c.string() << std::endl;
                }
            }
        }
    } else {
//        clock_gettime(CLOCK_MONOTONIC, &start);
//        gettimeofday(&startT, NULL);
        auto approximation = LossyFiniteAutomaton::downwardClosureDerivationTrees(equations, S_1);

//        gettimeofday(&endT, NULL);
//        clock_gettime(CLOCK_MONOTONIC, &end);

//        uint64_t timeElapsed = timespecDiff(&end, &start) / 1000000;
//        ret = getrusage(who, &usage);


//        std::string approx = approximation.string();


//        seconds  = endT.tv_sec  - startT.tv_sec;
//        useconds = endT.tv_usec - startT.tv_usec;
//
//        mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
//        unsigned long actualTime = (usage.ru_stime.tv_sec * 1000) + (usage.ru_stime.tv_usec /1000);
//        std::cout  /*<< "\t time:\t" << mtime  << "\tmemory used: " << usage.ru_maxrss << "KB\t" << vm["file"].as<std::string>() /*<< " closure:\t " << approx << std::endl;
    }
//                  std::cout << "Delossified regex (via string manipulation):\t" << LossyFiniteAutomaton::lossifiedRegex(approx) << std::endl;
//                  std::cout << "Delossified regex (via automaton):\t" << approximation.lossify().string() << std::endl;
  } */ else if(vm.count("lossyC")) {

      struct timespec start, end;
      auto equations = p.lossy_fa_parser(input_all);
      if (equations.empty()) return EXIT_FAILURE;

      VarId S_1 = equations[0].first;

  //              std::cout << "approximating language..." << std::endl;

      if (vm.count("file2")) {
          auto equations_2 = p.lossy_fa_parser(input_2_all);
          VarId S_2 = equations_2[0].first;

          int refinementDepth = 0;
          if(vm.count("refine")) {
              refinementDepth = vm["refine"].as<int>();
          }


          auto witness = LossyFiniteAutomaton::refineCourcelle(equations, S_1, equations_2, S_2, refinementDepth);

          if(witness != LossyFiniteAutomaton::null()) {
              std::cout << "Found witness: " << witness.string() << std::endl;
          }
//          std::cout << currentDateTime() << "\tstarting..." << std::endl;
//          auto approx_1 = LossyFiniteAutomaton::downwardClosureCourcelle(equations, S_1);
////          std::cout << currentDateTime() << "\tapprox 1 states:" << approx_1.size() << std::endl;
////          std::cout << "A1: " << approx_1.string() << "\t";
//          auto approx_2 = LossyFiniteAutomaton::downwardClosureCourcelle(equations_2, S_2);
////          std::cout << currentDateTime() << "\tapprox 2 states:" << approx_2.size() << std::endl;
////          std::cout << "A2: " << approx_2.string() << std::endl;
//
//          bool A1_subset_A2 = approx_2.contains(approx_1);
//          bool A2_subset_A1 = approx_1.contains(approx_2);
//          ret = getrusage(who, &usage);
//          auto memoryUsage = usage.ru_maxrss;
////          std::cout << memoryUsage << ",";
//          if(A1_subset_A2 && A2_subset_A1) {
//              std::cout << "equal"<< std::endl;
////              std::cout << "0,0,0,0,0" << std::endl;
////              std::cout << "Same closure: " << vm["file"].as<std::string>() << vm["file2"].as<std::string>() << std::endl;
//          } else {
//              Timer timer;
//              timer.Start();
//              int L21raw = 0, L12raw = 0, L21 = 0, L12 = 0;
//              LossyFiniteAutomaton L1_intersect_A2c = LossyFiniteAutomaton::null();
//              bool L1_intersect_A2c_changed = false;
//              LossyFiniteAutomaton L2_intersect_A1c = LossyFiniteAutomaton::null();
//              bool L2_intersect_A1c_changed = false;
//
//              if(!A1_subset_A2) {
//                  VarId startSymbol_1_2;
//                  auto A2c = approx_1.minus(approx_2);
////                  std::cout << "A2c size: " << A2c.size() << std::endl;
////                  std::cout << "A2c: " << A2c.string() << std::endl;
//                  auto intersectionGrammar = A2c.intersectionWithCFG(startSymbol_1_2, S_1, equations);
//                  L12raw = intersectionGrammar.size();
//                  std::queue<VarId> worklist;
//                  worklist.push(startSymbol_1_2);
//                  intersectionGrammar = NonCommutativePolynomial<LossyFiniteAutomaton>::cleanSystem(intersectionGrammar, worklist);
//                  L12 = intersectionGrammar.size();
////                  std::cout << "clean intersection grammar 1 x A2c:" << std::endl;
////                  for(auto &equation: intersectionGrammar) {
////                      std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
////                  }
//
//                  L1_intersect_A2c = NonCommutativePolynomial<LossyFiniteAutomaton>::shortestWord
//                          (intersectionGrammar, startSymbol_1_2);
//                  L1_intersect_A2c_changed = true;
//              }
//
//              if(!A2_subset_A1) {
//                  VarId startSymbol_2_1;
//                  auto A1c = approx_2.minus(approx_1);
////                  std::cout << "A1c size: " << A1c.size() << std::endl;
////                  std::cout << "A1c: " << A1c.string() << std::endl;
//                  auto intersectionGrammar_2 = A1c.intersectionWithCFG(startSymbol_2_1, S_2, equations_2);
//                  L21raw = intersectionGrammar_2.size();
//                  std::queue<VarId> worklist;
//                  worklist.push(startSymbol_2_1);
//                  intersectionGrammar_2 = NonCommutativePolynomial<LossyFiniteAutomaton>::cleanSystem(intersectionGrammar_2, worklist);
//                  L21 = intersectionGrammar_2.size();
////                  std::cout << "clean intersection grammar 2 x A1c:" << std::endl;
////                  for(auto &equation: intersectionGrammar_2) {
////                      std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
////                  }
//
//                  L2_intersect_A1c = NonCommutativePolynomial<LossyFiniteAutomaton>::shortestWord
//                          (intersectionGrammar_2, startSymbol_2_1);
//                  L2_intersect_A1c_changed = true;
//              }
//              timer.Stop();
//              ret = getrusage(who, &usage);
//              auto memoryUsage = usage.ru_maxrss;
////              std::cout << memoryUsage << ",";
//
//              std::cout << "different" << std::endl;
//              if(L1_intersect_A2c_changed && L2_intersect_A1c_changed) {
//                  auto first = L1_intersect_A2c.string();
//                  auto second = L2_intersect_A1c.string();
//
//
//                  if(first.size() <= second.size()) {
////                      std::cout << timer.GetMilliseconds().count() << "," << L12raw << "," << L12 << "," << first.size() << std::endl;
////                      std::cout << vm["file"].as<std::string>() << vm["file2"].as<std::string>() << " Shortest word in L1 \\ L2: " << first << std::endl;
////                      std::cout << "length of shortest word: " << first.size() << std::endl;
//                  } else {
////                      std::cout << vm["file"].as<std::string>() << vm["file2"].as<std::string>() << " Shortest word in L2 \\ L1: " << second << std::endl;
////                      std::cout << "length of shortest word: " << second.size() << std::endl;
////                      std::cout << timer.GetMilliseconds().count() << "," << L21raw << "," << L21 << "," << second.size() << std::endl;
////                      std::cout << "There is a word in L2 that is not in L1." << std::endl;
////                      std::cout << "Shortest such word:\t" << second << std::endl;
//                  }
//              } else {
////                  std::cout << currentDateTime() << "\tfound witness for difference" << std::endl;
//                  if(L1_intersect_A2c_changed) {
////                      std::cout << vm["file"].as<std::string>() << vm["file2"].as<std::string>() << " Shortest word in L1 \\ L2: " << L1_intersect_A2c.string() << std::endl;
////                      std::cout << "length of shortest word: " << L1_intersect_A2c.string().size() << std::endl;
////                      std::cout << "witness:\t" << L1_intersect_A2c.string() << std::endl;
////                      std::cout << timer.GetMilliseconds().count() << "," << L12raw << "," << L12 << "," << L1_intersect_A2c.string().size() << std::endl;
////                      std::cout << "There is a word in L1 that is not in L2." << std::endl;
////                      std::cout << "Shortest such word:\t" << L1_intersect_A2c.string() << std::endl;
//                  } else {
////                      std::cout << vm["file"].as<std::string>() << vm["file2"].as<std::string>() << " Shortest word in L1 \\ L2: " << L2_intersect_A1c.string() << std::endl;
////                      std::cout << "witness:\t" << L2_intersect_A1c.string() << std::endl;
////                      std::cout << timer.GetMilliseconds().count() << "," << L21raw << "," << L21 << "," << L2_intersect_A1c.string().size() << std::endl;
////                      std::cout << "length of shortest word: " << L2_intersect_A1c.string().size() << std::endl;
////                      std::cout << "There is a word in L2 that is not in L1." << std::endl;
////                      std::cout << "Shortest such word:\t" << L2_intersect_A1c.string() << std::endl;
//                  }
//              }
//          }
      } else {
//          clock_gettime(CLOCK_MONOTONIC, &start);
//          gettimeofday(&startT, NULL);

//          clock_t startTime = clock();
          auto approximation = LossyFiniteAutomaton::downwardClosureCourcelle(equations, S_1);
//          gettimeofday(&endT, NULL);
//          clock_t endTime = clock();
//          clock_gettime(CLOCK_MONOTONIC, &end);
//          auto totalTime = (endTime - startTime)/* / (CLOCKS_PER_SEC / 1000)*/;

//          ret = getrusage(who, &usage);
//          uint64_t timeElapsed = timespecDiff(&end, &start) / 1000000;
          std::string approx = approximation.string();
//          unsigned long actualTime = (usage.ru_stime.tv_sec * 1000) + (usage.ru_stime.tv_usec /1000);
//          seconds  = endT.tv_sec  - startT.tv_sec;
//          useconds = endT.tv_usec - startT.tv_usec;
//
//          mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
//          std::cout << "system:" << std::endl;
//          for(auto &equation: equations) {
//              std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
//          }
          std::cout /* << "\t\t\t" << "\t time:\t" << mtime << "\tmemory used: " << usage.ru_maxrss << "KB\t" *//*<< vm["file"].as<std::string>() */<< /*"\nclosure:\t" <<*/ approx << /*"\n\n\n"<<*/ std::endl;
      }
  } else if(vm.count("lossyIntersectTest")) {
      auto equations = p.lossy_fa_parser(input_all);
      if (equations.empty()) return EXIT_FAILURE;
      LossyFiniteAutomaton::intersectionTest(equations);
  } /*else if(vm.count("lossyComp")) {

      struct timespec start, end;
      auto equations = p.lossy_fa_parser(input_all);
      if (equations.empty()) return EXIT_FAILURE;

      VarId S_1 = equations[0].first;
      auto derivationTreeClosure = LossyFiniteAutomaton::downwardClosureDerivationTrees(equations, S_1).minimize();
      auto courcelleClosure = LossyFiniteAutomaton::downwardClosureCourcelle(equations, S_1).minimize();

//      FILE * file1;
//      file1 = fopen ("courcelle.dot","w");
//      courcelleClosure.write_dot_file(file1);
//      fclose (file1);
//
//      FILE * file2;
//      file2 = fopen ("derivation_tree.dot","w");
//      derivationTreeClosure.write_dot_file(file2);
//      fclose (file2);

      if(!(derivationTreeClosure == courcelleClosure)) {
          std::cout << "Closures not equal" << std::endl; //, grammar:\t" << std::endl; // << vm["file"].as<std::string>() << "\t size: " << equations.size();
//          for(auto &equation: equations) {
//              std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
//          }
          std::cout << "Derivation tree closure: " << derivationTreeClosure.string() << std::endl;
          std::cout << "Courcelle closure: " << courcelleClosure.string() << std::endl;
      } else {
          std::cout << vm["file"].as<std::string>() << "\t size: " << equations.size() << "\t same results" << std::endl;
//          std::cout << "closure derivation trees: " << derivationTreeClosure.lossify().string() << std::endl;
//          std::cout << "closure courcelle: " << courcelleClosure.string() << std::endl;
//          std::cout << "grammar:" << std::endl;
//          for(auto &equation: equations) {
//              std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//          }
      }
  }  else if(vm.count("lossyArbTest")) {
      LossyFiniteAutomaton der = LossyFiniteAutomaton("(b*a[ab]|b*[ab])(c|())|b*ac|b*c|b*a|b*");
      LossyFiniteAutomaton cour = LossyFiniteAutomaton("(b*a[ab]|b*[ab])(c|())|b*a(c|())|b*c|b*");
      std::cout << "same language:\t" << (cour == der) << std::endl;

  }*/ else if (vm.count("prefix")) {

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
