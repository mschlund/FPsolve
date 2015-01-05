/*
 * solver_utils.h
 *
 *  Created on: 20.12.2014
 *      Author: schlund
 */

#ifndef SOLVER_UTILS_H_
#define SOLVER_UTILS_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graphviz.hpp>

#include "../datastructs/equations.h"

#include "../polynomials/commutative_polynomial.h"
#include "../utils/timer.h"


template <typename SR>
struct VertexProp {
  std::string name;   // used for graphviz output
  VarId var;         // var and rex combines the equations in the vertex
  CommutativePolynomial<SR> rex;
};

// group the equations by SCCs
template <typename SR, template <typename> class Poly>
std::vector< GenericEquations<Poly, SR> >
group_by_scc(const GenericEquations<Poly, SR> &equations,
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
  std::vector< std::vector< std::pair< VarId,CommutativePolynomial<SR> > > >
    grouped_equations(num_comp);
  // grouped_equations.resize(num_comp);

  // iterate over all vertices (0 to n)
  // collect the necessary variables + equations for every component
  for (std::size_t j = 0; j != component.size(); ++j) {
    //std::cout << j << ", " << graph[j].var << " is in component " << component[j] << std::endl;
    grouped_equations[component[j]].push_back(std::pair<VarId,CommutativePolynomial<SR>>(graph[j].var, graph[j].rex));
  }

  return grouped_equations;
}

// apply solving method to the given input
template <template <typename> class SolverType,
          template <typename> class Poly,
          typename SR>
ValuationMap<SR> apply_solver(
    const GenericEquations<Poly, SR> &equations,
    bool scc, bool iteration_flag, std::size_t iterations, bool graphviz_output) {

  // TODO: sanity checks on the input!

  // generate an instance of the solver
  SolverType<SR> solver;

  // if we use the scc method, group the equations
  // the outer vector contains SCCs starting with a bottom SCC at 0
  std::vector<GenericEquations<Poly, SR>> equations2;
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
    GenericEquations<Poly, SR> tmp1;
    for (auto it = equations2[j].begin(); it != equations2[j].end(); ++it)
    { // it = (VarId, Polynomial[SR])
      auto tmp2 = it->second.partial_eval(solution);
      tmp1.push_back(std::pair<VarId, Poly<SR>>(it->first, tmp2));
    }
    // replace old equations with simplified ones
    equations2[j] = tmp1;

    // dynamic iterations
    if (!iteration_flag) {
      // for commutative and idempotent SRs Newton has converged after n+1 iterations, so use this number as default
      iterations = equations2[j].size() + 1;
    }

    //std::cout << "Iterations: " << iterations << std::endl;

    // do some real work here
    ValuationMap<SR> result = solver.solve_fixpoint(equations2[j], iterations);

    // copy the results into the solution map
    solution.insert(result.begin(), result.end());
  }

  timer.Stop();
  std::cout
      << "Solving time:\t" << timer.GetMilliseconds().count()
      << " ms" << " ("
    << timer.GetMicroseconds().count()
    << "us)" << std::endl;



  return solution;
}


#endif /* SOLVER_UTILS_H_ */
