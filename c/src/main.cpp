#include <iostream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/program_options.hpp>
#include "float-semiring.h"
#include "free-semiring.h"
#include "semilinSetExp.h"
#include "matrix.h"
#include "polynomial.h"
#include "newton.h"
#include "commutativeRExp.h"
#include "parser.h"

template <typename SR>
struct VertexProp {
	std::string name;   // used for graphviz output
	VarPtr var;         // var and rex combines the equations in the vertex
	Polynomial<SR> rex;
};

// group the equations to SCCs
template <typename SR>
std::vector<std::vector<std::pair<VarPtr, Polynomial<SR>>>> group_by_scc(std::vector<std::pair<VarPtr, Polynomial<SR>>> equations, bool graphviz_output)
{
	// create map of variables to [0..n]. this is used to enumerate important variables in a clean way from 0 to n during graph construction
	std::map<VarPtr, int> var_key;

	// build the graph
	boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProp<SR>> graph(equations.size());
	for(auto e_it = equations.begin(); e_it != equations.end(); ++e_it)
	{
		// if the variable is not yet in the map insert it together with the size of
		// the map. this way we get a unique identifier for the variable counting from 0 to n
		if(var_key.find(e_it->first) == var_key.end())
			var_key.insert(var_key.begin(), std::pair<VarPtr,int>(e_it->first,var_key.size()));
		int a = var_key.find(e_it->first)->second; // variable key
		graph[a].var = e_it->first; // store VarPtr to the vertex
		graph[a].name = e_it->first->string(); // store the name to the vertex
		graph[a].rex = e_it->second; // store the regular expression to the vertex

		auto v = e_it->second.get_variables(); // all variables of this rule;
		for(auto v_it = v.begin(); v_it != v.end(); ++v_it)
		{
			if(var_key.find(*v_it) == var_key.end()) // variable is not yet in the map
				var_key.insert(var_key.begin(), std::pair<VarPtr,int>(*v_it,var_key.size()));
			int b = var_key.find(*v_it)->second; // variable key
			boost::add_edge(a, b, graph);
		}
	} // graph is complete

	if(graphviz_output)
	{
		// output the created graph to graph.dot
		boost::dynamic_properties dp;
		dp.property("label", boost::get(&VertexProp<SR>::name, graph)); // vertex name is the name of the equation
		dp.property("node_id", get(boost::vertex_index, graph)); // this is needed
		std::ofstream outf("graph.dot");
		boost::write_graphviz_dp(outf, graph, dp);
	}

	// calculate strong connected components and store them in 'component'
	std::vector<int> component(boost::num_vertices(graph));
	boost::strong_components(graph,&component[0]);

	// group neccessary equations together
	int num_comp = *std::max_element(component.begin(), component.end()) + 1; // find the number of components
	std::vector<std::vector<std::pair<VarPtr,Polynomial<SR>>>> grouped_equations;
	grouped_equations.resize(num_comp);

	// iterate over all vertices (0 to n)
	// collect the necessary variables + equations for every component
	for (unsigned int j = 0; j != component.size(); ++j)
	{
		//std::cout << j << ", " << graph[j].var << " is in component " << component[j] << std::endl;
		grouped_equations[component[j]].push_back(std::pair<VarPtr,Polynomial<SR>>(graph[j].var, graph[j].rex));
	}

	return grouped_equations;
}

// apply the newton method to the given input
template <typename SR>
std::map<VarPtr, SR> apply_newton(std::vector<std::pair<VarPtr, Polynomial<SR>>> equations, bool scc, bool iteration_flag, int iterations, bool graphviz_output)
{
	// TODO: sanity checks on the input!

	// generate an instance of the newton solver
	Newton<SR> newton;

	// if we use the scc method, group the equations
	// the outer vector contains SCCs starting with a bottom SCC at 0
	std::vector<std::vector<std::pair<VarPtr,Polynomial<SR>>>> equations2;
	if(scc)
	{
		auto tmp = group_by_scc(equations, graphviz_output);
		equations2.insert(equations2.begin(), tmp.begin(), tmp.end());
	}
	else if (!scc)
	{
		equations2.push_back(equations);
	}

	// this holds the solution
	std::map<VarPtr, SR> solution;

	// the same loop is used for both the scc and the non-scc variant
	// in the non-scc variant, we just run once through the loop
	//for(auto it1 = equations2.begin(); it != equations2.end(); ++it)
	for(unsigned int j = 0; j != equations2.size(); ++j)
	{
		// use the solutions to get rid of variables in the remaining equations
		// does nothing in the first round
		std::vector<std::pair<VarPtr, Polynomial<SR>>> tmp1;
		for(auto it = equations2[j].begin(); it != equations2[j].end(); ++it)
		{ // it = (VarPtr, Polynomial[SR])
			auto tmp2 = it->second.partial_eval(solution);
			tmp1.push_back(std::pair<VarPtr, Polynomial<SR>>(it->first, tmp2));
		}
		// replace old equations with simplified ones
		equations2[j] = tmp1;

		// dynamic iterations
		if(!iteration_flag)
			iterations = equations2[j].size();

		// do some real work here
		std::map<VarPtr, SR> result = newton.solve_fixpoint(equations2[j], iterations);

		// copy the results into the solution map
		solution.insert(result.begin(), result.end());
	}

	return solution;
}

template <typename SR>
std::string result_string(std::map<VarPtr,SR> result)
{
	std::stringstream ss;
	for(auto r_it = result.begin(); r_it != result.end(); ++r_it)
		ss << r_it->first << " == " << r_it->second << std::endl;
	return ss.str();
}

int main(int argc, char* argv[])
{
	namespace po = boost::program_options;

	po::options_description desc("Allowed options");
	desc.add_options()
		( "scc", "apply newton method iteratively to strongly connected components of the equation graph" )
		( "help,h", "print this help message" )
		( "iterations,i", po::value<int>(), "specify the number of newton iterations. default is optimal number" )
		//( "verbose", "enable verbose output" )
		//( "debug", "enable debug output" )
		( "file,f", po::value<std::string>(), "input file" )
		( "float", "float semiring" )
		( "rexp", "commutative regular expression semiring" )
		( "graphviz", "create the file graph.dot with the equation graph" )
		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if(vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 0;
	}

	int iterations=0;
	if(vm.count("iterations"))
		iterations = vm["iterations"].as<int>();

	// check if we can do something useful
	if(!vm.count("float") && !vm.count("rexp")) // check for all compatible parameters
		return 0;


	std::vector<std::string> input;
	std::string line;
	if(vm.count("file"))
	{
		// we are reading the input from the given file
		std::ifstream file;
		file.open(vm["file"].as<std::string>(), std::ifstream::in);
		while(std::getline(file, line))
			input.push_back(line);
	}
	else
	{
		// we are reading from stdin
		while(std::getline(std::cin, line))
			input.push_back(line);
	}

	// join the input into one string
	std::string input_all = std::accumulate(input.begin(), input.end(), std::string(""));

	Parser p;

	if(vm.count("rexp"))
	{
		// parse the input into a list of (Var → Polynomial[SR])
		std::vector<std::pair<VarPtr, Polynomial<CommutativeRExp>>> equations(p.rexp_parser(input_all));
		if(equations.empty()) return -1;

		for(auto eq_it = equations.begin(); eq_it != equations.end(); ++eq_it)
		{
			std::cout << "* " << eq_it->first << " → " << eq_it->second << std::endl;
		}

		// apply the newton method to the equations
		auto result = apply_newton<CommutativeRExp>(equations, vm.count("scc"), vm.count("iterations"), iterations, vm.count("graphviz"));
		std::cout << result_string(result) << std::endl;
	}
	else if(vm.count("float"))
	{
		std::vector<std::pair<VarPtr, Polynomial<FloatSemiring>>> equations(p.float_parser(input_all));
		if(equations.empty()) return -1;

		for(auto eq_it = equations.begin(); eq_it != equations.end(); ++eq_it)
		{
			std::cout << "* " << eq_it->first << " → " << eq_it->second << std::endl;
		}

		auto result = apply_newton<FloatSemiring>(equations, vm.count("scc"), vm.count("iterations"), iterations, vm.count("graphviz"));
		std::cout << result_string(result) << std::endl;
	}

	return 0;
}
