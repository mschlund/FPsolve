#include <iostream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graphviz.hpp>
#include "float-semiring.h"
#include "free-semiring.h"
#include "semilinSetExp.h"
#include "matrix.h"
#include "polynomial.h"
#include "newton.h"
#include "commutativeRExp.h"
#include "parser.h"

void test_newton()
{
	/*
	Newton<FloatSemiring> newton;
	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("x"));
	variables.push_back(Var::getVar("y"));
	variables.push_back(Var::getVar("z"));
	std::cout << "- newton (float):" << std::endl;
	std::vector<Polynomial<FloatSemiring> > polynomials = get_newton_test_polynomials();
	Matrix<FloatSemiring> result = newton.solve_fixpoint(polynomials, variables, 2);
	std::cout << result << std::endl;

*/


/*	Newton<SemilinSetExp> newton;
	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("x1"));
	variables.push_back(Var::getVar("x2"));
	variables.push_back(Var::getVar("x3"));
	variables.push_back(Var::getVar("x4"));
	std::cout << "- newton (cnt-SR):" << std::endl;

	std::vector<Polynomial<SemilinSetExp> > polynomials;
	Polynomial<SemilinSetExp> f1 = Polynomial<SemilinSetExp>({
		Monomial<SemilinSetExp>(SemilinSetExp::one(), {Var::getVar("x1"),Var::getVar("x1")}),
		Monomial<SemilinSetExp>(SemilinSetExp::one(), {Var::getVar("x2"),Var::getVar("x3")})
		});
	Polynomial<SemilinSetExp> f2 = Polynomial<SemilinSetExp>({
		Monomial<SemilinSetExp>(SemilinSetExp::one(), {Var::getVar("x1"),Var::getVar("x2")}),
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("a")), {})
		});
	Polynomial<SemilinSetExp> f3 = Polynomial<SemilinSetExp>({
			Monomial<SemilinSetExp>(SemilinSetExp::one(), {Var::getVar("x4"),Var::getVar("x3")}),
			Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("b")), {})
		});
	Polynomial<SemilinSetExp> f4 = Polynomial<SemilinSetExp>({
		Monomial<SemilinSetExp>(SemilinSetExp::one(), {Var::getVar("x4"),Var::getVar("x4")}),
		Monomial<SemilinSetExp>(SemilinSetExp::one(), {Var::getVar("x2"),Var::getVar("x3")})
		});

	polynomials.push_back(f1);
	polynomials.push_back(f2);
	polynomials.push_back(f3);
	polynomials.push_back(f4);

	Matrix<SemilinSetExp> result = newton.solve_fixpoint(polynomials, variables, 2);
	std::cout << result << std::endl;
*/

	Newton<CommutativeRExp> newton;
	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("x"));
	variables.push_back(Var::getVar("y"));
	variables.push_back(Var::getVar("z"));
	std::cout << "- newton (counting-SR):" << std::endl;

	std::vector<Polynomial<CommutativeRExp> > polynomials;
	// define new polynomial axy+b
	Polynomial<CommutativeRExp> f1 = Polynomial<CommutativeRExp>({
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("a")), {Var::getVar("x"),Var::getVar("y")}),
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("b")), {}) });

	// define new polynomial cyz+dyx+e
	Polynomial<CommutativeRExp> f2 = Polynomial<CommutativeRExp>({
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("c")), {Var::getVar("y"),Var::getVar("z")}),
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("d")), {Var::getVar("y"),Var::getVar("x")}),
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("e")), {}) });

	// define new polynomial fx+g
	Polynomial<CommutativeRExp> f3 =  Polynomial<CommutativeRExp>({
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("f")), {Var::getVar("x")}),
		Monomial<CommutativeRExp>(CommutativeRExp(Var::getVar("g")), {}) });
	polynomials.push_back(f1);
	polynomials.push_back(f2);
	polynomials.push_back(f3);

	Matrix<CommutativeRExp> result = newton.solve_fixpoint(polynomials, variables, 4000);
	std::cout << "done!" << std::endl;
	//std::cout << result << std::endl;


/*	Newton<SemilinSetExp> newton;
	std::vector<VarPtr> v;
	v.push_back(Var::getVar("x"));
	std::cout << "- newton (cnt-SR):" << std::endl;

	std::vector<Polynomial<SemilinSetExp> > polys;
	Polynomial<SemilinSetExp> f0 = Polynomial<SemilinSetExp>({
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("s")), {Var::getVar("x"),Var::getVar("x")}),
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("s")), {}) });

	polys.push_back(f0);
	Matrix<SemilinSetExp> r = newton.solve_fixpoint(polys, v, 2);
	std::vector<SemilinSetExp> sol1 = r.getElements();

	std::cout << sol1[0] << std::endl;

//	std::cout << (SemilinSetExp(Var::getVar("r"))*SemilinSetExp(Var::getVar("s"))*sol1[0])<< std::endl;

	std::vector<VarPtr> variables;
	variables.push_back(Var::getVar("e"));
	std::vector<Polynomial<SemilinSetExp> > polynomials;
	Polynomial<SemilinSetExp> f1 = Polynomial<SemilinSetExp>({
		Monomial<SemilinSetExp>(SemilinSetExp::one(), {}),
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("l"))*SemilinSetExp(Var::getVar("s")), {Var::getVar("e")}),
		Monomial<SemilinSetExp>(SemilinSetExp(Var::getVar("r"))*SemilinSetExp(Var::getVar("s"))*sol1[0], {Var::getVar("e")}) });

	polynomials.push_back(f1);

	Matrix<SemilinSetExp> result = newton.solve_fixpoint(polynomials, variables, 1);
	std::cout << result << std::endl;
	std::cout << (SemilinSetExp(Var::getVar("l"))*SemilinSetExp(Var::getVar("s")) + SemilinSetExp(Var::getVar("r"))*SemilinSetExp(Var::getVar("s"))*sol1[0]).star() << std::endl;
*/

/*
	SemilinSetExp x = SemilinSetExp(Var::getVar("a"));
	SemilinSetExp y = SemilinSetExp(Var::getVar("b"));
	SemilinSetExp z = (x*x*x+y*y*y).star();
	std::cout << z << std::endl;
*/

}

int main(int argc, char* argv[])
{
	if(argc == 1) // no command-line arguments
		test_newton();
	else // there is an argument
	{
		Parser p;
		std::string input;
		if(std::string("-float").compare(argv[1]) == 0)
		{
			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_float(input);
				std::cout << result << std::endl;
			}
		}
		else if(std::string("-free").compare(argv[1]) == 0)
		{
			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_free(input);
				std::cout << result << std::endl;
			}
		}
		else if(std::string("-rexp").compare(argv[1]) == 0)
		{
			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_rexp(input);
				std::cout << result << std::endl;
			}
		}
		else if(std::string("-polyrexp").compare(argv[1]) == 0)
		{
			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_polyrexp(input);
				std::cout << result << std::endl;
			}
		}
		else if(std::string("-grammar").compare(argv[1]) == 0)
		{
			Newton<CommutativeRExp> newton;

			// accumulate the input in rules
			std::vector<Polynomial<CommutativeRExp> > rules;
			// and all the variables in vars
			std::vector<VarPtr> vars;

			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_grammar(input); // (string, polynomial)
				rules.push_back(result.second);
				vars.push_back(Var::getVar(result.first));

				std::cout << result.first << " → " << result.second << std::endl;
			}
			int i = 2; // default value
			if(argc > 2) // set the given iteration count
				i = std::atoi(argv[2]);

			Matrix<CommutativeRExp> result = newton.solve_fixpoint(rules, vars, i);
			std::cout << i << "-th newton iteration: " << std::endl;
			std::cout << result << std::endl;
		}
		else if(std::string("-grammar_scc").compare(argv[1]) == 0)
		{
			Newton<CommutativeRExp> newton;

			// accumulate all the rules with the variables
			std::vector<std::pair<VarPtr, Polynomial<CommutativeRExp>>> rules;

			while(std::cout << "> " && std::getline(std::cin, input))
			{
				auto result = p.parse_grammar(input); // (string, polynomial)
				rules.push_back(result);
				std::cout << result.first << " → " << result.second << std::endl;
			}
			int i = 2; // default value
			if(argc > 2) // set the given iteration count
				i = std::atoi(argv[2]);

			// create map of variables to [0..n]. this is used to enumerate important variables in a clean way from 0 to n
			std::map<VarPtr, int> var_key;
			// create graph of variable dependencies
			std::cout << "rules size: " << rules.size() << std::endl;

			struct VertexProp {
				std::string name;
				VarPtr var;
				Polynomial<CommutativeRExp> rex;
			};

			boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexProp> graph(rules.size());
			for(auto r_it = rules.begin(); r_it != rules.end(); ++r_it)
			{
				if(var_key.find(r_it->first) == var_key.end()) // variable is not yet in the map
					var_key.insert(var_key.begin(), std::pair<VarPtr,int>(r_it->first,var_key.size()));
				int a = var_key.find(r_it->first)->second; // variable key
				graph[a].var = r_it->first; // store VarPtr to the vertex
				graph[a].name = r_it->first->string(); // store the name to the vertex
				graph[a].rex = r_it->second; // store the regular expression to the vertex

				auto v = r_it->second.get_variables(); // all variables of this rule;
				for(auto v_it = v.begin(); v_it != v.end(); ++v_it)
				{
					if(var_key.find(*v_it) == var_key.end()) // variable is not yet in the map
						var_key.insert(var_key.begin(), std::pair<VarPtr,int>(*v_it,var_key.size()));
					int b = var_key.find(*v_it)->second; // variable key
					boost::add_edge(a, b, graph);
				}
			}

			// output the created graph to graph.dot
			boost::dynamic_properties dp;
			dp.property("label", boost::get(&VertexProp::name, graph)); // vertex name is the name of the equation
			dp.property("node_id", get(boost::vertex_index, graph)); // this is needed
			std::ofstream outf("graph.dot");
			boost::write_graphviz_dp(outf, graph, dp);

			// calculate strong connected components and store them in 'component'
			std::vector<int> component(boost::num_vertices(graph));
			boost::strong_components(graph,&component[0]);

			// group neccessary equations together
			int num_comp = *std::max_element(component.begin(), component.end()) + 1; // find the number of components
			std::vector<std::vector<VarPtr>> vars;
			std::vector<std::vector<Polynomial<CommutativeRExp>>> rexs;
			//std::cout << "number of components = " << num_comp << std::endl;
			vars.resize(num_comp);
			rexs.resize(num_comp);

			// iterate over all vertices (0 to n)
			// collect the necessary variables + equations for every component
			for (int j = 0; j != component.size(); ++j)
			{
				//std::cout << j << ", " << graph[j].var << " is in component " << component[j] << std::endl;
				vars[component[j]].push_back(graph[j].var);
				rexs[component[j]].push_back(graph[j].rex);
			}

			// calculate the solution
			std::map<VarPtr, CommutativeRExp> solution;
			for(int j = 0; j != num_comp; ++j)
			{
				// use the solutions to get rid of variables in the remaining equations
				std::vector<Polynomial<CommutativeRExp>> tmp1;
				for(auto r_it = rexs[j].begin(); r_it != rexs[j].end(); ++r_it)
				{
					auto tmp2 = r_it->partial_eval(solution);
					tmp1.push_back(tmp2);
				}
				// replace old equations with simplified ones
				rexs[j] = tmp1;

				// do some real work here
				Matrix<CommutativeRExp> result = newton.solve_fixpoint(rexs[j], vars[j], i);
				std::cout << i << "-th newton iteration: " << std::endl;
				// save the result in the valuation map
				std::vector<CommutativeRExp> result_vec = result.getElements();
				auto r_it = result_vec.begin();
				for(auto v_it = vars[j].begin(); v_it != vars[j].end(); ++v_it)
				{
					std::cout << *v_it << " == " << *r_it << std::endl;
					solution.insert(solution.begin(), std::pair<VarPtr, CommutativeRExp>(*v_it, *r_it));
					++r_it;
				}
			}

			// solution now holds the solution

		}
		else
		{
			std::cout << "unknown argument" << std::endl;
		}
	}

	return 0;
}
