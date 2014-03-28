#ifndef LOSSY_SEMIRING_H_
#define LOSSY_SEMIRING_H_

#include <string>
#include <memory>
#include <unordered_map>
#include <queue>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graphviz.hpp>

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/free-structure.h"
#include "../polynomials/non_commutative_polynomial.h"

#include "semiring.h"

template<typename LSR>
class Evaluator;

template <typename LSR>
struct Vertex {
    VarId var;         // var and rex combines the equations in the vertex
    NonCommutativePolynomial<LSR> rex;
};

template <typename LSR>
class LossySemiring: public StarableSemiring<LSR,
		Commutativity::NonCommutative, Idempotence::Idempotent> {
public:
    virtual ~LossySemiring(){};

    /*
     * Solves a polynomial system over the lossy semiring. For derivation of
     * the algorithm and proof of correctness, see "Esparza, Kiefer, Luttenberger:
     * Derivation Tree Analysis for Accelerated Fixed-Point Computation".
     *
     * If the system is not clean, i.e. there are polynomials in it with a fixpoint of 0, then you must set
     * systemNeedsCleaning to true in order for the algorithm to work as intended.
     *
     * If the system is not in QNF, then you need to set systemNeedsQNF to true.
     */
    static std::map<VarId, LSR> solvePolynomialSystem(
        std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations, bool systemNeedsCleaning, bool systemNeedsQNF) {

        std::map<VarId, LSR> solution;

//        auto components = group_by_scc(equations);
//        for(int i = 0; i < components.size(); i++) {
//            std::cout << "component " << i << ":\t" << std::endl;
//
//            for(auto &equation: components[i]) {
//                std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//            }
//        }

        // start out with a clean system
        if(systemNeedsCleaning) {
            std::queue<VarId> worklist;
            equations = cleanSystem(equations, worklist);
        }

        // bring the clean system into qnf, if necessary
        std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> qnf;

        if(systemNeedsQNF) {
            qnf = lossyQuadraticNormalForm(equations);
        } else {
            qnf = equations;
        }

        // build a zero vector for the evaluation of f
        std::map<VarId, LSR> zeroSystem;
        for(auto &equation: qnf) {
            zeroSystem.insert(std::make_pair(equation.first, LSR::null()));
        }

        // build the vectors f(0) and f^n(0), where n is the number of variables in the system
        int times = 1;
        std::map<VarId, LSR> f_0 = evaluateSystem(qnf, times, zeroSystem);
        times = qnf.size();
        std::map<VarId, LSR> f_n_0 = evaluateSystem(qnf, times, zeroSystem);

//        std::cout << "f(0):" << std::endl;
//        for(auto &mapping: f_0) {
//            std::cout << Var::GetVar(mapping.first).string() + ":\t\t" + mapping.second.string() << std::endl;
//        }
//        std::cout << "f^n(0):" << std::endl;
//        for(auto &mapping: f_n_0) {
//            std::cout << Var::GetVar(mapping.first).string() + ":\t\t" + mapping.second.string() << std::endl;
//        }

        // find the LossySemiring element in the "middle" of the expression
        LSR middle = LSR::null();
        for(auto &elem_mapping: f_0) {
            middle += elem_mapping.second;
        }

        // build the differential of the system
        std::map<VarId, NonCommutativePolynomial<LSR>> differential;
        for(auto &equation: qnf) {
            differential.insert(std::make_pair(equation.first, equation.second.differential_at(f_n_0)));
        }

        // sum all polynomials of the differential of the system; we don't need attribution
        // of each polynomial to the respective variable since we only care about the
        // leading and trailing coefficients of each monomial
        NonCommutativePolynomial<LSR> differential_sum = NonCommutativePolynomial<LSR>::null();
        for(auto &equation: differential) {
            differential_sum += equation.second;
        }

        // get the lefthand and righthand semiring element of the fixpoint
        LSR lefthandSum = differential_sum.getSumOfLeadingFactors();
        LSR righthandSum = differential_sum.getSumOfTrailingFactors();

        // fixpoint element
        LSR fixpoint = lefthandSum.star() * middle * righthandSum.star();

        for(auto &variable_mapping: qnf) {
            solution.insert(std::make_pair(variable_mapping.first, fixpoint));
        }

        return solution;
    }

    /*
     * Approximates the language defined by the grammar starting at S given in the system of equations by calculating
     * its downward closure.
     */
    static LSR downwardClosureDerivationTrees(const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations,
            const VarId &S) {

        std::queue<VarId> worklist;
        worklist.push(S);
        auto cleanEquations = cleanSystem(equations, worklist);

        // if the start symbol is unproductive, then it doesn't generate anything;
        // the downward closure of the empty set is the empty set
        bool SisProductive = false;
        for(auto &equation: cleanEquations) {
            if(equation.first == S) {
                SisProductive = true;
            }
        }

        if(!SisProductive) {
            return LSR::null();
        }

        auto cleanQNF = lossyQuadraticNormalForm(cleanEquations);
        auto components = group_by_scc(cleanQNF);

        std::map<VarId, LSR> knownValuations;
        for(int i = 0; i < components.size(); i++) {

            // if this is not the bottom component, substitute known valuations
            if(i != 0) {
                for(auto &equation: components[i]) {
                    equation.second = equation.second.partial_eval(knownValuations);
                }
            }

            auto newValuations = solvePolynomialSystem(components[i], false, false);

            for(auto &value: newValuations) {
                knownValuations.insert(std::make_pair(value.first, value.second));
            }
        }

        return knownValuations[S];
    }

private:

    /*
     * Brings a system into quadratic normal form.
     */
    static std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> lossyQuadraticNormalForm
            (const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations) {
            std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> systemBeingProcessed;

            // for all nonterminal productions, introduce new variables for the terminal factors that
            // appear in them as one would when calculating the CNF of a CFG; do not do this for linear terms since
            // QNF allows for linear terms
            std::map<VarId, LSR> variablesToConstants;
            systemBeingProcessed = NonCommutativePolynomial<LSR>::eliminateTerminalsInNonterminalProductions
                    (equations, variablesToConstants, false);

            // binarize all productions, i.e. nonterminal productions of length at least 3 will be split up until
            // there are only productions of length <= 2 left
            systemBeingProcessed = NonCommutativePolynomial<LSR>::binarizeProductions(systemBeingProcessed, variablesToConstants);

            // add 1 to each equation
            for(auto &equation: systemBeingProcessed) {
                equation.second = equation.second + NonCommutativePolynomial<LSR>::one();
            }

            // add the monomials of degree 1; since we are in an idempotent semiring, we have a+a = a and so
            // we don't need to worry about the number of occurrences of each variable
            std::set<VarId> vars;
            NonCommutativePolynomial<LSR> monomialsOfDegreeOne;
            for(auto &equation: systemBeingProcessed) {
                vars = equation.second.get_variables();
                monomialsOfDegreeOne = NonCommutativePolynomial<LSR>::null();

                for(VarId var: vars) {
                    monomialsOfDegreeOne += NonCommutativePolynomial<LSR>(var);
                }

                equation.second = equation.second + monomialsOfDegreeOne;
            }

//            std::cout << "lossy QNF:" << std::endl;
//            for(auto &equation: systemBeingProcessed) {
//                std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
//            }

            return systemBeingProcessed;
    }

    /*
     * Used to decompose a system into its strongly connected components.
     *
     * The resulting vector will have the components bottom-up, i.e. if the graph is {A,B,C} with edges A->B, B->C, C->B,
     * then the vector will have component {B,C} at its first entry, the second entry will be {A}. Components between
     * which there is no directed path may appear in any order.
     */
    static std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> group_by_scc
        (const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations) {

        // create map of variables to [0..n]. this is used to enumerate important variables
        // in a clean way from 0 to n during graph construction
        std::unordered_map<VarId, int> var_key;

        // build the graph
        boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex<LSR>> graph(equations.size());

        for (const auto &eq : equations) {
            // if the variable is not yet in the map insert it together with the size of
            // the map. this way we get a unique identifier for the variable counting from 0 to n
            if (var_key.find(eq.first) == var_key.end()) {
                var_key.insert(std::make_pair(eq.first, var_key.size()));
            }

            int a = var_key.find(eq.first)->second; // variable key
            graph[a].var = eq.first; // store VarId to the vertex
            graph[a].rex = eq.second; // store the regular expression to the vertex

            auto variables = eq.second.get_variables(); // all variables of this rule;

            for (const auto &var : variables) {
                if (var_key.find(var) == var_key.end()) { // variable is not yet in the map
                    var_key.insert(var_key.begin(), std::pair<VarId,int>(var,var_key.size()));
                }

                int b = var_key.find(var)->second; // variable key
                boost::add_edge(a, b, graph);
            }
        }

        // calculate strongly connected components and store them in 'component'
        std::vector<int> component(boost::num_vertices(graph));
        boost::strong_components(graph,&component[0]);

        // group neccessary equations together
        int num_comp = *std::max_element(component.begin(), component.end()) + 1; // find the number of components
        std::vector<std::vector<std::pair<VarId,NonCommutativePolynomial<LSR>>>> grouped_equations(num_comp);

        // iterate over all vertices (0 to n)
        // collect the necessary variables + equations for every component
        for (std::size_t j = 0; j != component.size(); ++j) {
            grouped_equations[component[j]].push_back(std::pair<VarId,NonCommutativePolynomial<LSR>>(graph[j].var, graph[j].rex));
        }

        return grouped_equations;
    }

    /*
     * Cleans the polynomial system, i.e. removes variables that are unproductive for all the grammars we get
     * if exactly the symbols in "startSymbols" are considered as initial symbols of a grammar.
     *
     * In the context of fixpoints, if you only want to eliminate symbols that have a 0 fixpoint, hand this function an
     * empty queue "startSymbols".
     */
    static std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> cleanSystem
        (const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations, std::queue<VarId> &startSymbols) {

        if(startSymbols.empty()) {
            for(auto &equation: equations) {
                startSymbols.push(equation.first);
            }
        } else {
            std::queue<VarId> temp;

            while(!startSymbols.empty()) {
                auto var = startSymbols.front();
                startSymbols.pop();
                temp.push(var);
            }

            swap(temp,startSymbols);
        }

        // just hand the system to the cleaning procedure for grammars and use the set of all
        // variables for which there exist productions as the set of initial nonterminals
        return NonCommutativePolynomial<LSR>::cleanSystem(equations, startSymbols);
    }

    /*
     * Evaluates a polynomial system at a given vector.
     *
     * Not passing valuation as a reference is purposeful, since we need to copy it anyway as we don't want to change it.
     */
    static std::map<VarId, LSR> evaluateSystem (std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations,
        int &times, std::map<VarId, LSR> valuation) {

        std::map<VarId, LSR> tempValuation;
        VarId variable;

        // iterate the desired number of times
        for(int i = 0; i < times; i++) {

            // evaluate each polynomial and map the appropriate variable to the result
            for(auto &equation: equations) {
                tempValuation.insert(std::make_pair(equation.first, equation.second.eval(valuation)));
            }

            // prepare next iteration
            valuation.swap(tempValuation);
            tempValuation.clear();
        }

        return valuation;
    }
};

#endif
