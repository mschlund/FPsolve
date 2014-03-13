#ifndef LOSSY_SEMIRING_H_
#define LOSSY_SEMIRING_H_

#include <string>
#include <memory>
#include <unordered_map>
#include <queue>

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/free-structure.h"
#include "../polynomials/non_commutative_polynomial.h"

#include "semiring.h"

template<typename LSR>
class Evaluator;

template <typename LSR>
class LossySemiring: public StarableSemiring<LSR,
		Commutativity::NonCommutative, Idempotence::Idempotent> {
public:
    virtual ~LossySemiring(){};

    /*
     * Solves a polynomial system over the lossy semiring. For derivation of
     * the algorithm and proof of correctness, see "Esparza, Kiefer, Luttenberger:
     * Derivation Tree Analysis for Accelerated Fixed-Point Computation".
     */
    static ValuationMap<LSR> solvePolynomialSystem(
        const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>&equations) {

        ValuationMap<LSR> solution;
        std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> normalForm = quadraticNormalForm(equations);

        // build a zero vector
        std::map<VarId, LSR> zeroSystem;
        for(auto &equation: normalForm) {
            zeroSystem.insert(std::pair<VarId, LSR>(equation.first, LSR::null()));
        }

        // build the vectors f(0) and f^n(0), where n is the number of variables in the system
        int times = 1;
        std::map<VarId, LSR> f_0 = evaluateSystem(normalForm, times, zeroSystem);
        times = normalForm.size();
        std::map<VarId, LSR> f_n_0 = evaluateSystem(normalForm, times, zeroSystem);

        // find the LossySemiring element in the "middle" of the expression
        LSR middle = LSR::null();
        for(auto &elem_mapping: f_0) {
            middle += elem_mapping.second;
        }

        // build the differential of the system
        std::map<VarId, NonCommutativePolynomial<LSR>> differential;
        for(auto &equation: normalForm) {
            differential.insert(std::make_pair(equation.first, equation.second.differential_at(f_0)));
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

        for(auto &variable_mapping: normalForm) {
            solution.insert(std::make_pair(variable_mapping.first, fixpoint));
        }

        return solution;
    }

private:

    /*
     * Brings a system into quadratic normal form.
     */
    static std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> quadraticNormalForm
            (const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations) {
            std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> normalFormSystem;

            // clean the system
            normalFormSystem = cleanSystem(equations);

//            std::cout << "clean system:" << std::endl;
//            for(auto &equation: normalFormSystem) {
//                std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
//            }

            // bring the clean system into Chomsky Normal Form
            normalFormSystem = NonCommutativePolynomial<LSR>::chomskyNormalForm(normalFormSystem);

            std::cout << "CNF system:" << std::endl;
            for(auto &equation: normalFormSystem) {
                std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
            }


            // add 1 to each equation
            for(auto &equation: normalFormSystem) {
                equation.second = equation.second + NonCommutativePolynomial<LSR>::one();
            }

            std::cout << "CNF + 1:" << std::endl;
            for(auto &equation: normalFormSystem) {
                std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
            }

            // add the monomials of degree 1; since we are in a lossy semiring, we have 1+1 = 1 and so
            // we don't need to worry about the number of occurrences of each variable
            std::set<VarId> vars;
            NonCommutativePolynomial<LSR> monomialsOfDegreeOne;
            for(auto &equation: normalFormSystem) {
                vars = equation.second.get_variables();
                monomialsOfDegreeOne = NonCommutativePolynomial<LSR>::null();

                for(VarId var: vars) {
                    monomialsOfDegreeOne += NonCommutativePolynomial<LSR>(var);
                }

                equation.second = equation.second + monomialsOfDegreeOne;
            }

            std::cout << "QNF:" << std::endl;
            for(auto &equation: normalFormSystem) {
                std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
            }

            return normalFormSystem;
    }

    /*
     * Cleans the polynomial system, i.e. removes variables that are unproductive.
     */
    static std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> cleanSystem
        (const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations) {

        std::queue<VarId> worklist;

        for(auto &equation: equations) {
            worklist.push(equation.first);
        }

        // just hand the system to the cleaning procedure for grammars and use the set of all
        // variables for which there exist productions as the set of initial nonterminals
        return NonCommutativePolynomial<LSR>::cleanSystem(equations, worklist);
    }

    /*
     * Evaluates a polynomial system at a given vector.
     */
    static std::map<VarId, LSR> evaluateSystem
        (std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations,
        int &times, std::map<VarId, LSR> &initialValuation) {

        std::map<VarId, LSR> valuation, tempValuation;
        valuation = initialValuation;
        VarId variable;

        // iterate the desired number of times
        for(int i = 0; i < times; i++) {

            // evaluate each polynomial and map the appropriate variable to the result
            for(auto &equation: equations) {
                tempValuation.insert(std::pair<VarId, LSR>(equation.first, equation.second.eval(valuation)));
            }

            // prepare next iteration
            valuation.swap(tempValuation);
            tempValuation.clear();
        }

        return valuation;
    }
};

#endif
