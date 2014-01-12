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

template<typename SR>
class Evaluator;

class LossySemiring: public StarableSemiring<LossySemiring,
		Commutativity::NonCommutative, Idempotence::Idempotent> {
public:
	/* Default constructor creates zero element. */
	LossySemiring() {
		node_ = factory_.GetEmpty();
	}

	LossySemiring(const VarId var) {
		node_ = factory_.NewElement(var);
		// node_ = factory_.NewAddition(factory_.NewElement(var), factory_.GetEpsilon());
	}

	static LossySemiring null() {
		return LossySemiring { factory_.GetEmpty() };
	}

	static LossySemiring one() {
		return LossySemiring { factory_.GetEpsilon() };
	}

	LossySemiring star() const {
		return LossySemiring { factory_.NewStar(node_) };
	}

	LossySemiring operator+(const LossySemiring &x) {
		return LossySemiring { factory_.NewAddition(node_, x.node_) };
	}

	LossySemiring& operator+=(const LossySemiring &x) {
		node_ = factory_.NewAddition(node_, x.node_);
		return *this;
	}

	LossySemiring operator*(const LossySemiring &x) {
		return LossySemiring { factory_.NewMultiplication(node_, x.node_) };
	}

	LossySemiring& operator*=(const LossySemiring &x) {
		node_ = factory_.NewMultiplication(node_, x.node_);
		return *this;
	}

	bool operator==(const LossySemiring &x) const {
		return node_ == x.node_;
	}

	std::string string() const {
		return NodeToString(*node_);
	}

	std::string RawString() const {
		return NodeToRawString(*node_);
	}

	template<typename SR>
	SR Eval(const ValuationMap<SR> &valuation) const;

	template<typename SR>
	SR Eval(Evaluator<SR> &evaluator) const;

	void PrintDot(std::ostream &out) {
		factory_.PrintDot(out);
	}

	void PrintStats(std::ostream &out = std::cout) {
		factory_.PrintStats(out);
	}

	/*
	 * Solves a polynomial system over the lossy semiring. For derivation of
	 * the algorithm and proof of correctness, see "Esparza, Kiefer, Luttenberger:
	 * Derivation Tree Analysis for Accelerated Fixed-Point Computation".
	 */
	static ValuationMap<LossySemiring> solvePolynomialSystem(
			const std::vector<std::pair<VarId, NonCommutativePolynomial<LossySemiring>>>&equations) {

				ValuationMap<LossySemiring> solution;
				std::vector<std::pair<VarId, NonCommutativePolynomial<LossySemiring>>> workEquations = cleanSystem(equations);

				// build a zero vector
				std::map<VarId, LossySemiring> zeroSystem;
				for(auto &equation: workEquations) {
					zeroSystem.insert(std::pair<VarId, LossySemiring>(equation.first, LossySemiring::null()));
				}

				// build the vectors f(0) and f^n(0), where n is the number of variables in the system
				int times = 1;
				std::map<VarId, LossySemiring> f_0 = evaluateSystem(workEquations, times, zeroSystem);
				times = workEquations.size();
				std::map<VarId, LossySemiring> f_n_0 = evaluateSystem(workEquations, times, zeroSystem);

				// find the LossySemiring element in the "middle" of the expression
				LossySemiring middle = LossySemiring::null();
				for(auto &elem_mapping: f_0) {
					middle += elem_mapping.second;
				}

				// build the differential of the system
				std::map<VarId, NonCommutativePolynomial<LossySemiring>> differential;
				for(auto &equation: workEquations) {
					differential.insert(std::make_pair(equation.first, equation.second.differential_at(f_0)));
				}

				// sum all polynomials of the differential of the system; we don't need attribution
				// of each polynomial to the respective variable since we only care about the
				// leading and trailing coefficients of each monomial
				NonCommutativePolynomial<LossySemiring> differential_sum;
				for(auto &equation: differential) {
					differential_sum += equation.second;
				}
				std::cout << differential_sum << std::endl;
				// get the lefthand and righthand semiring element of the fixpoint
				LossySemiring lefthandSum = differential_sum.getSumOfLeadingFactors();
				LossySemiring righthandSum = differential_sum.getSumOfTrailingFactors();

				// fixpoint element
				LossySemiring fixpoint = lefthandSum.star() * middle * righthandSum.star();
				for(auto &variable_mapping:workEquations) {
					solution.insert(std::make_pair(variable_mapping.first, fixpoint));
				}

				//std::map<VarId, NonCommutativePolynomial<LossySemiring>>
				return solution;
			}

			void GC() {
				factory_.GC();
			}

private:
	LossySemiring(NodePtr n) : node_(n) {}

	NodePtr node_;
	static NodeFactory factory_;

	friend struct std::hash<LossySemiring>;

	/*
	 * Cleans the polynomial system, i.e. removes variables that are unproductive.
	 */
	static std::vector<std::pair<VarId, NonCommutativePolynomial<LossySemiring>>> cleanSystem
		(const std::vector<std::pair<VarId, NonCommutativePolynomial<LossySemiring>>> &equations) {

		std::queue<VarId> temp1, temp2; // to store the variables that still need checking
		std::map<VarId, bool> productiveVariables;
		std::map<VarId, NonCommutativePolynomial<LossySemiring>> polyMap;

		// build the data structures that will store the info about which variables still
		// need checking and which variables are known to be productive
		for(auto equation: equations) {
			temp1.push(equation.first);
			productiveVariables.insert(std::make_pair(equation.first, false));
			polyMap.insert(std::make_pair(equation.first, equation.second));
		}

		// keep checking the variables that are left until no further variable becomes known
		// to be productive
		bool update = true;
		VarId var;
		while(update) {
			update = false;

			while(!temp1.empty()) {
				var = temp1.front();
				temp1.pop();

				// if the variable was found to be productive, remember the update
				// and change the flag for the variable
				if(polyMap[var].isProductive(productiveVariables)) {
					productiveVariables[var] = true;
					update = true;
				} else { // if the variable isn't productive, put it in the queue
					 temp2.push(var);
				}
			}

			temp1.swap(temp2);
		}

		// build the clean system: only use productive variables, and remove unproductive monomials
		std::vector<std::pair<VarId, NonCommutativePolynomial<LossySemiring>>> cleanEquations;
		for(auto equation: equations) {
			if(productiveVariables[var]) {
				cleanEquations.push_back(std::make_pair(var, polyMap[var].removeUnproductiveMonomials(productiveVariables)));
			}
		}

		return cleanEquations;
	}

	/*
	 * Evaluates a polynomial system at a given vector.
	 */
	static std::map<VarId, LossySemiring> evaluateSystem
		(std::vector<std::pair<VarId, NonCommutativePolynomial<LossySemiring>>> &equations,
		int& times, std::map<VarId, LossySemiring> &initialValuation) {

		std::map<VarId, LossySemiring> valuation, tempValuation;
		valuation = initialValuation;
		VarId variable;

		// iterate the desired number of times
		for(int i = 0; i < times; i++) {

			// evaluate each polynomial and map the appropriate variable to the result
			for(auto &equation: equations) {
				tempValuation.insert(std::pair<VarId, LossySemiring>(equation.first, equation.second.eval(valuation)));
			}

			// prepare next iteration
			valuation.swap(tempValuation);
			tempValuation.clear();
		}

		return valuation;
	}

};

#endif
