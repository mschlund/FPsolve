#pragma once

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <list>
#include <map>
#include <queue>
#include <string>
#include <climits>

#include "../semirings/semiring.h"
#include "../semirings/free-semiring.h"

#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/var_degree_map.h"

#include "non_commutative_monomial.h"
#include "polynomial.h"

template<typename SR>
class NonCommutativePolynomial: public Semiring<NonCommutativePolynomial<SR>,
		Commutativity::NonCommutative, SR::GetIdempotence()> {
private:
	std::map<NonCommutativeMonomial<SR>, std::uint_fast16_t> monomials_;

	static void InsertMonomial(
			std::map<NonCommutativeMonomial<SR>, std::uint_fast16_t> &monomials,
			NonCommutativeMonomial<SR> monomial, std::uint_fast16_t coeff) {
	    //std::cout << "NPC<SR>::IM2" << std::endl;
		auto iter = monomials.find(monomial);
		if (iter == monomials.end()) {
			monomials.insert( { monomial, 1 });
		} else {
			auto tmp = *iter;
			monomials.erase(iter);
			monomials.insert( { tmp.first, tmp.second + coeff });
		}
	}

	static void InsertMonomial(
			std::map<NonCommutativeMonomial<SR>, std::uint_fast16_t> &monomials,
			std::vector<std::pair<elemType, int>> &idx,
			std::vector<VarId> &variables, std::vector<SR> &srs) {
	  //  std::cout << "NPC<SR>::IM1" << std::endl;
		auto tmp_monomial = NonCommutativeMonomial<SR>(idx, variables, srs);
		InsertMonomial(monomials, tmp_monomial, 1);
	}
	;


    /*
     * Binarizes this polynomial, as one would when calculating a CNF of a grammar. Make sure all monomials
     * either have no terminals (i.e. semiring elements) or no nonterminals (i.e. variables) as factors.
     *
     * If a monomial consists of exactly one variable, it will not be changed (chain productions are not eliminated).
     *
     * Make sure there are no productions with both terminals and nonterminals as factors; those will be eliminated if
     * they have at least 3 variables as factors. Call eliminateTerminalsInNonterminalProductions(..) first to avoid this.
     *
     * The variables and productions that are introduced during that process will be stored in binarizationVariables
     * and binarizationVariablesEquations, respectively.
     */
    NonCommutativePolynomial<SR> binarize(std::map<std::string, VarId> &binarizationVariables,
            std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &binarizationVariablesEquations,
            std::map<VarId, SR> &variablesToConstants) const {
        NonCommutativePolynomial<SR> binarizedPoly = null();

        // delegate to the monomials
        for(auto &monomial: monomials_) {
            if(monomial.first.get_degree() > 2) { // only binarize productions that have at least 3 variables in them
                binarizedPoly += monomial.first.binarize(binarizationVariables, binarizationVariablesEquations, variablesToConstants);
            } else {
                InsertMonomial(binarizedPoly.monomials_, monomial.first, monomial.second);
            }
        }

        return binarizedPoly;
    }

public:
	NonCommutativePolynomial() = default;

	NonCommutativePolynomial(const NonCommutativePolynomial &p) = default;
	NonCommutativePolynomial(NonCommutativePolynomial &&p) = default;

	/* Create a 'constant' polynomial. */
	NonCommutativePolynomial(const SR &elem) {
//        std::cout << "NCP(cSR)" << std::endl;
//        std::cout << "elem:\t" + elem.string() << std::endl;
		std::vector<std::pair<elemType, int>> idx = { {SemiringType, 0}};
		std::vector<VarId> variables = {};
		std::vector<SR> srs = {elem};
		if(!(elem == SR::null()))
		{
			InsertMonomial(monomials_, idx, variables, srs);
		} // else this is an empty polynomial
	}

	NonCommutativePolynomial(SR &elem) {
	    //std::cout << "NCP(SR)" << std::endl;
		std::vector<std::pair<elemType, int>> idx = { {SemiringType, 0}};
		std::vector<VarId> variables = {};
		std::vector<SR> srs = {elem};
		if(!(elem == SR::null()))
		{
			InsertMonomial(monomials_, idx, variables, srs);
		} // else this is an empty polynomial
	}

	/* Create a polynomial which consists only of one variable. */
	NonCommutativePolynomial(const VarId var) {
		std::vector<std::pair<elemType, int>> idx = { {Variable, 0}};
		std::vector<VarId> variables = {var};
		std::vector<SR> srs = {};
		InsertMonomial(monomials_, idx, variables, srs);
	}

	/* initialize with some monomials */
	NonCommutativePolynomial(std::initializer_list<NonCommutativeMonomial<SR>> monomials) {
		for(auto monomial : monomials) {
			InsertMonomial(monomials_, monomial, 1);
		}
	}

	NonCommutativePolynomial<SR>& operator=(const NonCommutativePolynomial<SR> &p) = default;
	NonCommutativePolynomial<SR>& operator=(NonCommutativePolynomial<SR> &&p) = default;

	NonCommutativePolynomial<SR>& operator+=(const NonCommutativePolynomial<SR> &polynomial) {
		if( *this == null())
		{
			*this = polynomial;
			return *this;
		}
		for (const auto &monomial : polynomial.monomials_) {
			// InsertMonomial handles monomials which are already in 'this' polynomial
			InsertMonomial(monomials_, monomial.first, monomial.second);
		}

		return *this;
	}

	NonCommutativePolynomial<SR>& operator+=(const NonCommutativeMonomial<SR> &monomial) {
		InsertMonomial(monomials_, monomial, 1);
		return *this;
	}

	NonCommutativePolynomial<SR>& operator+(const NonCommutativeMonomial<SR> &monomial) const {
		auto result = *this;
		result += monomial;
		return result;
	}

	NonCommutativePolynomial<SR>& operator*=(const NonCommutativePolynomial<SR> &rhs) {
		if (*this == one())
		{
			*this = rhs;
			return *this;
		} else if (rhs == one())
		{
			return *this;
		}
		if (*this == null())
		{
			return *this;
		} else if (rhs == null())
		{
			*this = rhs;
			return *this;
		}
		if (monomials_.empty()) {
			return *this;
		} else if (rhs.monomials_.empty()) {
			monomials_.clear();
			return *this;
		}

		std::map<NonCommutativeMonomial<SR>, std::uint_fast16_t> tmp_monomials;

		// iterate over both lhs und rhs
		for (const auto &lhs_monomial : monomials_) {
			for (const auto &rhs_monomial : rhs.monomials_) {
				auto tmp_monomial = lhs_monomial.first * rhs_monomial.first;
				auto tmp_coeff = lhs_monomial.second * rhs_monomial.second;

				InsertMonomial(tmp_monomials, tmp_monomial, tmp_coeff);
			}
		}
		monomials_ = std::move(tmp_monomials);
		return *this;
	}

	// multiplying a polynomial with a variable
	NonCommutativePolynomial<SR> operator*(const VarId &var) const {
		NonCommutativePolynomial result_polynomial;
		for (auto &monomial : monomials_) {
			InsertMonomial(result_polynomial.monomials_, monomial.first * var, monomial.second);
		}
		return result_polynomial;
	}

	friend NonCommutativePolynomial<SR> operator*(const SR &elem,
			const NonCommutativePolynomial<SR> &polynomial) {
		NonCommutativePolynomial result_polynomial;
		for (auto &monomial : polynomial.monomials_) {
			result_polynomial.monomials_ += monomial * elem;
		}
		return result_polynomial;
	}

	bool operator==(const NonCommutativePolynomial<SR> &polynomial) const {
		return monomials_ == polynomial.monomials_;
	}

	/* convert this non-commutative-polynomial to a commutative one */
	Polynomial<SR> make_commutative() const {
		Polynomial<SR> result;

		for(auto const &monomial : monomials_)
		{
			result *= monomial.make_commutative();
		}

		return result;
	}

	/* calculate the delta for this polynomial with the given data at iteration d
	 * de2 is data from newton iterand for d-2
	 * dl1 is data from newton update for d-1 */
	SR calculate_delta(
			const std::map<VarId, SR> &de2,
			const std::map<VarId, SR> &dl1) const {
		SR result = SR::null();
		for(auto const &monomial : monomials_) {
			result += monomial.first.calculate_delta(de2, dl1) * monomial.second;
		}

		return result;
	}

	/* calculates the derivative of this polynomial
	 * return a schema and the used mapping from X to X[d-1]*/
	std::pair<NonCommutativePolynomial<SR>,std::map<VarId,VarId>> derivative() const {
		NonCommutativePolynomial<SR> result;

		/* get a mapping for all variables X s.t. X -> X[d-1]
		 * X[d-1] has to be a fresh variable*/
		auto vars = get_variables();
		std::map<VarId, VarId> mapping;
		for(auto const &var : vars) {
			auto tmp = Var::GetVarId();
			mapping.insert( {var, tmp});
		}

		/* use the mapping and sum up all derivatives of all monomials */
		for(auto const &monomial : monomials_) {
			result += monomial.first.derivative(mapping) * monomial.second;
		}

		return {result, mapping};
	}

	NonCommutativePolynomial<SR> differential_at(const std::map<VarId,SR> valuation) const {
		NonCommutativePolynomial<SR> result;
		for(auto const &monomial : monomials_) {
			result += monomial.first.differential_at(valuation) * monomial.second;
		}
		return result;
	}

	SR eval(const std::map<VarId, SR> &values) const {
		SR result = SR::null();
		for (const auto &monomial : monomials_) {
			result = result + (monomial.first.eval(values) * monomial.second);
		}


//        std::cout << "poly result of eval:\t\t" + result.string() << std::endl;
		return result;
	}

	/* Evaluate as much as possible (depending on the provided values) and
	 * return the remaining (i.e., unevaluated) monomial. */
	NonCommutativePolynomial<SR> partial_eval(const std::map<VarId, SR> &values) const {
		NonCommutativePolynomial<SR> result = NonCommutativePolynomial<SR>::null();
		for(const auto &monomial : monomials_) {
			result.monomials_.insert(std::make_pair(monomial.first.partial_eval(values), monomial.second));
		}
		return result;
	}

	/* Variable substitution. */
	NonCommutativePolynomial<SR> subst(const std::map<VarId, VarId> &mapping) const {
		NonCommutativePolynomial result_polynomial;

		for (const auto &monomial : monomials_) {
			result_polynomial.monomials_.insert( {monomial.first.subst(mapping), monomial.second});
		}
		return result_polynomial;
	}

	static Matrix<SR> eval(const Matrix< NonCommutativePolynomial<SR> > &poly_matrix,
			const std::map<VarId, SR> &values) {
		const std::vector<NonCommutativePolynomial<SR>> &tmp_polynomials = poly_matrix.getElements();
		std::vector<SR> result;
		for (const auto &polynomial : tmp_polynomials) {
			result.emplace_back(polynomial.eval(values));
		}
		return Matrix<SR> {poly_matrix.getRows(), std::move(result)};
	}

	static Matrix<NonCommutativePolynomial<SR> > eval(Matrix<NonCommutativePolynomial<SR> > poly_matrix,
			std::map<VarId,NonCommutativePolynomial<SR> > values) {
		const std::vector<NonCommutativePolynomial<SR>> &tmp_polynomials = poly_matrix.getElements();
		std::vector<NonCommutativePolynomial<SR>> result;
		for (const auto &polynomial : tmp_polynomials) {
			result.emplace_back(polynomial.eval(values));
		}
		return Matrix<NonCommutativePolynomial<SR>> {poly_matrix.getRows(), result};
	}

	/* Convert this polynomial to an element of the free semiring.  Note that
	 * the provided valuation might be modified with new elements. */
	FreeSemiring make_free(std::unordered_map<SR, VarId, SR> *valuation) const {
		assert(valuation);

		auto result = FreeSemiring::null();
		// convert this polynomial by adding all converted monomials
		for (const auto &monomial : monomials_) {
			result += monomial.first.make_free(valuation) * monomial.second;
		}
		return result;
	}

	/* Same as make_free but for matrix form. */
	static Matrix<FreeSemiring> make_free(
			const Matrix< NonCommutativePolynomial<SR> > &poly_matrix,
			std::unordered_map<SR, VarId, SR> *valuation) {

		assert(valuation);

		const std::vector< NonCommutativePolynomial<SR> > &tmp_polynomials = poly_matrix.getElements();
		std::vector<FreeSemiring> result;

		for (const auto &polynomial : tmp_polynomials) {
			result.emplace_back(polynomial.make_free(valuation));
		}
		return Matrix<FreeSemiring> {poly_matrix.getRows(), std::move(result)};
	}

	Degree get_degree() {
		Degree degree = 0;
		for (auto &monomial : monomials_) {
			degree = std::max(degree, monomial.first.get_degree());
		}
		return degree;
	}

	// FIXME: get rid of this...
	// if you use this, make sure, you cache it...
	std::set<VarId> get_variables() const {
		std::set<VarId> vars;
		for(auto const &monomial : monomials_)
		{
			auto tmp = monomial.first.get_variables();
			vars.insert(tmp.begin(), tmp.end());
		}
		return vars;
	}

	/*
	 * Returns the sum of the leading constant factors of all monomials in this polynomial.
	 */
	SR getSumOfLeadingFactors() {
		SR sum = SR::null();

		for(auto &monomial: monomials_) {
			sum += monomial.first.getLeadingSR();
		}

		return sum;
	}

	/*
	 * Returns the sum of the leading constant factors of all monomials in this polynomial.
	 */
	SR getSumOfTrailingFactors() {
		SR sum = SR::null();

		for(auto &monomial: monomials_) {
			sum += monomial.first.getTrailingSR();
		}

		return sum;
	}

	/*
	 * Finds all the constants appearing in nonterminal monomials of this polynomial and appends each new one to the vector.
	 *
	 * Only checks linear monomials if checkLinearTerms is true.
	 */
	void findConstants(std::set<SR> &constants, bool checkLinearTerms) const {

		// delegate to the monomials
		for(auto &monomial: monomials_) {
			monomial.first.findConstants(constants, checkLinearTerms);
		}
	}

	/*
	 * Replaces all constants in the nonterminal monomials of the polynomial with their respective VarId mapping.
	 * Used in transformation to Chomsky Normal Form.
	 *
	 * Only replaces constants in linear monomials if replaceInLinearMonomials is true.
	 */
	NonCommutativePolynomial<SR> replaceConstants(std::map<SR, VarId> &constantsToVariables,
	        bool replaceInLinearMonomials) const {
		NonCommutativePolynomial<SR> temp = null();

		// delegate to the monomials
		for(auto &monomial: monomials_) {

		    // we only replace constants in nonterminal monomials, all terminal monomials are added to the polynomial as is
		    // only replace constants in linear monomials if replaceInLinearMonomials is true
		    if((monomial.first.get_degree() != 0)
		            && (replaceInLinearMonomials || monomial.first.get_degree() != 1)) {
	            temp += monomial.first.replaceConstants(constantsToVariables);
		    } else {
                InsertMonomial(temp.monomials_, monomial.first, monomial.second);
		    }
		}

		return temp;
	}

	/*
	 * Removes all epsilon monomials from a polynmoial.
	 * An epsilon monomial is a monomial consisting only of epsilon.
	 */
	NonCommutativePolynomial<SR> removeEpsilonMonomials() const {
	    NonCommutativePolynomial<SR> temp = null();

	    for(auto monomial: monomials_) {

	        // if the monomial has at least one variable or it doesn't have a variable but isn't
	        // equal to epsilon, then it's not an epsilon monomial
	        if(!monomial.first.isEpsilonMonomial()) {
	            temp += monomial.first * monomial.second;
	        }
	    }

	    return temp;
	}

	/*
	 * Removes all epsilon productions from the system.
	 */
	static void removeEpsilonProductions(std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &equations) {
        for(auto equation: equations) {
            equation.second = equation.second.removeEpsilonMonomials();
        }
	}

    /*
     * Cleans the polynomial system, i.e. removes variables that are unproductive or unreachable
     * from the set of nonterminals given in the worklist; the worklist needs to contain the
     * variables we want to start derivations from, i.e. the set of initial nonterminals.
     */
    static std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> cleanSystem
        (const std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &equations,
                std::queue<VarId> &worklist) {

//        std::cout << "system before cleaning:" << std::endl;
//        for(auto &equation: equations) {
//            std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
//        }

        std::queue<VarId> temp; // to store the variables that still need checking
        std::map<VarId, int> encounteredVariables; // to remember whether a variable was encountered and check for
                                                   // for being productive yet;
                                                   // 0 means "never encountered",
                                                   // 1 means "encountered but not checked",
                                                   // 2 means "checked at least once"
        std::map<VarId, bool> productiveVariables; // to map variables to "is this variable known to be productive?"
        std::map<VarId, NonCommutativePolynomial<SR>> productions; // to map variables to their productions

        VarId var;
        std::set<VarId> vars;
        // build the data structures that will store the info about which variables still
        // need checking and which variables are known to be productive
        for(auto &equation: equations) {
            vars = equation.second.get_variables();

            // just to make sure we catch all variables in case someone inputs a grammar with a terminal
            // that only appears in a lefthand side of a production
            for(auto &prodVar: vars) {
                encounteredVariables.insert(std::make_pair(prodVar, 0));
                productiveVariables.insert(std::make_pair(prodVar, false));
            }

            encounteredVariables.insert(std::make_pair(equation.first, 0));
            productiveVariables.insert(std::make_pair(equation.first, false));
            productions.insert(std::make_pair(equation.first, equation.second));

            // prepare the info for the initial variables
            while(!worklist.empty()) {
                var = worklist.front();
                worklist.pop();
                encounteredVariables[var] = 1;
                temp.push(var);
            }
            worklist.swap(temp);
        }

        // keep checking the variables that are not yet known to be productive
        // until no further variable becomes known to be productive; we do this
        // using BFS so unreachable variables will be removed as well since they
        // can never be found to be productive this way
        bool update;
        do {
            update = false;

            // terminates once there is nothing more to do; specifically, this terminates at the latest once
            // all reachable variables are known to be productive because then they won't be put in the
            // worklist again; it terminates at the earliest once all reachable variables have been checked
            // without finding a productive one since the worklist is being filled with all reachable variables
            // before the end of the first iteration
            while(!worklist.empty()) {
                var = worklist.front();
                worklist.pop();

                // if we have not yet checked this variable, remember that we now have and put all
                // the unencountered variables in its productions into the worklist; this way, we only
                // put reachable variables into the worklist once
                if(encounteredVariables[var] == 1) {
                    encounteredVariables[var] = 2;

                    for(auto &monomial: productions[var].monomials_) {
                        for(auto &reachableVar: monomial.first.get_variables()) {

                            // only put variables into the worklist that have not been seen before
                            if(encounteredVariables[reachableVar] == 0) {
                                worklist.push(reachableVar);
                                encounteredVariables[reachableVar] = 1;
                            }
                        }
                    }
                }

                // if the variable was found to be productive, remember the update and change the
                // flag for the variable; we know a variable is productive once one of the productions
                // associated with it becomes known to be productive; since there is no need to recheck a
                // productive variable, don't put it in the worklist again
                if(productions[var].isProductive(productiveVariables)) {
                    productiveVariables[var] = true;
                    update = true;
                } else { // if the variable isn't known to be productive yet, put it in the worklist
                         // for the next iteration
                     temp.push(var);
                }
            }

            // don't do work we don't need
            if(update) {
                worklist.swap(temp);
            }
        } while(update);

        // build the clean system: only use variables that are both productive and reachable;
        // remove unproductive monomials
        std::cout << "clean system:" << std::endl;
        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> cleanEquations;
        for(auto &equation: equations) {
            if(productiveVariables[equation.first]) {

                // at this point we know that the polynomial we are associating with the
                // variable is not the null polynomial since at least one of its monomials was
                // found to be productive
                cleanEquations.push_back(std::make_pair(equation.first,
                        equation.second.removeUnproductiveMonomials(productiveVariables)));
                std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
            }
        }

        return cleanEquations;
    }

    /*
     * Eliminates chain productions while creating the Chomsky Normal Form of a system. Puts the cleaned
     * system into "equationsWithoutChainProductions".
     */
    static void eliminateChainProductions (const std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &equations,
            std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &equationsWithoutChainProductions,
            std::map<VarId, NonCommutativePolynomial<SR>> &polyMap) {

        // used to mark the variables for which chain productions have already been eliminated
        std::set<VarId> hadChainsEliminated;
        std::map<VarId, NonCommutativePolynomial<SR>> eliminatedChainsPolyMap;

        for(auto &equation: equations) {
            NonCommutativePolynomial<SR> chainfreeProduction = NonCommutativePolynomial<SR>::null();
            std::queue<VarId> worklist;
            worklist.push(equation.first);
            std::set<VarId> visitedVars;
            visitedVars.insert(equation.first);

            while(!worklist.empty()) {
                VarId currVar = worklist.front();
                worklist.pop();

                // check if we already know a chainfree polynomial for the current variable
                // if no, check it for chain productions; if yes, just add that polynomial
                if(hadChainsEliminated.find(currVar) == hadChainsEliminated.end()) {

                    // check all productions the current variable produces for being chain productions
                    for(auto &monomial: polyMap.at(currVar).monomials_) {

                        // if the current production is a chain production, insert it into the work list
                        // if we have not already seen it; otherwise, add it to the cleaned productions
                        if(monomial.first.isChainProduction()) {
                            VarId newVar = *(monomial.first.get_variables().begin());

                            if(visitedVars.find(newVar) == visitedVars.end()) {
                                worklist.push(newVar);
                                visitedVars.insert(newVar);
                            }
                        } else {
                            InsertMonomial(chainfreeProduction.monomials_, monomial.first, monomial.second);
                        }
                    }
                } else {
                    assert(eliminatedChainsPolyMap.find(currVar) != eliminatedChainsPolyMap.end());
                    chainfreeProduction += eliminatedChainsPolyMap[currVar];
                }
            }

            // remember that we eliminated the chain productions for this variable, also remember the "clean" production
            hadChainsEliminated.insert(equation.first);
            eliminatedChainsPolyMap.insert(std::make_pair(equation.first, chainfreeProduction));
            equationsWithoutChainProductions.push_back(std::make_pair(equation.first, chainfreeProduction));
        }
    }

    /*
     * Equations can contain productions that are sentential forms, and variablefiedEquations will afterwards hold
     * productions that are either terminal or only contain nonterminals (i.e. the monomials in each polynomial
     * are either semiring elements or only contain variables).
     *
     * Linear terms (i.e. such with exactly one variable in them) will have their terminals replaced by variables
     * iff eliminateInLinearTerms is true. This way, we can produce a quadratic normal form with needlessly blowing
     * up linear terms, for example.
     */
    static std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> eliminateTerminalsInNonterminalProductions
                (const std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &equations,
                       std::map<VarId, SR> &variablesToConstants, bool eliminateInLinearTerms) {

        std::set<SR> constants; // will hold all terminals appearing in nonterminal productions

        // find all constants in nonterminal productions in the system
        for(auto &equation: equations) {
            equation.second.findConstants(constants, eliminateInLinearTerms);
        }

        // return value
        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> variablefiedEquations;

        // introduce a new variable for each constant found, add the respective equation to the system
        VarId var;
        std::map<SR, VarId> constantsToVariables; // will hold a mapping from terminals to variables that produce them
        for(auto &constant: constants) {
            var = Var::GetVarId();
            variablefiedEquations.push_back(std::make_pair(var, NonCommutativePolynomial<SR>(constant)));
            constantsToVariables.insert(std::make_pair(constant, var));
            variablesToConstants.insert(std::make_pair(var,constant));
            std::cout << "new production during constant elimination:\t" << Var::GetVar(var).string() << " -> " << constant.string() << std::endl;
        }

        // replace the constants in each nonterminal production with the new constant variables
        NonCommutativePolynomial<SR> allVariablesPoly;
        for(auto &equation: equations) {
            allVariablesPoly = equation.second.replaceConstants(constantsToVariables, eliminateInLinearTerms);
            variablefiedEquations.push_back(std::make_pair(equation.first, allVariablesPoly));
        }

        return variablefiedEquations;
    }

    /*
     * Binarizes nonterminal productions the way it is done when constructing a CNF of a grammar.
     */
    static std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> binarizeProductions
                (const std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &equations,
                       std::map<VarId, SR> &variablesToConstants) {

        // stores mappings between suffixes of monomials and the variables that are introduced
        // during binarization to produce those suffixes
        std::map<std::string, VarId> binarizationVariables;

        // stores the productions (i.e. equations) that are introduced while producing the
        // Chomsky Normal Form
        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> binarizationVariablesEquations;

        // to store the result
        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> binarizedEquations;

        // determine the necessary new variables to give all non-terminal productions that need binarizing
        // (i.e. monomials of degree > 2) the form "X = YZ"
        NonCommutativePolynomial<SR> binarizedPoly;
        for(auto &equation: equations) {
            binarizedPoly =
                    equation.second.binarize(binarizationVariables, binarizationVariablesEquations, variablesToConstants);
            binarizedEquations.push_back(std::make_pair(equation.first, binarizedPoly));
        }

        // finally, add the productions to the system that were introduced during the last step
        for(auto &equation: binarizationVariablesEquations) {
            binarizedEquations.push_back(std::make_pair(equation.first, equation.second));
        }

        return binarizedEquations;
    }

    /*
     * Brings a polynomial system into CNF. The resulting system uses a superset of the nonterminals
     * used before the normalization, i.e.: if the system was a CFG beforehand with a designated start
     * symbol, then that start symbol can still be used.
     *
     * WARNING: the result will have a higher number of appearances of some monomials in the productions
     * (e.g. you might get 3*(SabcXVad) instead of 2*(SabcXVad) that). If you want to fix this,
     * have a look at eliminateChainProductions(..) and remove the memoization of "clean" productions.
     */
//    static std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> chomskyNormalForm
//                (const std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &equations) {
//
//        // only the completely done stuff goes in here
//        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> chomskyNormalFormEquations;
//
//        std::map<SR, VarId> constantsToVariables;
//        std::map<VarId, SR> variablesToConstants;
//        std::set<SR> constants;
//
//        // map for quick access to productions
//        std::map<VarId, NonCommutativePolynomial<SR>> polyMap;
//        for(auto &equation: equations) {
//            polyMap.insert(std::make_pair(equation.first, equation.second));
//        }
//
//        // first point of order: chain productions
//        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> noChains;
//        eliminateChainProductions(equations, noChains, polyMap);
//
////        std::cout << "chains eliminated:" << std::endl;
////        for(auto &equation: noChains) {
////            std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
////        }
//
//        // find all constants in nonterminal productions in the system
//        for(auto &equation: noChains) {
//            equation.second.findConstantsInNonterminalMonomials(constants);
//        }
//
//        // introduce a new variable for each constant found, add the respective equation to the system
//        VarId var;
//        for(auto &constant: constants) {
//            var = Var::GetVarId();
//            noChains.push_back(std::make_pair(var, NonCommutativePolynomial<SR>(constant)));
//            constantsToVariables.insert(std::make_pair(constant, var));
//            variablesToConstants.insert(std::make_pair(var,constant));
//        }
//
//        // replace the constants in each nonterminal production with the new constant variables
//        NonCommutativePolynomial<SR> allVariablesPoly;
//        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> variablefiedEquations;
//        for(auto &equation: noChains) {
//            allVariablesPoly = equation.second.replaceConstants(constantsToVariables);
//            variablefiedEquations.push_back(std::make_pair(equation.first, allVariablesPoly));
//        }
//
////        std::cout << "done replacing constants, system so far:" << std::endl;
////        for(auto &equation: variablefiedEquations) {
////            std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
////        }
//
//        // stores mappings between suffixes of monomials and the variables
//        // that are introduced to produce those suffixes in the process of generating
//        // the Chomsky Normal Form
//        std::map<std::string, VarId> chomskyVariables;
//
//        // stores the productions (i.e. equations) that are introduced while producing the
//        // Chomsky Normal Form
//        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> chomskyVariableEquations;
//
//        // determine the necessary new variables to give all non-terminal productions
//        // (i.e. monomials of degree > 2) the form "X = YZ"
//        NonCommutativePolynomial<SR> chomskyPoly;
//        for(auto &equation: variablefiedEquations) {
//            chomskyPoly =
//                    equation.second.binarize(chomskyVariables, chomskyVariableEquations, variablesToConstants);
//            chomskyNormalFormEquations.push_back(std::make_pair(equation.first, chomskyPoly));
//        }
//
//        // finally, add the productions to the system that were introduced during the last step
//        for(auto &equation: chomskyVariableEquations) {
//            chomskyNormalFormEquations.push_back(std::make_pair(equation.first, equation.second));
//        }
//
//        // we are done
//        return chomskyNormalFormEquations;
//    }



//    /*
//     * Cleans the polynomial system, i.e. removes variables that are unproductive.
//     */
//    static std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> cleanSystem
//        (const std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &equations) {
//
//        std::queue<VarId> worklist, temp; // to store the variables that still need checking
//        std::map<VarId, bool> productiveVariables;
//        std::map<VarId, NonCommutativePolynomial<SR>> polyMap;
//
//        // build the data structures that will store the info about which variables still
//        // need checking and which variables are known to be productive
//        for(auto equation: equations) {
//            worklist.push(equation.first);
//            productiveVariables.insert(std::make_pair(equation.first, false));
//            polyMap.insert(std::make_pair(equation.first, equation.second));
//        }
//
//        // keep checking the variables that are not yet known to be productive
//        // until no further variable becomes known to be productive
//        bool update = true;
//        VarId var;
//        while(update) {
//            update = false;
//
//            while(!worklist.empty()) {
//                var = worklist.front();
//                worklist.pop();
//
//                // if the variable was found to be productive, remember the update
//                // and change the flag for the variable
//                if(polyMap[var].isProductive(productiveVariables)) {
//                    productiveVariables[var] = true;
//                    update = true;
//                } else { // if the variable isn't productive, put it in the queue for the next iteration
//                     temp.push(var);
//                }
//            }
//
//            worklist.swap(temp);
//        }
//
//        // build the clean system: only use productive variables and remove unproductive monomials
//        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> cleanEquations;
//        for(auto equation: equations) {
//            if(productiveVariables[equation.first]) {
//                cleanEquations.push_back(std::make_pair(equation.first,
//                        equation.second.removeUnproductiveMonomials(productiveVariables)));
//            }
//        }
//
//        return cleanEquations;
//    }

	/*
	 * Checks whether this polynomial is productive, depending on the set of variables
	 * that are already known to be productive.
	 */
	bool isProductive(const std::map<VarId, bool> &productiveVariables) const {

		// check if any monomial in this polynomial is productive; that will be enough
		// for the polynomial to be productive
		for(auto &monomial: monomials_) {
			if(monomial.first.isProductive(productiveVariables)) {
				return true;
			}
		}

		// if there is no productive monomial, then the polynomial isn't known to be productive so far
		return false;
	}

	/*
	 * Used during the generation of the intersection between a CFG and a FiniteAutomaton.
	 * Generates all productions for a triple in states x states x nonterminals, represented
	 * by their indices in the vectors "states" and "oldGrammar",
	 * where "oldGrammar[nonterminalIndex].first" is the nonterminal we care about.
	 *
	 * See FiniteAutomaton::intersectionWithCFG(..) for details.
	 *
	 * WARNING: if you use an SR where the string representations of elements are NOT strings over
	 * the English alphabet (e.g. they may contain 1s to represent epsilons, or they may use a
	 * different character set), then this function will possibly not do what you want it to do.
	 * The problem we run into there is that we have to process all productions of the form
	 * A -> X_1 X_2 X_3 ... X_m where the (X_i)s are all elements of either the terminal or the
	 * nonterminal alphabet, so having any symbols in the grammar that don't belong to either
	 */
	NonCommutativePolynomial<SR> intersectionPolynomial(std::vector<unsigned long> &states,
	        std::map<unsigned long, std::map<unsigned char, std::forward_list<unsigned long>>> &transitionTable,
	        unsigned long &startState,
	        unsigned long &targetState,
	        std::map<unsigned long, unsigned long> &statesToIndices,
	        std::map<VarId, unsigned long> &oldVariablesToIndices,
	        std::vector<std::vector<std::vector<VarId>>> &newVariables) const {

	    NonCommutativePolynomial<SR> result = NonCommutativePolynomial<SR>::null();

	    // delegate to the monomials
	    for(auto &monomial: monomials_) {

	        // if the nonterminal can produce epsilon, then the new grammar can produce epsilon only without
	        // the FA changing state (since the FA does not have epsilon transitions)
	        if(monomial.first.isEpsilonMonomial()) {
	            if(startState == targetState) {
	                result += NonCommutativePolynomial<SR>::one();
	            }
	        } else { // if this is a non-epsilon production, we calculate the new productions that derive from it
	            result += monomial.first.intersectionPolynomial
	                    (states, transitionTable, startState, targetState,
	                            statesToIndices, oldVariablesToIndices, newVariables);
	        }
	    }

	    return result;
	}

	/*
	 * Finds a shortest word in the language described by the grammar starting from "startSymbol" and having
	 * productions "productions". Will be null if the language is empty. Undefined behavior if startSymbol is
	 * not productive in the grammar defined by productions - "productions" must define a clean grammar starting
	 * from "startSymbol".
	 */
	static SR shortestWord(std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &productions, VarId &startSymbol) {

	    // check that the grammar contains the start symbol
        bool startSymbolIsProductive = false;
        for(auto &equation: productions) {
            if(equation.first == startSymbol) {
                startSymbolIsProductive = true;
            }
        }

        if(!startSymbolIsProductive) {
            return SR::null();
        }

	    std::map<VarId, unsigned long> lengthOfShortestWords;
	    std::map<VarId, NonCommutativeMonomial<SR>> productionsForShortestWords;

	    for(auto &production: productions) {
	        lengthOfShortestWords.insert(std::make_pair(production.first, ULONG_MAX));
	    }

	    // while some shorter derivable word was found, check whether we can now find some other shorter
	    // derivable word; terminates once we know a shortest derivable word for every nonterminal
	    bool update;
	    do {
	        update = false;

	        for(auto &equation: productions) {
	            for(auto &monomial: equation.second.monomials_) {
	                update |= monomial.first.findLengthOfDerivableStrings(lengthOfShortestWords, productionsForShortestWords, equation.first);
	            }
	        }
	    } while(update);

	    assert(productionsForShortestWords.find(startSymbol) != productionsForShortestWords.end());
	    return productionsForShortestWords[startSymbol].shortestDerivableTerminal(productionsForShortestWords);
	}

	/*
	 * Returns a cleaned version of this polynomial, i.e. a version that
	 * had all monomials with unproductive variables eliminated.
	 */
	NonCommutativePolynomial<SR> removeUnproductiveMonomials(const std::map<VarId, bool> &productiveVariables) const {
		NonCommutativePolynomial<SR> cleanPoly = null();

		// add all monomials to the clean polynomial that only contain productive variables
		for(auto &monomial: monomials_) {
			if(monomial.first.isProductive(productiveVariables)) {
				InsertMonomial(cleanPoly.monomials_, monomial.first, monomial.second);
			}
		}

		return cleanPoly;
	}

	static NonCommutativePolynomial<SR> null() {
		return NonCommutativePolynomial<SR> {};
	}

	static NonCommutativePolynomial<SR> one() {
	    //std::cout << "NCP<SR>::one()" << std::endl;
		return NonCommutativePolynomial<SR> {SR::one()};
	}

	static bool is_idempotent;
	static bool is_commutative;

	std::string string() const {
		// TODO: implement this
		std::stringstream ss;
		for(auto monomial = monomials_.begin(); monomial != monomials_.end(); monomial++) {
			if(monomial != monomials_.begin())
			ss << " + ";
//    	ss << "[degree: " << monomial->first.get_degree() << "]";
			/*ss << monomial->second << " * ";*/
			ss << monomial->first;
		}
		return ss.str();
	}

	/*
	 * Gives a POSIX regex that represents this polynomial.
	 */
	std::string posixString() const {
        std::stringstream ss;
        std::string posixMonomial;
        bool firstMonomial = true;

        for(auto monomial: monomials_) {
            posixMonomial = monomial.first.posixString();

            if(posixMonomial.size() > 0) {
                if(!firstMonomial) {
                    ss << "|";
                } else {
                    firstMonomial = false;
                }

                ss << posixMonomial;
            }
        }

        return ss.str();
	}

};

template<typename SR> bool NonCommutativePolynomial<SR>::is_commutative = false;
template<typename SR> bool NonCommutativePolynomial<SR>::is_idempotent = false;

/*template <typename SR>
 std::ostream& operator<<(std::ostream& os, const std::map<VarId, SR>& values) {
 for (auto value = values.begin(); value != values.end(); ++value) {
 os << value->first << "→" << value->second << ";";
 }
 return os;
 }*/

template<typename SR>
std::ostream& operator<<(std::ostream& os,
		const std::map<VarId, NonCommutativePolynomial<SR> >& values) {
	for (const auto &key_value : values) {
		os << key_value->first << "→" << key_value->second << ";";
	}
	return os;
}

