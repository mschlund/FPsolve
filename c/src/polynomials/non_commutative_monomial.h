#pragma once

#include <iosfwd>
#include <map>
#include <climits>
#include <ctype.h>
#include <forward_list>

#include "../datastructs/var_degree_map.h"
#include "../semirings/free-semiring.h"

template<typename SR>
class NonCommutativePolynomial;
enum elemType {
	Variable, SemiringType
};

template<typename SR>
class NonCommutativeMonomial {
private:
	friend class NonCommutativePolynomial<SR> ;

	/* a monomial is represented by three vectors
	 * idx_ holds pairs of (type, index), type is either variable or sr
	 *   the index gives absolute position inside the chosen type vector
	 * variables_ holds all variables
	 * srs_ holds all semiring elements */
	std::vector<std::pair<elemType, int>> idx_;
	std::vector<VarId> variables_;
	std::vector<SR> srs_;

	/* Private constructor to not leak the internal data structure. */
	//NonCommutativeMonomial(VarDegreeMap &&vs) : variables_(std::move(vs)) {}
	NonCommutativeMonomial(std::vector<std::pair<elemType, int>> &idx,
			std::vector<VarId> &variables, std::vector<SR> &srs) :
			idx_(idx), variables_(variables), srs_(srs) {
	}

	/*
	 * evaluate all variables except the one at the exceptional position
	 * returns a monomial of the form aXb where a,b are in SR.
	 * TODO: check, that all variables are present in the valuation!
	 */
	NonCommutativeMonomial<SR> eval_all_except(const int except_pos,
			const std::map<VarId, SR>& valuation) const {
		SR a = SR::one();
		SR b = SR::one();

		assert(except_pos < variables_.size());

		/*
		 * go through the idx_-vector,
		 * - evaluate all variables encountered as long as its position is less than except_pos
		 * - on the way: multiply everything (evaluated vars and sr-elems) onto a
		 * - after except_pos has been found: multiply everything encountered onto b
		 */

		SR* prod = &a;
		for (auto p : idx_) {
			if (p.first == elemType::Variable && except_pos == p.second) {
				prod = &b;
				continue;
			}
			if (p.first == elemType::Variable)
				(*prod) *= valuation.at(variables_[p.second]);
			else
				(*prod) *= srs_[p.second];
		}

		std::vector < std::pair<elemType, int>
				> info { std::make_pair(elemType::SemiringType, 0),
						std::make_pair(elemType::Variable, 0), std::make_pair(
								elemType::SemiringType, 1), };
		std::vector<VarId> vars { variables_[except_pos] };
		std::vector<SR> sr { a, b };

		return NonCommutativeMonomial(info, vars, sr);
	}

	/*
	 * Used to recursively generate all monomials for an intersection grammar between a CFG
	 * and a FiniteAutomaton. See intersectionPolynomial(..) for details.
	 *
	 * currentState determines the state in which the FA "currently" would be, i.e. if
	 * we are currently generating all nonterminals <s_i, A, s_{i+1}> of the new grammar, then
	 * currentState would refer us to s_i.
	 *
	 * See Nederhof & Satta, "Probabilistic Parsing", 2008 for details on the generated nonterminals.
	 */
    NonCommutativePolynomial<SR> generateIntersectionMonomials(std::vector<unsigned long> &states,
            std::map<unsigned long,
                std::map<unsigned char, std::forward_list<unsigned long>>> &transitionTable,
            unsigned long currentState,
            unsigned long targetState,
            std::map<unsigned long, unsigned long> &statesToIndices,
            std::map<VarId, unsigned long> &oldVariablesToIndices,
            std::vector<std::vector<std::vector<VarId>>> &newVariables,
            unsigned long monomialFactorIndex) const {

        NonCommutativePolynomial<SR> result = NonCommutativePolynomial<SR>::null();

        // if the current position in the monomial is a semiring element, we see where the transitions using the element
        // starting at currentStartIndex can take us and recurse from there
        // otherwise, we replace the variable with all replacement variables and recursive from there
        if(idx_[monomialFactorIndex].first == SemiringType) { // we have a semiring element

            // check if we hit the last factor of the monomial; we stop the recursion here
            if(monomialFactorIndex == idx_.size() - 1) {
//                std::cout << "last index SR" << std::endl;

                // if the last element in the monomial is epsilon, then the only way the generated monomial is valid
                // is if the target state is equal to the one we are currently at
                if(srs_[idx_[monomialFactorIndex].second] == SR::one()) {
                    if(currentState == targetState) {
                        result = NonCommutativePolynomial<SR>::one();
                    }
                } else { // otherwise, see if we can reach a target state

                    // get a string representation of the semiring element, see where that string takes us in the FA
                    // -----------------------------------------------------
                    // -----------------------------------------------------
                    // -----------------------------------------------------
                    // note that due to this bit here, we can only work with grammars where we have elements of SIGMA*
                    // as semiring factors; if we want to use regular expressions over SIGMA here, the following
                    // needs to be adjusted
                    // -----------------------------------------------------
                    // -----------------------------------------------------
                    // -----------------------------------------------------
                    std::string transitionsToDo = srs_[idx_[monomialFactorIndex].second].string();

                    // for later debug purposes
//                    std::cout << "transitionsToDo during intersection:\t" + transitionsToDo << std::endl;

                    std::set<unsigned long> reachable, temp;
                    reachable.insert(currentState);

                    for(int i = 0; i < transitionsToDo.size(); i++) {
                        for(auto &intermediateState: reachable) {
//                            std::cout << "intermediate state:\t" << intermediateState << std::endl;
//                            std::cout << "transitionsToDo[i]:\t" << transitionsToDo[i] << std::endl;
                            for(unsigned long &nextIntermediateState: transitionTable[intermediateState][transitionsToDo[i]]) {
                                temp.insert(nextIntermediateState);
                            }
                        }

                        reachable.swap(temp);
                        temp.clear();
                    }

                    // if we can reach the intended state with the given transitions, we allow the semiring element
                    if(reachable.count(targetState) != 0) {
                        result = NonCommutativePolynomial<SR>(srs_[idx_[monomialFactorIndex].second]);
                    }
                }
            } else { // if this is not the last factor of the monomial, continue the recursion
//                std::cout << "not last index SR" << std::endl;

                // if the semiring element is 1, then the FA cannot make a move since it doesn't have epsilon transitions,
                // so we advance in the monomial but leave the currentState untouched
                if(srs_[idx_[monomialFactorIndex].second] == SR::one()) {
                    result = generateIntersectionMonomials(states, transitionTable, currentState, targetState,
                                 statesToIndices, oldVariablesToIndices, newVariables, (monomialFactorIndex + 1));
                } else { // if the semiring element is not the 1 element, get a string representation of it
                         // and see where that string takes us in the FAstateTransitions

                    std::string transitionsToDo = srs_[idx_[monomialFactorIndex].second].string();

                    // for later debug purposes
//                    std::cout << "transitionsToDo during intersection:\t" + transitionsToDo << std::endl;

                    std::set<unsigned long> reachable, temp;
                    reachable.insert(currentState);

                    for(int i = 0; i < transitionsToDo.size(); i++) {
                        for(auto &intermediateState: reachable) {
//                            std::cout << "intermediate state:\t" << intermediateState << std::endl;
//                            std::cout << "transitionsToDo[i]:\t" << transitionsToDo[i] << std::endl;
                            for(unsigned long &nextIntermediateState: transitionTable[intermediateState][transitionsToDo[i]]) {
                                temp.insert(nextIntermediateState);
                            }
                        }

                        reachable.swap(temp);
                        temp.clear();
                    }

                    // we proceed with generating the monomial starting at the next factor and all reachable states
                    NonCommutativePolynomial<SR> suffixPolynomial = NonCommutativePolynomial<SR>::null();
                    for(auto &nextState: reachable) {
                        suffixPolynomial += generateIntersectionMonomials(states, transitionTable, nextState, targetState,
                                                statesToIndices, oldVariablesToIndices, newVariables, (monomialFactorIndex + 1));
                    }

                    result = NonCommutativePolynomial<SR>(srs_[idx_[monomialFactorIndex].second]) * suffixPolynomial;
                }
            }
        } else { // we have a variable at the current index

            // check if we hit the last factor of the monomial; we stop the recursion here
            if(monomialFactorIndex == idx_.size() - 1) {
//                std::cout << "last index variable" << std::endl;
//                std::cout << "statesToIndices[currentState]: " << statesToIndices[currentState] << std::endl;
                result = NonCommutativePolynomial<SR>(newVariables[statesToIndices[currentState]][statesToIndices[targetState]]
                            [oldVariablesToIndices.at(variables_[idx_[monomialFactorIndex].second])]);
            } else { // if this is not the last factor of the monomial, continue the recursion
//                std::cout << "not last index variable" << std::endl;
                for(auto &nextState: states) {

                    // we want to replace some nonterminal X by <currentState, X, nextState> for all nextState;
                    // after that we need to append all possible monomials generated from the suffix of this monomial
                    // that starts one symbol after X and uses nextState
//                    std::cout << "statesToIndices[currentState]: " << statesToIndices[currentState] << std::endl;
                    result += NonCommutativePolynomial<SR>(newVariables[statesToIndices[currentState]][statesToIndices[nextState]]
                                  [oldVariablesToIndices.at(variables_[idx_[monomialFactorIndex].second])]) *
                              generateIntersectionMonomials(states, transitionTable, nextState, targetState,
                                  statesToIndices, oldVariablesToIndices, newVariables, (monomialFactorIndex + 1));
                }
            }
        }

        return result;
    }

public:

	/* Since we don't manage any resources, we can simply use the default
	 * constructors and assignment operators. */
	NonCommutativeMonomial() = default;
	NonCommutativeMonomial(const NonCommutativeMonomial &m) = default;
	NonCommutativeMonomial(NonCommutativeMonomial &&m) = default;
	NonCommutativeMonomial& operator=(const NonCommutativeMonomial &m) = default;
	NonCommutativeMonomial& operator=(NonCommutativeMonomial &&m) = default;

	NonCommutativeMonomial(std::initializer_list<std::pair<elemType,int>> idx,
			std::initializer_list<VarId> variables,
			std::initializer_list<SR> srs) :
	idx_(idx), variables_(variables), srs_(srs_) {}

	/* std::vector seems to be a neutral data type and does not leak internal
	 * data structure. */
	/*    NonCommutativeMonomial(std::vector< std::pair<elemType, int> > vs) {
	 for (auto var_degree : vs) {
	 variables_.Insert(var_degree.first, var_degree.second);
	 }
	 }
	 */

	/* Multiply two monomials and return a result in normalform if both monomials
	 * already are in normal form */
	NonCommutativeMonomial operator*(const NonCommutativeMonomial &monomial) const {
		auto tmp_idx = idx_;
		auto tmp_variables = variables_;
		auto tmp_srs = srs_;
		bool normalform = true;
		unsigned int offset_variables = tmp_variables.size();
		unsigned int offset_srs = tmp_srs.size();

		// assume both monomials are in normal form
		// then the only position where this assumption might
		// be violated is where both are concatenated
		if(tmp_idx.back().first == SemiringType && monomial.idx_.front().first == SemiringType)
		{
			tmp_idx.pop_back();
			offset_srs--;
			normalform = false;
		}

		for (auto &p : monomial.idx_) {
			if(p.first == Variable)
			tmp_idx.push_back( {Variable, p.second + offset_variables});
			else if (p.first == SemiringType)
			tmp_idx.push_back( {SemiringType, p.second + offset_srs});
		}

		for (auto v : monomial.variables_)
		tmp_variables.push_back(v);

		//for (auto s : monomial.srs_)
		for(auto s = monomial.srs_.begin(); s != monomial.srs_.end(); s++)
		if(!normalform && s == monomial.srs_.begin())
		{
			auto elem = tmp_srs.back();
			tmp_srs.pop_back();
			tmp_srs.push_back(elem * (*s));
			normalform = true;
		} else {
			tmp_srs.push_back(*s);
		}

		//return NonCommutativeMonomial(std::move(tmp_idx), std::move(tmp_variables), std::move(tmp_srs));
		return NonCommutativeMonomial(tmp_idx, tmp_variables, tmp_srs);
	}

	/* Multiply a monomial with a variable. */
	NonCommutativeMonomial operator*(const VarId &var) const {
		auto tmp_idx = idx_;
		auto tmp_variables = variables_;
		auto tmp_srs = srs_;
		tmp_idx.push_back(std::pair<elemType, int>(Variable, tmp_variables.size()));
		tmp_variables.push_back(var);

		return NonCommutativeMonomial(tmp_idx, tmp_variables, tmp_srs);
	}

	/* Multiply a monomial with a semiring element and return it in normal form. */
	NonCommutativeMonomial operator*(const SR &sr) const {
		auto tmp_idx = idx_;
		auto tmp_variables = variables_;
		auto tmp_srs = srs_;

		// assume the monomial is in normal form
		if(tmp_idx.back().first() == SemiringType)
		{
			auto elem = tmp_srs.back();
			tmp_srs.pop_back();
			tmp_srs.push_back(elem * sr);
		} else {
			tmp_idx.push_back(std::pair<elemType, int>(SemiringType, tmp_srs.size()));
			tmp_srs.push_back(sr);
		}

		return NonCommutativeMonomial(tmp_idx, tmp_variables, tmp_srs);
	}

	/* convert this monomial into an commutative one
	 * return a Polynomial, because the semiring element is not saved
	 * in the commutative version of the monomial but in the polynomial */
	NonCommutativePolynomial<SR> make_commutative() const {
		NonCommutativePolynomial<SR> result_polynomial;

		for(auto const &sr : srs_) {
			result_polynomial *= sr;
		}

		for(auto const &var : variables_) {
			result_polynomial *= var;
		}

		return result_polynomial;
	}

	/* derivation function which is used in the polynomial derivative function.
	 * for the variables for the d-1-th iterand we use the given map 'substitution' */
	NonCommutativePolynomial<SR> derivative(const std::map<VarId, VarId> &substitution) const {
		NonCommutativePolynomial<SR> result; // empty polynomial
		auto subst_monomial = this->subst(substitution);// substitute all variables
		for(unsigned int position = 0; position < variables_.size(); position++) {
			// the variable at the position is the variable which should not be touched
			auto tmp = subst_monomial;
			tmp.variables_.at(position) = this->variables_.at(position);// therefore restore this one...
			result += tmp;
		}
		return result;
	}

	/*
	 * compute the linearization of the polynomial at the given point "valuation"
	 * to this end, we linearize every monomial and sum them up
	 */
	NonCommutativePolynomial<SR> differential_at(const std::map<VarId, SR> &valuation) const {
		NonCommutativePolynomial<SR> result; // empty polynomial = 0 (constant monomials have differential f(x)=0)
        for(unsigned int position = 0; position < variables_.size(); position++) {
            // the variable at the position is the variable which should not be touched, all others will be evaluated
            result += this->eval_all_except(position, valuation);
        }
		return result;
	}

	SR calculate_delta_helper(
			const std::vector<bool> &permutation,
			const std::map<VarId, SR> &de2, // [d-2], true
			const std::map<VarId, SR> &dl1// (d-1), false
	) const {
		SR tmp = SR::one();
		for(auto const &p : idx_) {
			if(p.first == Variable) {
				if(permutation.at(p.second) == true) { // use [d-2]
					auto value_iter = de2.find(variables_.at(p.second));
					assert(value_iter != de2.end());
					tmp *= value_iter->second;
				} else { // if (permutation.at(p.second) == false) // use (d-1)
					auto value_iter = dl1.find(variables_.at(p.second));
					assert(value_iter != dl1.end());
					tmp *= value_iter->second;
				}
			} else if (p.first == SemiringType)
			tmp *= srs_.at(p.second);
		}
		return tmp;
	}

	SR calculate_delta(
			const std::map<VarId, SR> &de2, // [d-2], true
			const std::map<VarId, SR> &dl1// (d-1), false
	) const {
		SR result = SR::null();

		/* outer loop handles the different cases (trees with exactly n-times dim == d-1 )
		 * start with n = 2, which means, exactly 2 children have dimensions exactly d-1 */
		for(unsigned int n = 2; n <= variables_.size(); n++)
		{
			/* order of a vector of bools is [false, true] < [true, false] */
			std::vector<bool> permutation(n, false); // these are the (d-1) elements
			std::vector<bool> permutation2(variables_.size()-n, true);// these are the [d-2] elements
			permutation.insert(permutation.end(), permutation2.begin(), permutation2.end());
			do {
				result += calculate_delta_helper(permutation, de2, dl1);
			}while(std::next_permutation(permutation.begin(), permutation.end()));
		}

		return result;
	}

	/* Evaluate the monomial given the map from variables to values. */
	SR eval(const std::map<VarId, SR> &values) const {
		auto result = SR::one();

		for (const auto &p : idx_)
		{
			if (p.first == Variable)
			{
				auto value_iter = values.find(variables_.at(p.second));
				/* All variables should be in the values map. */
				assert(value_iter != values.end());
				result *= value_iter->second;
			}
			else if (p.first == SemiringType)
			result *= srs_.at(p.second);
		}
//		std::cout << "monomial eval result:\t\t" << result.string() << std::endl;

		return result;
	}

	/* Partially evaluate the monomial. */
	NonCommutativeMonomial partial_eval(
			const std::map<VarId, SR> &values) const {

		NonCommutativeMonomial result_monomial;

		for(auto p : idx_)
		{
			if (p.first == Variable)
			{
				auto value_iter = values.find(variables_.at(p.second));
				if (value_iter == values.end())
				{ /* Variable not found in the mapping, so keep it. */
					result_monomial.idx_.push_back(std::pair<elemType, int>(Variable, result_monomial.variables_.size()));
					result_monomial.variables_.push_back(variables_.at(p.second));
				} else
				{ /* Variable found, use it for evaluation */
					if(result_monomial.idx_.size() > 0 && result_monomial.idx_.back().first == SemiringType)
					{
						// multiplying sr-elements to achieve normal form
						auto elem = result_monomial.srs_.back();
						result_monomial.srs_.pop_back();
						result_monomial.srs_.push_back(elem * value_iter->second);

					} else {
						result_monomial.idx_.push_back(std::pair<elemType, int>(SemiringType, result_monomial.srs_.size()));
						result_monomial.srs_.push_back(value_iter->second);
					}
				}

			} else if (p.first == SemiringType)
			{
				if(result_monomial.idx_.size() > 0 && result_monomial.idx_.back().first == SemiringType)
				{
					// multiplying sr-elements to achieve normal form
					auto elem = result_monomial.srs_.back();
					result_monomial.srs_.pop_back();
					result_monomial.srs_.push_back(elem * srs_.at(p.second));

				} else {
					result_monomial.idx_.push_back(std::pair<elemType, int>(SemiringType, result_monomial.srs_.size()));
					result_monomial.srs_.push_back(srs_.at(p.second));
				}
			}
		}

		return result_monomial;
	}

	/* Variable substitution. */
	NonCommutativeMonomial subst(const std::map<VarId, VarId> &mapping) const {
		VarDegreeMap tmp_variables;

		auto result_monomial = *this; // copy it to work with it

		for(auto p : result_monomial.idx_)
		{
			if(p.first == Variable)
			{
				auto old_new_iter = mapping.find(variables_.at(p.second));
				if(old_new_iter != mapping.end())
				{ // substitute
					result_monomial.variables_.at(p.second) = old_new_iter->second;
				} else
				{ // do nothing
					continue;
				}
			} else
			{ // do nothing
				continue;
			}
		}

		return result_monomial;
	}

	/* Convert this monomial to an element of the free semiring. */
	FreeSemiring make_free(std::unordered_map<SR, VarId, SR> *valuation) const {
		FreeSemiring result = FreeSemiring::one();

		for (auto p : idx_) {
			if (p.first == Variable) {
				result *= FreeSemiring(variables_.at(p.second));
			} else if (p.first == SemiringType) {
				{
					auto tmp_sr = srs_.at(p.second);
					if(tmp_sr == SR::null()) {
						assert(false); //coefficients in the monomial are always != 0.. so this should not happen :)
					} else if (tmp_sr == SR::one()) {
						// does not do anything
					} else {
						auto value_iter = valuation->find(tmp_sr);
						if (value_iter == valuation->end()) {
							/* Use a fresh constant - the constructor of Var::getVar() will take
							 * care of this. */
							VarId tmp_var = Var::GetVarId();
							FreeSemiring tmp_var_free {tmp_var};
							valuation->emplace(tmp_sr, tmp_var);
							result *= tmp_var_free;
						} else {
							// there is already a variable for this element, use it
							result *= value_iter->second;
						}
					}
				}
			}
		}

		return result;
	}

	bool operator<(const NonCommutativeMonomial &rhs) const {
		// lexicographic ordering
		if(idx_ != rhs.idx_) return idx_ < rhs.idx_;
		if(variables_ != rhs.variables_) return variables_ < rhs.variables_;
		if(srs_ != rhs.srs_) return srs_ < rhs.srs_;

		// they are equal
		return false;
	}

	bool operator==(const NonCommutativeMonomial &rhs) const {
		return
		idx_ == rhs.idx_ &&
		variables_ == rhs.variables_ &&
		srs_ == rhs.srs_;
	}

	Degree get_degree() const {
		return variables_.size();
	}

	// FIXME: modify or remove
	std::set<VarId> get_variables() const {
		std::set<VarId> set;
		for (auto var : variables_) {
			set.insert(var);
		}
		return set;
	}

	/*
	 * Extracts the terminal letters used in this monomial. Currently only recognizes alphanumeric characters.
	 */
	std::set<unsigned char> get_terminals() const {
	    std::set<unsigned char> terminals;

	    for(auto &semiring_element: srs_) {
	        if(!(semiring_element == SR::null()) && !(semiring_element == SR::one())) {
	            std::string regex = semiring_element.string();

	            for(int i = 0; i < regex.size(); i++) {
	                if(isalnum(regex[i])) {
	                    terminals.insert(regex[i]);
	                }
	            }
	        }
	    }

	    return terminals;
	}

	/*
	 * If the monomial has form xYz where x is an element of the semiring,
	 * Y is a monomial over the semiring, and z is an element of the semiring,
	 * this method gives back x wrapped in a monomial.
	 */
	SR getLeadingSR() const {

	    if(get_degree() == 0) {
	        return srs_[0];
	    }

		SR leadingSR = SR::one();

		// give me a copy of the semiring factor on the leading side of the monomial
		for(int i = 0; i < idx_.size(); i++) {

			// multiply as long as there was no variable encountered
			if(idx_.at(i).first == SemiringType) {
				leadingSR = leadingSR * srs_.at(idx_.at(i).second);
			} else {
				break; // break on the first variable
			}
		}

		return leadingSR;
	}

	/*
	 * If the monomial has form xYz where x is an element of the semiring,
	 * Y is a monomial over the semiring, and z is an element of the semiring,
	 * this method gives back z wrapped in a monomial.
	 */
	SR getTrailingSR() const {

	    if(get_degree() == 0) {
	        return srs_[0];
	    }

		SR trailingSR = SR::one();

		// give me a copy of the semiring factor on the trailing side of the monomial
		for(int i = idx_.size() - 1; i >= 0; i--) {

			// multiply as long as no variable is encountered
			if(idx_.at(i).first == SemiringType) {
				trailingSR = srs_.at(idx_.at(i).second) * trailingSR;
			} else {
				break; // break on the first variable
			}
		}

		return trailingSR;
	}

	/*
	 * Quick check to see whether this monomial is just an epsilon.
	 */
	bool isEpsilonMonomial() const {
	    return (get_degree() == 0) && (getLeadingSR() == SR::one());
	}

	/*
	 * Finds all the constants appearing in this monomial if it has degree at least 1 if checkLinearMonomials is true,
	 * otherwise if it has degree at least 2.
	 */
	void findConstantsInNonterminalMonomials(std::set<SR> &constants, bool checkLinearMonomials) const {


	    // we only replace the constants in monomials of length >= 2 with new variables; we also ignore linear
	    // monomials depending on checkLinearMonomials
	    if(get_degree() == 0 || (!checkLinearMonomials && get_degree() == 1)
	            || (get_degree() == 1 && getLeadingSR() == SR::one() && getTrailingSR() == SR::one())) {
	        return;
	    }

		for(int i = idx_.size() - 1; i >= 0; i--) {
			//we only want to extract the constants from the monomial; we also don't care about epsilon
			if((idx_.at(i).first == SemiringType) && !(srs_.at(idx_.at(i).second) == SR::one())) {

				//std::set.insert(..) does its own checking whether the element already is in the set
				constants.insert(srs_.at(idx_.at(i).second));
			}
		}

		return;
	}

    void findLowerLinearTerms(std::set<NonCommutativeMonomial<SR>> &linearMonomialsOfLowerOrder,
            std::map<VarId, int> &varToComponent,
            int component) const {

        if((get_degree() != 1) || (varToComponent[variables_[0]] == component)) {
            return;
        }

        linearMonomialsOfLowerOrder.insert(*this);
    }

	/*
	 * Replaces all constants in the polynomial with their respective VarId mapping.
	 *
	 * WARNING: This function will fail if there is an SR element in this monomial that does
	 * not have a mapped VarId. Will not change the monomial consisting exclusively of epsilon.
	 */
	NonCommutativePolynomial<SR> replaceConstants(std::map<SR, VarId> &constantsToVariables) const {
		NonCommutativePolynomial<SR> temp;

		bool firstFactorIsSR = (idx_.at(0).first == SemiringType);
		bool firstFactorIsEpsilon = false;

		if(firstFactorIsSR) {
		    firstFactorIsEpsilon = (srs_.at(idx_.at(0).second) == SR::one());
		}

		if((idx_.size() == 1) && firstFactorIsSR && firstFactorIsEpsilon) {
		    return NonCommutativePolynomial<SR>::one();
		}

		int startIndex = 1;
		if(firstFactorIsSR) {

		    // we don't want unnecessary epsilon-factors
		    if(firstFactorIsEpsilon) {
		        temp = NonCommutativePolynomial<SR>(variables_.at(idx_.at(1).second));
		        startIndex = 2;
		    } else {
	            temp = NonCommutativePolynomial<SR>(constantsToVariables[srs_.at(idx_.at(0).second)]);
		    }
		} else {
			temp = NonCommutativePolynomial<SR>(variables_.at(idx_.at(0).second));
		}

		for(int i = startIndex; i < idx_.size(); i++) {

			// multiply as long as there was no variable encountered
			if(idx_.at(i).first == SemiringType) {

			    // ignore all epsilon factors
			    if(!(srs_.at(idx_.at(i).second) == SR::one())) {
	                temp *= NonCommutativePolynomial<SR>(constantsToVariables[srs_.at(idx_.at(i).second)]);
			    }
			} else {
				temp *= variables_.at(idx_.at(i).second);
			}
		}

//		std::cout << "poly in CNF:\t"+ temp.string() << std::endl;
		return temp;
	}

	/*
	 * This produces the polynomial (that derives from this monomial) that is generated while
	 * building the intersection grammar of a CFG and an FA.
	 *
	 * See NonCommutativePolynomial<SR>::intersectionPolynomial(..) for further documentation.
	 *
     * WARNING: if you use an SR where the string representations of elements are NOT strings over
     * the English alphabet (e.g. they may contain 1s to represent epsilons, or they may use a
     * different character set), then this function will most likely not do what you want it to do.
	 */
	NonCommutativePolynomial<SR> intersectionPolynomial(std::vector<unsigned long> &states,
	            std::map<unsigned long, std::map<unsigned char, std::forward_list<unsigned long>>> &transitionTable,
	            unsigned long startState,
	            unsigned long targetState,
	            std::map<unsigned long, unsigned long> &statesToIndices,
	            std::map<VarId, unsigned long> &oldVariablesToIndices,
	            std::vector<std::vector<std::vector<VarId>>> &newVariables) const {

        return (NonCommutativePolynomial<SR>::one() *
                generateIntersectionMonomials(states, transitionTable, startState,
                        targetState, statesToIndices, oldVariablesToIndices, newVariables, 0));
	}

	/*
	 * Binarizes this monomial, as one would when calculating a CNF of a grammar. Only works on monomials
	 * that have no terminals (i.e. semiring elements) as factors; for all other monomials (i.e. those with
	 * terminals as factors), this will return NonCommutativePolynomial<SR>::null().
	 *
	 * If the monomial consists of exactly one or exactly two variables, the resulting polynomial will not differ from it.
	 *
     * The variables and productions that are introduced during that process will be stored in binarizationVariables
     * and binarizationVariablesEquations, respectively.
	 */
	NonCommutativePolynomial<SR> binarize(std::map<std::string, VarId> &binarizationVariables,
			std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &binarizationVariablesEquations,
			std::map<VarId, SR> &variablesToConstants) const {

	    // make sure there are no terminal factors in the monomial
	    if(idx_.size() != variables_.size()) {
	        return NonCommutativePolynomial<SR>::null();
	    }

		NonCommutativePolynomial<SR> temp = NonCommutativePolynomial<SR>::one();
		NonCommutativePolynomial<SR> suffix;

		// we have at least two variables; we need to shorten that to exactly two by
		// introducing new variables where needed
		int variablesProcessed = 0;
		int index = 0;
		VarId var;

		// loop over all suffixes to see which need new variables introduced for them
		while(variablesProcessed < get_degree() - 1) {

			// if this is the first iteration, we start with the leftmost variable
			if(variablesProcessed == 0) {
				index = variables_.size() - 1;

				temp = NonCommutativePolynomial<SR>(variables_.at(index));
				suffix = NonCommutativePolynomial<SR>(variables_.at(index));
			} else { /* if this is not the first iteration, that means we already remembered a variable
					  * to produce the current suffix of the monomial; multiply the next variable from left,
					  * check if some variable already maps to the current suffix, if none does introduce
					  * a new variable that produces the product of the now two variables we have, and
					  * make this new variable the starting point for the next iteration; also remember
					  * the suffix so far so we can check whether maybe we already have a variable that
					  * produces that suffix in the next iteration
					  */
				temp = NonCommutativePolynomial<SR>(variables_.at(index)) * temp; // multiplication from left is relevant -> noncommutative!
				suffix = NonCommutativePolynomial<SR>(variables_.at(index)) * suffix;

				// check if we already have a variable for the current suffix; if we do, proceed from there;
				// else, introduce a new one
				if(binarizationVariables.find(suffix.string()) != binarizationVariables.end()) {
					temp = NonCommutativePolynomial<SR>(binarizationVariables[suffix.string()]);
				} else {
					var = Var::GetVarId();
					binarizationVariablesEquations.push_back(std::make_pair(var, temp));
					binarizationVariables.insert(std::make_pair(suffix.string(), var));
					temp = NonCommutativePolynomial<SR>(var);
				}
			}

			index--; // we process the monomial from right to left
			variablesProcessed++;
		}


		return NonCommutativePolynomial<SR>(variables_.at(0)) * temp;
	}

	/*
	 * Used when finding the shortest word derivable from a CFG. This method tells us whether
	 * the algorithm already knows a terminal string that derives from this monomial.
	 */
	bool findLengthOfDerivableStrings(std::map<VarId, unsigned long> &lengthOfShortestTerminal,
	        std::map<VarId, NonCommutativeMonomial<SR>> &productionsForShortestWords,
	        VarId &lhsOfProduction) const {

	    // if this is a terminal, remember its length if we don't already know it; otherwise, see if we already
	    // know the length of some derivable string and update it if necessary
	    if(variables_.size() == 0){

	        // if this is the shortest terminal the lefthand side of the production can produce, update
	        if(getLeadingSR() != SR::null() && lengthOfShortestTerminal[lhsOfProduction] > getLeadingSR().string().size()) {
	            lengthOfShortestTerminal[lhsOfProduction] = getLeadingSR().string().size();
	            productionsForShortestWords[lhsOfProduction] = *this;
	            return true;
	        }

	        return false;
	    } else {

	        unsigned long shortestMonomialLength = 0;

	        // check for all variables whether we know how to derive a terminal
	        for(VarId var: variables_) {
	            if(lengthOfShortestTerminal.at(var) == ULONG_MAX) {
	                shortestMonomialLength = ULONG_MAX;
	                break;
	            }

	            shortestMonomialLength += lengthOfShortestTerminal.at(var);
	        }

	        // if we found a shorter word than previously known for the nonterminal this production
	        // belongs to, update
	        if(shortestMonomialLength < ULONG_MAX) {

	            for(auto terminal: srs_) {
	                shortestMonomialLength += terminal.string().size();
	            }

                if(shortestMonomialLength < lengthOfShortestTerminal.at(lhsOfProduction)) {
                    // update the map with the lengths of shortest terminals
                    productionsForShortestWords[lhsOfProduction] = *this;
                    lengthOfShortestTerminal[lhsOfProduction] = shortestMonomialLength;
                    return true;
	            }
	        }

	        return false;
	    }
	}

	/*
	 * Derives the shortest derivable terminal for this monomial, assuming the info in productionsForShortestWords
	 * is correct.
	 */
	SR shortestDerivableTerminal(std::map<VarId, NonCommutativeMonomial<SR>> &productionsForShortestWords) const {
	    SR terminal = SR::one();

	    if(get_degree() == 0) {
	        terminal = getLeadingSR();
	    } else {
	        for(auto &indexpair: idx_) {
	            if(indexpair.first == SemiringType) {
	                terminal *= srs_.at(indexpair.second);
	            } else {
	                assert(productionsForShortestWords.find(variables_.at(indexpair.second)) != productionsForShortestWords.end());
	                terminal *= productionsForShortestWords[variables_.at(indexpair.second)]
	                                                        .shortestDerivableTerminal(productionsForShortestWords);
	            }
	        }
	    }

	    return terminal;
	}

	/*
	 * Checks whether this monomial is productive, depending on the set of variables
	 * that are already known to be productive.
	 */
	bool isProductive(const std::map<VarId, bool> &productiveVariables) const {

		// if this monomial has no variables, then it is productive since
		// it represents an element of the semiring
		if(get_degree() == 0) {
			return true;
		}

		// check if all variables of this monomial are productive
		for(auto &variable: variables_) {

			// if there is any variable that is not known to be productive in
			// the current monomial, then it isn't productive
			if(!productiveVariables.at(variable)) {
				return false;
			}
		}

		// if all variables are productive, then so is the monomial
		return true;
	}

	/*
	 * See NonCommutativePolynomial<SR>::componentIsSquarable(..) for commentary.
	 */
    bool componentIsSquarable(std::map<VarId, int> &varToComponent, int component) const {
        int count = 0;

        // if the monomial contains any two variables from the same component,
        // then we can duplicate all variables in the component
        for(auto var: variables_) {
            if(varToComponent[var] == component) {
                count++;

                if(count >= 2) {
                    return true;
                }
            }
        }

        return false;
    }

    void mapQuadraticLHStoRHS(std::map<int, std::map<int, std::set<int>>> &quadraticLHStoRHS,
                std::map<VarId, int> &varToComponent, int component) const {
        if(get_degree() != 2) {
            return;
        }

        int lhsComponent = varToComponent[variables_[0]];
        int rhsComponent = varToComponent[variables_[1]];

        if((lhsComponent != component) && (rhsComponent != component)) {
            quadraticLHStoRHS[component][lhsComponent].insert(rhsComponent);
//            for(auto entry: quadraticLHStoRHS[component][lhsComponent]) {
//                std::cout << "entry in quadraticLHStoRHS[" << component << "][" << lhsComponent << "]: " << rhsComponent << std::endl;
//            }
        }
    }

    void calculateLowerComponentVariables(
                std::map<int, std::set<int>> &lhsLowerComponentVariables,
                std::map<int, std::set<int>> &rhsLowerComponentVariables,
                std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>>> &components,
                std::map<VarId, int> &varToComponent,
                int component) const {

        if(get_degree() != 2) {
            return;
        }

        int lhsComponent = varToComponent[variables_[0]];
        int rhsComponent = varToComponent[variables_[1]];

        if(lhsComponent == component) {
            rhsLowerComponentVariables[component].insert(rhsComponent);
        } else if(rhsComponent == component) {
            lhsLowerComponentVariables[component].insert(lhsComponent);
        }
    }

    void calculateSameComponentLetters(
            std::map<int, std::set<unsigned char>> &lhsSameComponentLetters,
            std::map<int, std::set<unsigned char>> &rhsSameComponentLetters,
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>>> &components,
            std::map<VarId, int> &varToComponent,
            int component) const {

        if(get_degree() != 1) {
            return;
        }

        int varComponent = varToComponent[variables_[0]];

        if(varComponent == component) {
            std::string lhsString = getLeadingSR().string();
            std::string rhsString = getTrailingSR().string();

            for(int i = 0; i < lhsString.size(); i++) {
                if(isalnum(lhsString[i])) {
                    lhsSameComponentLetters[component].insert(lhsString[i]);
                } else if (lhsString[i] == '[') {
                    unsigned char letter, startOfRange, endOfRange;
                    i++;
                    while(lhsString[i] != ']') {
                        if(lhsString[i+1] == '-') {
                            startOfRange = lhsString[i];
                            endOfRange = lhsString[i+2];
                            i += 3;

                            for(letter = startOfRange; letter <= endOfRange; letter++) {
                                lhsSameComponentLetters[component].insert(letter);
                            }
                        } else {
                            lhsSameComponentLetters[component].insert(lhsString[i]);
                            i++;
                        }
                    }
                }
            }

            for(int i = 0; i < rhsString.size(); i++) {
                if(isalnum(rhsString[i])) {
                    rhsSameComponentLetters[component].insert(rhsString[i]);
                } else if (rhsString[i] == '[') {
                    unsigned char letter, startOfRange, endOfRange;
                    i++;
                    while(rhsString[i] != ']') {
                        if(rhsString[i+1] == '-') {
                            startOfRange = rhsString[i];
                            endOfRange = rhsString[i+2];
                            i += 3;

                            for(letter = startOfRange; letter <= endOfRange; letter++) {
                                rhsSameComponentLetters[component].insert(letter);
                            }
                        } else {
                            rhsSameComponentLetters[component].insert(rhsString[i]);
                            i++;
                        }
                    }
                }
            }
        }
    }

	/*
	 * Checks if this monomial is a chain production, that is if it contains
	 * exactly one variable and nothing else.
	 */
	bool isChainProduction() const {
	    return (get_degree() == 1 && getLeadingSR() == SR::one() && getTrailingSR() == SR::one());
	}

	std::string string() const {

		std::stringstream ss;
		//for(auto &p : idx_) {
		for(auto p = idx_.begin(); p != idx_.end(); p++) {
			if(p != idx_.begin())
			ss << " * ";

			if(p->first == Variable)
			ss << variables_.at(p->second);
			else
			ss << srs_.at(p->second);
		}
		return std::move(ss.str());
	}

	std::string posixString() const {
        std::stringstream ss;
        for(auto p = idx_.begin(); p != idx_.end(); p++) {
            if(p->first == Variable) {
                ss << variables_.at(p->second);
            } else {
                ss << srs_.at(p->second);
            }
        }
        return std::move(ss.str());
	}
};

template<typename SR>
std::ostream& operator<<(std::ostream &out,
		const NonCommutativeMonomial<SR> &monomial) {
	return out << monomial.string();
}
