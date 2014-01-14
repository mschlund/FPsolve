#pragma once

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <list>
#include <map>
#include <string>

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
		auto tmp_monomial = NonCommutativeMonomial<SR>(idx, variables, srs);
		InsertMonomial(monomials, tmp_monomial, 1);
	}
	;

public:
	NonCommutativePolynomial() = default;

	NonCommutativePolynomial(const NonCommutativePolynomial &p) = default;
	NonCommutativePolynomial(NonCommutativePolynomial &&p) = default;

	/* Create a 'constant' polynomial. */
	NonCommutativePolynomial(const SR &elem) {
		std::vector<std::pair<elemType, int>> idx = { {SemiringType, 0}};
		std::vector<VarId> variables = {};
		std::vector<SR> srs = {elem};
		if(!(elem == SR::null()))
		{
			InsertMonomial(monomials_, idx, variables, srs);
		} // else this is an empty polynomial
	}

	NonCommutativePolynomial(SR &elem) {
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
			result += monomial.first.eval(values) * monomial.second;
		}
		return result;
	}

	/* Evaluate as much as possible (depending on the provided values) and
	 * return the remaining (i.e., unevaluated) monomial. */
	NonCommutativePolynomial<SR> partial_eval(const std::map<VarId, SR> &values) const {
		NonCommutativePolynomial<SR> result = NonCommutativePolynomial<SR>::null();
		for(const auto &monomial : monomials_) {
			result.monomials_.insert(monomial.partial_eval(values));
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
	 * Finds all the constants appearing in this polynomial and appends each new one to the vector.
	 */
	void findAllConstants(std::set<NonCommutativePolynomial<SR>> constants) const {

		// delegate to the monomials
		for(auto &monomial: monomials_) {
			monomial.first.findAllConstants(constants);
		}
	}

	/*
	 * Replaces all constants in the polynomial with their respective VarId mapping.
	 * Used in transformation to Chomsky Normal Form.
	 */
	NonCommutativePolynomial<SR> replaceConstantsWithVariables(std::map<SR, VarId> constantsVariables) const {
		NonCommutativePolynomial<SR> temp = null();

		// delegate to the monomials
		for(auto &monomial: monomials_) {
			temp += monomial.first.replaceConstantsWithVariables(constantsVariables);
		}

		return temp;
	}

	/*
	 * Transforms a polynomial into its Chomsky Normal Form.
	 * Only works for polynomials in which no constants exists (neither as factors nor as summands).
	 */
	NonCommutativePolynomial<SR> chomskyNormalForm(std::map<std::string, VarId> &chomskyVariables,
			std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &chomskyVariableEquations,
			std::map<VarId, SR> variablesToConstants) const {
		NonCommutativePolynomial<SR> temp = null();

		// delegate to the monomials
		for(auto &monomial: monomials_) {
			if(monomial.first.get_degree() != 2) {
				temp += monomial.first.chomskyNormalForm(chomskyVariables, chomskyVariableEquations, variablesToConstants);
			}
		}

		return temp;
	}

	/*
	 * Checks whether this polynomial is productive, depending on the set of variables
	 * that are already known to be productive.
	 */
	bool isProductive(const std::map<VarId, bool> &productiveVariables) const {

		// check if any monomial in this polynomial is productive; that will be enough
		// for the polynomial to be productive
		for(auto monomial: monomials_) {
			if(monomial.first.isProductive(productiveVariables)) {
				return true;
			}
		}

		// if there is no productive monomial, then the polynomial isn't known to be productive so far
		return false;
	}

	/*
	 * Returns a cleaned version of this polynomial, i.e. a version that
	 * had all monomials with unproductive variables eliminated.
	 */
	NonCommutativePolynomial<SR> removeUnproductiveMonomials(const std::map<VarId, bool> &productiveVariables) const {
		NonCommutativePolynomial<SR> cleanPoly = null();

		// add all monomials to the clean polynomial that only contain productive variables
		for(auto monomial: monomials_) {
			if(monomial.first.containsOnlyProductiveVariables(productiveVariables)) {
				cleanPoly += monomial.first * monomial.second;
			}
		}

		return cleanPoly;
	}

	static NonCommutativePolynomial<SR> null() {
		return NonCommutativePolynomial<SR> {};
	}

	static NonCommutativePolynomial<SR> one() {
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
			ss << monomial->second << " * ";
			ss << monomial->first;
		}
		return ss.str();
	}

	std::string posixString() const {
        std::stringstream ss;
        std::string posixMonomial;

        for(auto monomial = monomials_.begin(); monomial != monomials_.end(); monomial++) {
            posixMonomial = monomial.first.posixString();

            if(posixMonomial.size() > 0) {
                if(monomial != monomials_.begin()) {
                    ss << "|";
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

