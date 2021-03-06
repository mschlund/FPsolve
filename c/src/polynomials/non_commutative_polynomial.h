#pragma once

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <list>
#include <map>
#include <string>
#include <queue>
#include <climits>


#include "../semirings/semiring.h"
#include "../semirings/free-semiring.h"

#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/var_degree_map.h"

#include "../utils/string_util.h"

#include "non_commutative_monomial.h"
#include "commutative_polynomial.h"

template <typename SR> class NonCommutativePolynomial;

template <typename SR>
class NonCommutativePolynomialBase : public Semiring<NonCommutativePolynomialBase<SR>,
Commutativity::NonCommutative,
SR::GetIdempotence()> {

protected:
  template <typename S> friend class NonCommutativePolynomialBase;
  friend class NonCommutativePolynomial<SR>;
  std::map<NonCommutativeMonomialBase<SR>,std::uint_fast16_t> monomials_;

  static void InsertMonomial(
      std::map<NonCommutativeMonomialBase<SR>, std::uint_fast16_t> &monomials,
      NonCommutativeMonomialBase<SR> monomial,
      std::uint_fast16_t coeff)
  {
    auto iter = monomials.find(monomial);
    if(iter == monomials.end()) {
      monomials.insert({monomial, 1});
    } else {
      auto tmp = *iter;
      monomials.erase(iter);
      monomials.insert({tmp.first, tmp.second + coeff});
    }
  }

  static void InsertMonomial(
      std::map<NonCommutativeMonomialBase<SR>, std::uint_fast16_t> &monomials,
      std::vector<std::pair<elemType,int>> &idx,
      std::vector<VarId> &variables,
      std::vector<SR> &srs
  )
  {
    auto tmp_monomial = NonCommutativeMonomialBase<SR>(idx, variables, srs);
    InsertMonomial(monomials, tmp_monomial, 1);
  };

  // protected to hide internal structure
  NonCommutativePolynomialBase(std::map<NonCommutativeMonomialBase<SR>,std::uint_fast16_t> && mon) : monomials_(std::move(mon)) {};


public:

  NonCommutativePolynomialBase() = default;

  NonCommutativePolynomialBase(const NonCommutativePolynomialBase &p) = default;
  NonCommutativePolynomialBase(NonCommutativePolynomialBase &&p) = default;

  /* Create a 'constant' polynomial. */
  NonCommutativePolynomialBase(const SR &elem) {
    std::vector<std::pair<elemType, int>> idx = {{SemiringType, 0}};
    std::vector<VarId> variables = {};
    std::vector<SR> srs = {elem};
    if(!(elem == SR::null()))
    {
      InsertMonomial(this->monomials_, idx, variables, srs);
    } // else this is an empty polynomial
  }

  NonCommutativePolynomialBase(SR &elem) {
    std::vector<std::pair<elemType, int>> idx = {{SemiringType, 0}};
    std::vector<VarId> variables = {};
    std::vector<SR> srs = {elem};
    if(!(elem == SR::null()))
    {
      InsertMonomial(this->monomials_, idx, variables, srs);
    } // else this is an empty polynomial
  }

  /* Create a polynomial which consists only of one variable. */
  NonCommutativePolynomialBase(const VarId var) {
    std::vector<std::pair<elemType, int>> idx = {{Variable, 0}};
    std::vector<VarId> variables = {var};
    std::vector<SR> srs = {};
    InsertMonomial(this->monomials_, idx, variables, srs);
  }

  /* initialize with some monomials */
  NonCommutativePolynomialBase(std::initializer_list<NonCommutativeMonomialBase<SR>> monomials) {
    for(auto monomial : monomials) {
      InsertMonomial(this->monomials_, monomial, 1);
    }
  }


  NonCommutativePolynomialBase<SR>& operator=(const NonCommutativePolynomialBase<SR> &p) = default;
  NonCommutativePolynomialBase<SR>& operator=(NonCommutativePolynomialBase<SR> &&p) = default;

  NonCommutativePolynomialBase<SR>& operator+=(const NonCommutativePolynomialBase<SR> &polynomial) {
    if( *this == null())
    {
      *this = polynomial;
      return *this;
    }
    for (const auto &monomial : polynomial.monomials_) {
      // InsertMonomial handles monomials which are already in 'this' polynomial
      InsertMonomial(this->monomials_, monomial.first, monomial.second);
    }

    return *this;
  }

  NonCommutativePolynomialBase<SR>& operator+=(const NonCommutativeMonomialBase<SR> &monomial) {
    InsertMonomial(this->monomials_, monomial, 1);
    return *this;
  }

  NonCommutativePolynomialBase<SR>& operator+(const NonCommutativeMonomialBase<SR> &monomial) const {
    auto result = *this;
    result += monomial;
    return result;
  }

  NonCommutativePolynomialBase<SR>& operator*=(const NonCommutativePolynomialBase<SR> &rhs) {
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
    if (this->monomials_.empty()) {
      return *this;
    } else if (rhs.monomials_.empty()) {
      this->monomials_.clear();
      return *this;
    }

    std::map<NonCommutativeMonomialBase<SR>, std::uint_fast16_t> tmp_monomials;

    // iterate over both lhs und rhs
    for (const auto &lhs_monomial : this->monomials_) {
      for (const auto &rhs_monomial : rhs.monomials_) {
        auto tmp_monomial = lhs_monomial.first * rhs_monomial.first;
        auto tmp_coeff = lhs_monomial.second * rhs_monomial.second;

        InsertMonomial(tmp_monomials, tmp_monomial, tmp_coeff);
      }
    }
    this->monomials_ = std::move(tmp_monomials);
    return *this;
  }

  // multiplying a polynomial with a variable
  NonCommutativePolynomialBase<SR> operator*(const VarId &var) const {
    NonCommutativePolynomialBase result_polynomial;
    for (auto &monomial : this->monomials_) {
      InsertMonomial(result_polynomial.monomials_, monomial.first * var, monomial.second);
    }
    return result_polynomial;
  }

  friend NonCommutativePolynomialBase<SR> operator*(const SR &elem,
      const NonCommutativePolynomialBase<SR> &polynomial) {
    NonCommutativePolynomialBase result_polynomial;
    for (auto &monomial : polynomial.monomials_) {
      result_polynomial.monomials_ += monomial * elem;
    }
    return result_polynomial;
  }

  bool operator==(const NonCommutativePolynomialBase<SR> &polynomial) const {
    return this->monomials_ == polynomial.monomials_;
  }

  /* convert this non-commutative-polynomial to a commutative one */
  CommutativePolynomial<SR> make_commutative() const {
    CommutativePolynomial<SR> result = CommutativePolynomial<SR>::null();

    for(auto const &monomial : this->monomials_)
    {
      result += monomial.first.make_commutative() * monomial.second;
    }

    return result;
  }

  /* calculate the delta for this polynomial with the given data at iteration d
   * de2 is data from newton iterand for d-2
   * dl1 is data from newton update for d-1 */
  SR calculate_delta(
      const ValuationMap<SR> &de2,
      const ValuationMap<SR> &dl1) const {
    SR result = SR::null();
    for(auto const &monomial : this->monomials_) {
      result += monomial.first.calculate_delta(de2, dl1) * monomial.second;
    }

    return result;
  }

  /* calculates the derivative of this polynomial
   * return a schema and the used mapping from X to X[d-1]*/
  std::pair<NonCommutativePolynomialBase<SR>,SubstitutionMap> derivative() const {
    NonCommutativePolynomialBase<SR> result;

    /* get a mapping for all variables X s.t. X -> X[d-1]
     * X[d-1] has to be a fresh variable*/
    auto vars = get_variables();
    SubstitutionMap mapping;
    for(auto const &var : vars) {
      auto tmp = Var::GetVarId();
      mapping.insert({var, tmp});
    }

    /* use the mapping and sum up all derivatives of all monomials */
    for(auto const &monomial : this->monomials_) {
      result += monomial.first.derivative(mapping) * monomial.second;
    }

    return {result, mapping};
  }

  NonCommutativePolynomialBase<SR> differential_at(const ValuationMap<SR>& valuation) const {
    NonCommutativePolynomialBase<SR> result;
    for(auto const &monomial : this->monomials_) {
      result += monomial.first.differential_at(valuation) * monomial.second;
    }
    return result;
  }

  SR eval(const ValuationMap<SR> &values) const {
    SR result = SR::null();
    for (const auto &monomial : this->monomials_) {
      result += monomial.first.eval(values) * monomial.second;
    }
    return result;
  }

  static constexpr Commutativity GetCommutativity() { return Commutativity::NonCommutative; }

  /* Evaluate as much as possible (depending on the provided values) and
   * return the remaining (i.e., unevaluated) monomial. */
  NonCommutativePolynomialBase<SR> partial_eval(const ValuationMap<SR> &values) const {
    NonCommutativePolynomialBase<SR> result = NonCommutativePolynomialBase<SR>::null();
    for(const auto &monomial : this->monomials_) {
      result.monomials_.insert(monomial.partial_eval(values));
    }
    return result;
  }

  /* Variable substitution. */
  NonCommutativePolynomialBase<SR> subst(const SubstitutionMap &mapping) const {
    NonCommutativePolynomialBase result_polynomial;

    for (const auto &monomial : this->monomials_) {
      result_polynomial.monomials_.insert({monomial.first.subst(mapping), monomial.second});
    }
    return result_polynomial;
  }

  static Matrix<SR> eval(const Matrix< NonCommutativePolynomialBase<SR> > &poly_matrix,
      const ValuationMap<SR> &values) {
    const std::vector<NonCommutativePolynomialBase<SR>> &tmp_polynomials = poly_matrix.getElements();
    std::vector<SR> result;
    for (const auto &polynomial : tmp_polynomials) {
      result.emplace_back(polynomial.eval(values));
    }
    return Matrix<SR>{poly_matrix.getRows(), std::move(result)};
  }

  static Matrix<NonCommutativePolynomialBase<SR> > eval(Matrix<NonCommutativePolynomialBase<SR> > poly_matrix,
      ValuationMap<NonCommutativePolynomialBase<SR> > values) {
    const std::vector<NonCommutativePolynomialBase<SR>> &tmp_polynomials = poly_matrix.getElements();
    std::vector<NonCommutativePolynomialBase<SR>> result;
    for (const auto &polynomial : tmp_polynomials) {
      result.emplace_back(polynomial.eval(values));
    }
    return Matrix<NonCommutativePolynomialBase<SR>>{poly_matrix.getRows(), result};
  }

  /* Convert this polynomial to an element of the free semiring.  Note that
   * the provided valuation might be modified with new elements. */
  FreeSemiring make_free(std::unordered_map<SR, VarId, SR> *valuation) const {
    assert(valuation);

    auto result = FreeSemiring::null();
    // convert this polynomial by adding all converted monomials
    for (const auto &monomial : this->monomials_) {
      result += monomial.first.make_free(valuation) * monomial.second;
    }
    return result;
  }

  /* Same as make_free but for matrix form. */
  static Matrix<FreeSemiring> make_free(
      const Matrix< NonCommutativePolynomialBase<SR> > &poly_matrix,
      std::unordered_map<SR, VarId, SR> *valuation) {

    assert(valuation);

    const std::vector< NonCommutativePolynomialBase<SR> > &tmp_polynomials = poly_matrix.getElements();
    std::vector<FreeSemiring> result;

    for (const auto &polynomial : tmp_polynomials) {
      result.emplace_back(polynomial.make_free(valuation));
    }
    return Matrix<FreeSemiring>{poly_matrix.getRows(), std::move(result)};
  }

  Degree get_degree() {
    Degree degree = 0;
    for (auto &monomial : this->monomials_) {
      degree = std::max(degree, monomial.first.get_degree());
    }
    return degree;
  }

  // FIXME: get rid of this...
  // if you use this, make sure, you cache it...
  std::set<VarId> get_variables() const {
    std::set<VarId> vars;
    for(auto const &monomial : this->monomials_)
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

    for(auto &monomial: this->monomials_) {
      sum += monomial.first.getLeadingSR();
    }

    return sum;
  }

  /*
   * Returns the sum of the trailing constant factors of all monomials in this polynomial.
   */
  SR getSumOfTrailingFactors() {
    SR sum = SR::null();

    for(auto &monomial: this->monomials_) {
      sum += monomial.first.getTrailingSR();
    }

    return sum;
  }

  static NonCommutativePolynomialBase<SR> null() {
    return NonCommutativePolynomialBase<SR>{};
  }

  static NonCommutativePolynomialBase<SR> one() {
    return NonCommutativePolynomialBase<SR>{SR::one()};
  }

  static bool is_idempotent;
  static bool is_commutative;

  std::string string() const {
    std::stringstream ss;
    for(auto monomial = this->monomials_.begin(); monomial != this->monomials_.end(); monomial++) {
      if(monomial != this->monomials_.begin())
        ss << " + ";
      //      ss << "[degree: " << monomial->first.get_degree() << "]";
      ss << "[";
      ss << monomial->second << " * ";
      ss << monomial->first;
      ss << "]";
    }
    return ss.str();
  }

  template <typename F>
  auto Map(F fun) const -> NonCommutativePolynomialBase<typename std::result_of<F(SR)>::type> {
    typedef typename std::result_of<F(SR)>::type SR2;

    std::map<NonCommutativeMonomialBase<SR2>,std::uint_fast16_t> result_monomials;

    //TODO: map every monomial and insert it into the result
    for(auto mon_mult : monomials_) {
      NonCommutativePolynomialBase<SR2>::InsertMonomial(result_monomials, mon_mult.first.Map( [&fun](const SR& s){return fun(s);} ), mon_mult.second);
    }


    //FIXME: if fun(pair.second) == zero-element then do not add the monomial!!
    /*std::transform(
           monomials_.begin(), monomials_.end(),
           std::inserter(result_monomials, result_monomials.begin()),
           [&fun](const std::pair< CommutativeMonomial, SR > &pair) {
             return std::make_pair(pair.first, fun(pair.second));
           });
     */

    return NonCommutativePolynomialBase<SR2>(std::move(result_monomials));
  }

};

template <typename SR> bool NonCommutativePolynomialBase<SR>::is_commutative = false;
template <typename SR> bool NonCommutativePolynomialBase<SR>::is_idempotent = false;


template <typename SR>
class NonCommutativePolynomial : public NonCommutativePolynomialBase<SR> {
public:
  using NonCommutativePolynomialBase<SR>::NonCommutativePolynomialBase;

  NonCommutativePolynomial() = default;

  // constructor for implicit conversions
  NonCommutativePolynomial(const NonCommutativePolynomialBase<SR>& p) {
    this->monomials_ = p.monomials_;
  }
};
