#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <list>
#include <map>
#include <string>
#include <unordered_map>
#include <tuple>

#include <cppunit/extensions/HelperMacros.h>
#include "../semirings/semiring.h"
#include "../semirings/free-semiring.h"

#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/var_degree_map.h"

#include "../utils/deriv_generator.h"
#include "../utils/string_util.h"

#include "commutative_monomial.h"

template <typename SR> using MonomialMap = std::map<CommutativeMonomial, SR> ;

template <typename SR>
class CommutativePolynomial : public Semiring<CommutativePolynomial<SR>,
                                   Commutativity::Commutative,
                                   SR::GetIdempotence()> {
  private:

     // FIXME: as a necessary intermediate step
     // (when calling make_commutative() for non-comm. polys we can have polys over non-comm. SR)
     //static_assert(SR::IsCommutative(),
     //    "The semiring SR must be commutative to be used with Polynomial!");


    /* Invariant:  The map is never empty.  In particular the 0 element of
     * Polynomial is represented as singleton map with empty monomial pointing
     * to the 0 element of the semiring. */

     MonomialMap<SR> monomials_;

    /* The maximum degree of each of the variables that appear in the
     * polynomial. */
    VarDegreeMap variables_;

    template <typename SR2>
    friend class CommutativePolynomial;

    static void InsertMonomial(MonomialMap<SR> &map, const CommutativeMonomial &m,
        const SR &c) {
      if (c == SR::null()) {
        return;
      }
      // FIXME: GCC 4.7 is missing emplace
      // map.emplace(m, c);
      map.insert(std::make_pair(m, c));
    }

    void InsertMonomial(MonomialMap<SR> &map, CommutativeMonomial &&m, SR &&c) {
      if (c == SR::null()) {
        return;
      }
      // FIXME: GCC 4.7 is missing emplace
      // map.emplace(std::move(m), std::move(c));
      map.insert(std::make_pair(std::move(m), std::move(c)));
    }

    void InsertMonomial(const CommutativeMonomial &m, const SR &c) {
      InsertMonomial(monomials_, m, c);
    }

    void InsertMonomial(CommutativeMonomial &&m, SR &&c) {
      InsertMonomial(monomials_, std::move(m), std::move(c));
    }

    CommutativePolynomial(MonomialMap<SR> &&ms, VarDegreeMap &&vs)
        : monomials_(std::move(ms)), variables_(std::move(vs)) {
      assert(SanityCheck());
    }

    /* Private constructor to hide the internal data structure. */
    CommutativePolynomial(const MonomialMap<SR> &ms) : monomials_(ms) {
      for (auto &monomial_coeff : monomials_) {
        variables_.Merge(monomial_coeff.first.variables_);
      }
      assert(SanityCheck());
    }

    bool SanityCheck() const {
      if (monomials_.empty()) {
        return variables_.empty();
      }
      VarDegreeMap tmp_variables;
      for (const auto &monomial_coeff : monomials_) {
        tmp_variables.Merge(monomial_coeff.first.variables_);
      }
      return tmp_variables == variables_;
    }

  public:
    CommutativePolynomial() = default;

    CommutativePolynomial(SR &&c, CommutativeMonomial &&m) {
      variables_ = m.variables_;
      InsertMonomial(std::move(m), std::move(c));

      assert(SanityCheck());
    }

    CommutativePolynomial(std::initializer_list< std::pair<SR, CommutativeMonomial> > init_list) {
      for (const auto &coeff_monomial : init_list) {
        auto iter = monomials_.find(coeff_monomial.second);
        if (iter == monomials_.end()) {
          variables_.Merge(coeff_monomial.second.variables_);
          InsertMonomial(coeff_monomial.second, coeff_monomial.first);
        } else {
          /* FIXME: can we use here std::move? */
          auto tmp_monomial_coeff = *iter;
          monomials_.erase(iter);
          InsertMonomial(tmp_monomial_coeff.first,
                         tmp_monomial_coeff.second + coeff_monomial.first);
        }
      }
      assert(SanityCheck());
    }

    CommutativePolynomial(const CommutativePolynomial &p) = default;
    CommutativePolynomial(CommutativePolynomial &&p) = default;

    /* Create a 'constant' polynomial. */
    CommutativePolynomial(const SR &elem) { InsertMonomial(CommutativeMonomial{}, elem); }
    CommutativePolynomial(SR &&elem) { InsertMonomial(CommutativeMonomial{}, std::move(elem)); }

    /* Create a polynomial which consists only of one variable. */
    CommutativePolynomial(const VarId var) {
      InsertMonomial(CommutativeMonomial{var}, SR::one());
      variables_.Insert(var);
      assert(SanityCheck());
    }

    CommutativePolynomial<SR>& operator=(const CommutativePolynomial<SR> &p) = default;
    CommutativePolynomial<SR>& operator=(CommutativePolynomial<SR> &&p) = default;

    CommutativePolynomial<SR>& operator+=(const CommutativePolynomial<SR> &polynomial) {
      for (const auto &monomial_coeff : polynomial.monomials_) {
        auto iter = monomials_.find(monomial_coeff.first);
        if (iter == monomials_.end()) {
          variables_.Merge(monomial_coeff.first.variables_);
          InsertMonomial(monomial_coeff.first, monomial_coeff.second);
        } else {
          /* FIXME: can we use here std::move? */
          auto tmp_monomial_coeff = *iter;
          monomials_.erase(iter);
          InsertMonomial(tmp_monomial_coeff.first,
                         tmp_monomial_coeff.second + monomial_coeff.second);
        }
      }

      assert(SanityCheck());
      return *this;
    }

    CommutativePolynomial<SR>& operator*=(const CommutativePolynomial<SR> &rhs) {
      if (monomials_.empty()) {
        return *this;
      } else if (rhs.monomials_.empty()) {
        monomials_.clear();
        variables_.clear();
        return *this;
      }

      MonomialMap<SR> tmp_monomials;
      VarDegreeMap tmp_variables;

      // iterate over both this and the poly polynomial
      for (const auto &lhs_monomial_coeff : monomials_) {
        for (const auto &rhs_monomial_coeff : rhs.monomials_) {
          auto tmp_monomial = lhs_monomial_coeff.first * rhs_monomial_coeff.first;
          auto tmp_coeff = lhs_monomial_coeff.second * rhs_monomial_coeff.second;

          auto iter = tmp_monomials.find(tmp_monomial);
          if (iter == tmp_monomials.end()) {
            /* The monomial is not in the list: update the variables and insert
             * the monomial. */
            tmp_variables.Merge(tmp_monomial.variables_);
            InsertMonomial(tmp_monomials, std::move(tmp_monomial),
                           std::move(tmp_coeff));
          } else {
            /* The monomial is already in the list, so add the coefficients. */
            iter->second = tmp_coeff + iter->second;
          }
        }
      }

      monomials_ = std::move(tmp_monomials);
      variables_ = std::move(tmp_variables);
      assert(SanityCheck());
      return *this;
    }

    // multiplying a polynomial with a variable
    CommutativePolynomial<SR> operator*(const VarId &var) const {
      MonomialMap<SR> tmp_monomials;
      VarDegreeMap tmp_variables;
      for (const auto &monomial_coeff : monomials_) {
        auto tmp_monomial = monomial_coeff.first * var;
        tmp_variables.Merge(tmp_monomial.variables_);
        InsertMonomial(tmp_monomials, std::move(tmp_monomial),
                                      monomial_coeff.second);
      }

      return CommutativePolynomial{std::move(tmp_monomials)};
    }

    friend CommutativePolynomial<SR> operator*(const SR &elem,
        const CommutativePolynomial<SR> &polynomial) {

      MonomialMap<SR> tmp_monomials;

      for (const auto &monomial_coeff : polynomial.monomials_) {
        InsertMonomial(tmp_monomials, monomial_coeff.first,
                       elem * monomial_coeff.second);
      }

      return CommutativePolynomial{std::move(tmp_monomials)};
    }

    bool operator==(const CommutativePolynomial<SR> &polynomial) const {

      return variables_ == polynomial.variables_ &&
             monomials_ == polynomial.monomials_;

    }

    CommutativePolynomial<SR> derivative(const VarId &var) const {
      return derivative_binom(VarDegreeMap{std::make_pair(var,1)});
    }

    CommutativePolynomial<SR> derivative_binom(const VarDegreeMap &vars) const {
      MonomialMap<SR> tmp_monomials;
      VarDegreeMap tmp_variables;

      for (const auto &monomial_coeff : monomials_) {
        /* Take the derivative of every monomial and add it to the result. */
        auto count_derivative = monomial_coeff.first.derivative_binom(vars);
        auto iter = tmp_monomials.find(count_derivative.second);
        SR tmp_coeff = monomial_coeff.second;
        tmp_coeff *= count_derivative.first;
        if (iter == tmp_monomials.end()) {
          tmp_variables.Merge(count_derivative.second.variables_);
          InsertMonomial(tmp_monomials, count_derivative.second,
                         std::move(tmp_coeff));
        } else {
          iter->second += std::move(tmp_coeff);
        }
      }
      return CommutativePolynomial{std::move(tmp_monomials), std::move(tmp_variables)};
    }

    SR AllNewtonDerivatives(const ValuationMap<SR> &previous_newton,
                            const ValuationMap<SR> &newton_update) const;


    SR DerivativeBinomAt(const std::unordered_map<VarId, Degree> &deriv_variables,
                         const ValuationMap<SR> &valuation) const;

    /*
     * Unfold each monomial and sum up everything
     * returns mappings from variables X_i to X_i^{<}
     */
    std::tuple<CommutativePolynomial<SR>,
               std::unordered_map<VarId, VarId>,
               std::unordered_map<VarId, VarId> > HeightUnfolding() const {

      std::unordered_map<VarId,VarId> prev_var_map;
      std::unordered_map<VarId,VarId> var_map;

      std::set<VarId> vars = get_variables();
      for (VarId v : vars) {
        VarId tmp = Var::GetVarId(); // get a fresh variable
        var_map.insert({v,tmp});

        VarId tmp2 = Var::GetVarId();
        prev_var_map.insert({v,tmp2});
      }

      //std::cout << "(poly) X^{<h}: "<< prev_var_map << std::endl;
      //std::cout << "(poly) X^{<h+1}: "<< var_map << std::endl;

      CommutativePolynomial<SR> result = CommutativePolynomial<SR>::null();

      for (const auto &monomial_coeff : monomials_) {
        result += monomial_coeff.first.HeightUnfolding(monomial_coeff.second, prev_var_map, var_map);
      }
      return std::make_tuple(result, prev_var_map, var_map);
    }

    static Matrix< CommutativePolynomial<SR> > jacobian(
        const std::vector< CommutativePolynomial<SR> > &polynomials,
        const std::vector<VarId> &variables) {
    	Matrix< CommutativePolynomial<SR> > result(polynomials.size(),variables.size());

    	std::unordered_map<VarId, unsigned int> varpos;
    	for(unsigned int j=0; j<variables.size(); ++j){
    		varpos[variables[j]] = j;
    	}

    	for (unsigned int i=0; i<polynomials.size(); ++i) {
    		for(const auto &monomial_coeff : polynomials[i].monomials_) {
    			for(const auto &v : monomial_coeff.first.variables_){
    				//variables_ is a var_degree_map
    				int j = varpos[v.first];
    				auto mult_mon = monomial_coeff.first.derivative(v.first);
    				result.At(i,j) += CommutativePolynomial(monomial_coeff.second*mult_mon.first, std::move(mult_mon.second));
    			}
    		}
    	}
    	return result;
    };

    SR eval(const ValuationMap<SR> &values) const {
      SR result = SR::null();
      for (const auto &monomial_coeff : monomials_) {
        result += monomial_coeff.second * monomial_coeff.first.eval(values);
      }
      return result;
    }

    static constexpr Commutativity GetCommutativity() { return Commutativity::Commutative; }

    /* Evaluate as much as possible (depending on the provided values) and
     * return both the concrete result and the remaining (i.e., unevaluated)
     * monomial. */
    CommutativePolynomial<SR> partial_eval(const ValuationMap<SR> &values) const {
      CommutativePolynomial<SR> result = CommutativePolynomial<SR>::null();
      for (const auto &monomial_coeff : monomials_) {
        auto tmp_coeff_monomial = monomial_coeff.first.partial_eval(values);
        result += CommutativePolynomial(monomial_coeff.second * tmp_coeff_monomial.first,
                             std::move(tmp_coeff_monomial.second));
      }
      return result;
    }

    /* Variable substitution. */
    CommutativePolynomial<SR> subst(const std::unordered_map<VarId, VarId> &mapping) const {
      MonomialMap<SR> tmp_monomials;

      for (const auto &monomial_coeff : monomials_) {
        InsertMonomial(tmp_monomials, monomial_coeff.first.subst(mapping),
                       monomial_coeff.second);
      }

      return CommutativePolynomial<SR>{std::move(tmp_monomials)};
    }

    static Matrix<SR> eval(const Matrix< CommutativePolynomial<SR> > &poly_matrix,
        const ValuationMap<SR> &values) {
      const std::vector< CommutativePolynomial<SR> > &tmp_polynomials = poly_matrix.getElements();
      std::vector<SR> result;
      for (const auto &polynomial : tmp_polynomials) {
        result.emplace_back(polynomial.eval(values));
      }
      return Matrix<SR>{poly_matrix.getRows(), std::move(result)};
    }

    static Matrix<CommutativePolynomial<SR> > eval(Matrix<CommutativePolynomial<SR> > poly_matrix,
        const std::unordered_map<VarId,CommutativePolynomial<SR> >& values) {
      const std::vector< CommutativePolynomial<SR> > &tmp_polynomials = poly_matrix.getElements();
      std::vector< CommutativePolynomial<SR> > result;
      for (const auto &polynomial : tmp_polynomials) {
        result.emplace_back(polynomial.eval(values));
      }
      return Matrix< CommutativePolynomial<SR> >{poly_matrix.getRows(), result};
    }

    /* Convert this polynomial to an element of the free semiring.  Note that
     * the provided valuation might be modified with new elements. */
    FreeSemiring make_free(std::unordered_map<SR, VarId, SR> *valuation) const {
      assert(valuation);

      auto result = FreeSemiring::null();
      // convert this polynomial by adding all converted monomials
      for (const auto &monomial_coeff : monomials_) {
        if (monomial_coeff.second == SR::null()) {
          assert(false); //coefficients in the monomial are always != 0.. so this should not happen :)
        } else if (monomial_coeff.second == SR::one()) {
         result += monomial_coeff.first.make_free();
        }
      else {
          auto value_iter = valuation->find(monomial_coeff.second);
          if (value_iter == valuation->end()) {
            /* Use a fresh constant - the constructor of Var::getVar() will take
             * care of this. */
            VarId tmp_var = Var::GetVarId();
            FreeSemiring tmp_var_free{tmp_var};
            valuation->emplace(monomial_coeff.second, tmp_var);
            result += tmp_var_free * monomial_coeff.first.make_free();
          } else {
            result += value_iter->second * monomial_coeff.first.make_free();
          }
        }
      }

      return result;
    }

    /* Same as make_free but for matrix form. */
    static Matrix<FreeSemiring> make_free(
        const Matrix< CommutativePolynomial<SR> > &poly_matrix,
        std::unordered_map<SR, VarId, SR> *valuation) {

      assert(valuation);

      const std::vector< CommutativePolynomial<SR> > &tmp_polynomials = poly_matrix.getElements();
      std::vector<FreeSemiring> result;

      for (const auto &polynomial : tmp_polynomials) {
        result.emplace_back(polynomial.make_free(valuation));
      }
      return Matrix<FreeSemiring>{poly_matrix.getRows(), std::move(result)};
    }

    // FIXME: this is inefficient!
    Degree get_degree() const {
      Degree degree = 0;
      for (auto &monomial_coeff : monomials_) {
        degree = std::max(degree, monomial_coeff.first.get_degree());
      }
      return degree;
    }

    Degree GetMaxDegreeOf(const VarId var) const {
      auto var_degree_iter = variables_.find(var);
      if (var_degree_iter == variables_.end()) {
        return 0;
      }
      return var_degree_iter->second;
    }

    const VarDegreeMap& GetVarDegreeMap() const {
      return variables_;
    }

    std::set<VarId> get_variables() const {
      std::set<VarId> vars;
      for (auto var_degree : variables_) {
        vars.insert(var_degree.first);
      }
      return vars;
    }

    static CommutativePolynomial<SR> null() {
      return CommutativePolynomial<SR>{};
    }

    static CommutativePolynomial<SR> one() {
      return CommutativePolynomial<SR>{SR::one()};
    }

    static bool is_idempotent;
    static bool is_commutative;

    std::string string() const {
      std::stringstream ss;
      for (auto iter = monomials_.begin(); iter != monomials_.end(); ++iter) {
        if (iter != monomials_.begin())
          ss << " + ";
        ss << iter->second << " * " << iter->first;
      }
    /*  ss << " degree info: ";
      for (auto &var_degree : variables_) {
        ss << var_degree.first << " |-> " << var_degree.second;
      }
*/
      return ss.str();
    }

    template <typename F>
    auto Map(F fun) const -> CommutativePolynomial<typename std::result_of<F(SR)>::type> {
      typedef typename std::result_of<F(SR)>::type SR2;

      /* Variables don't change, so just copy them over. */
      VarDegreeMap result_variables = variables_;
      MonomialMap<SR2> result_monomials;

      std::transform(
        monomials_.begin(), monomials_.end(),
        std::inserter(result_monomials, result_monomials.begin()),
        [&fun](const std::pair< CommutativeMonomial, SR > &pair) {
          return std::make_pair(pair.first, fun(pair.second));
        });
      return CommutativePolynomial<SR2>{std::move(result_monomials),
                             std::move(result_variables)};
    }
};



template <typename SR>
SR CommutativePolynomial<SR>::AllNewtonDerivatives(
    const ValuationMap<SR> &previous_newton,
    const ValuationMap<SR> &newton_update) const {

  std::unordered_map<VarId, Degree> current_max_degree;

  SR result = SR::null();

  for (const auto &monomial_coeff : monomials_) {
    if (monomial_coeff.first.variables_.size() == 0) {
      /* No variables, so the derivative over anything will be zero. */
      continue;
    }

    SR monomial_value = monomial_coeff.second;
    assert(!(monomial_value == SR::null()));

    current_max_degree.clear();
    Degree monomial_degree = 0;
    for (auto &var_degree : monomial_coeff.first.variables_) {
      monomial_degree += var_degree.second;
      current_max_degree.insert(var_degree);
    }

    if (monomial_degree < 2) {
      /* We look at only 2nd derivatives and in this case the value of the
       * monomial will be zero. */
      continue;
    }

    Generator generator{current_max_degree, 2, monomial_degree};

    while (generator.NextCombination()) {
      const auto &deriv_variables = generator.GetMap();

      /*for(auto &kv : deriv_variables) {
    	  std::cout << "[" << kv.first << ":" <<kv.second<<"]";
      }
      std::cout << std::endl;
       */

      /* First consider the variables that are in deriv_variables. */
      for (const auto &deriv_variable_degree : deriv_variables) {
        const auto variable = deriv_variable_degree.first;
        const auto deriv_degree = deriv_variable_degree.second;

        if (deriv_degree == 0) {
          continue;
        }

        auto degree = monomial_coeff.first.GetDegreeOf(variable);

        if (deriv_degree > degree) {
          monomial_value = SR::null();
          break;
        }

        auto binomial_coeff_d =
          boost::math::binomial_coefficient<double>(degree, deriv_degree);
        /* Check if we don't overflow. */
        assert(static_cast<std::uint_fast64_t>(binomial_coeff_d) <=
               static_cast<std::uint_fast64_t>(
                    std::numeric_limits<Degree>::max()));
        monomial_value *= static_cast<Degree>(binomial_coeff_d);

        if (degree > deriv_degree) {
          auto value_lookup = previous_newton.find(variable);
          assert(value_lookup != previous_newton.end());
          for (Degree c = 0; c < degree - deriv_degree; ++c) {
            monomial_value *= value_lookup->second;
          }
        }
      }

      /* If the current monomial_value is 0 (and thus it'll remain to be 0),
       * we can continue with the next monomial. */
      if (monomial_value == SR::null()) {
        continue;
      }

      /* Finally consider all the variables that are *not* in deriv_variables. */
      for (const auto &variable_degree : monomial_coeff.first) {
        const auto variable = variable_degree.first;
        const auto degree = variable_degree.second;

        auto lookup = deriv_variables.find(variable);
        if (lookup != deriv_variables.end() && lookup->second > 0) {
          /* Already considered in the previous loop. */
          continue;
        }

        auto value_lookup = previous_newton.find(variable);
        assert(value_lookup != previous_newton.end());
        for (Degree c = 0; c < degree; ++c) {
          monomial_value *= value_lookup->second;
        }
      }

      SR prod = SR::one();

      for (const auto &var_val : generator.GetMap()) {
        for (unsigned int j = 0; j < var_val.second; ++j) {
          auto lookup = newton_update.find(var_val.first);
          assert(lookup != newton_update.end());
          prod *= lookup->second;
        }
      }

      result += monomial_value * prod;
    }
  }

  // DMSG("DerivativeBinomAt:");
  // DMSG(result);
  // DMSG("eval . derivative_binom:");
  // DMSG(derivative_binom(deriv_variables).eval(valuation));

  return result;
  // for every (coefficient, monomial):
  //   for every (variable, degree) in monomial:
  //     find the corresponding deriv_degree in deriv_variables
  //     if (devir_degree >= degree)
  //       the derivative is equal to zero, so continue with the next
  //       monomial
  //     otherwise
  //       tmp = coefficient * BinomCoeff(deegre, deriv_degree)
  //     if it doesn't exist return lookup the variable in valuation
}


template <typename SR>
SR CommutativePolynomial<SR>::DerivativeBinomAt(
    const std::unordered_map<VarId, Degree> &deriv_variables,
    const ValuationMap<SR> &valuation) const {

  SR result = SR::null();
  for (const auto &monomial_coeff : monomials_) {
    SR monomial_value = monomial_coeff.second;
    assert(!(monomial_value == SR::null()));

    /* First consider the variables that are in deriv_variables. */
    for (const auto &deriv_variable_degree : deriv_variables) {
      const auto variable = deriv_variable_degree.first;
      const auto deriv_degree = deriv_variable_degree.second;

      if (deriv_degree == 0) {
        continue;
      }

      auto degree = monomial_coeff.first.GetDegreeOf(variable);

      if (deriv_degree > degree) {
        monomial_value = SR::null();
        break;
      }

      auto binomial_coeff_d =
        boost::math::binomial_coefficient<double>(degree, deriv_degree);
      /* Check if we don't overflow. */
      assert(static_cast<std::uint_fast64_t>(binomial_coeff_d) <=
             static_cast<std::uint_fast64_t>(
                  std::numeric_limits<Degree>::max()));
      monomial_value *= static_cast<Degree>(binomial_coeff_d);

      if (degree > deriv_degree) {
        auto value_lookup = valuation.find(variable);
        assert(value_lookup != valuation.end());
        for (Degree c = 0; c < degree - deriv_degree; ++c) {
          monomial_value *= value_lookup->second;
        }
      }
    }

    /* If the current monomial_value is 0 (and thus it'll remain to be 0),
     * we can continue with the next monomial. */
    if (monomial_value == SR::null()) {
      continue;
    }

    /* Finally consider all the variables that are *not* in deriv_variables. */
    for (const auto &variable_degree : monomial_coeff.first) {
      const auto variable = variable_degree.first;
      const auto degree = variable_degree.second;

      auto lookup = deriv_variables.find(variable);
      if (lookup != deriv_variables.end() && lookup->second > 0) {
        /* Already considered in the previous loop. */
        continue;
      }

      auto value_lookup = valuation.find(variable);
      assert(value_lookup != valuation.end());
      for (Degree c = 0; c < degree; ++c) {
        monomial_value *= value_lookup->second;
      }
    }
    result += monomial_value;
  }

  // DMSG("DerivativeBinomAt:");
  // DMSG(result);
  // DMSG("eval . derivative_binom:");
  // DMSG(derivative_binom(deriv_variables).eval(valuation));

  return result;
  // for every (coefficient, monomial):
  //   for every (variable, degree) in monomial:
  //     find the corresponding deriv_degree in deriv_variables
  //     if (devir_degree >= degree)
  //       the derivative is equal to zero, so continue with the next
  //       monomial
  //     otherwise
  //       tmp = coefficient * BinomCoeff(deegre, deriv_degree)
  //     if it doesn't exist return lookup the variable in valuation
}



template <typename SR>
std::ostream& operator<<(std::ostream& os,
    const std::unordered_map<VarId, CommutativePolynomial<SR> >& values) {
  for (const auto &key_value : values) {
    os << key_value->first << "→" << key_value->second << ";";
  }
  return os;
}


#endif
