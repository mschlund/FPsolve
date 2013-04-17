#pragma once

#include <iosfwd>

#include "../semirings/free-semiring.h"
#include "../datastructs/var_degree_map.h"

#include <boost/math/special_functions/binomial.hpp>


template <typename SR>
class Polynomial;

class Monomial {
  private:
    /* Maps each variable to its degree. */
    VarDegreeMap variables_;

    template <typename SR>
    friend class Polynomial;

    /* Private constructor to not leak the internal data structure. */
    Monomial(VarDegreeMap &&vs) : variables_(std::move(vs)) {}

  public:
    typedef typename VarDegreeMap::iterator iterator;
    typedef typename VarDegreeMap::const_iterator const_iterator;


    /* Since we don't manage any resources, we can simply use the default
     * constructors and assignment operators. */
    Monomial() = default;
    Monomial(const Monomial &m) = default;
    Monomial(Monomial &&m) = default;
    Monomial& operator=(const Monomial &m) = default;
    Monomial& operator=(Monomial &&m) = default;

    Monomial(std::initializer_list<VarId> vs) {
      for (auto var : vs) {
        variables_.Insert(var);
      }
    }

    /* std::vector seems to be a neutral data type and does not leak internal
     * data structure. */
    Monomial(std::vector< std::pair<VarId, Degree> > vs) {
      for (auto var_degree : vs) {
        variables_.Insert(var_degree.first, var_degree.second);
      }
    }

    /* Multiply two monomials. */
    Monomial operator*(const Monomial &monomial) const {
      auto tmp_variables = variables_;
      for (auto var_degree : monomial.variables_) {
        tmp_variables.Insert(var_degree.first, var_degree.second);
      }
      return Monomial(std::move(tmp_variables));
    }

    /* Multiply a monomial with a variable. */
    Monomial operator*(const VarId &var) const {
      auto tmp_variables = variables_;
      tmp_variables.Insert(var);
      return Monomial{std::move(tmp_variables)};
    }

    /* Commutative version of derivative. */
    std::pair<Degree, Monomial> derivative(const VarId &var) const {

      auto var_degree_iter = variables_.find(var);

      /* If the variable does not appear in the monomial, the derivative
       * must be 0. */
      if (var_degree_iter == variables_.end()) {
        return { 0, Monomial{} };
      }

      /* Remove one of these by removing the first of them and then "multiply"
       * the coefficient with degree_before. */
      auto tmp_variables = variables_;
      tmp_variables.Erase(var);

      return { var_degree_iter->second, Monomial{std::move(tmp_variables)} };
    }

    /*
     * take the derivative w.r.t. a "multiindex" of variables (e.g. d/dx^i = d/dx^2y^3 for i=(2,3))
     *  and divide by i! (= 2!*3! in the example),
     *  effectively this deletes as many variables as indicated by i and multiplies
     *  with binom{D}{i} where D is the vector of variable degrees
     *  i.e. if above operator d/dx^2y^3 is applied to the monomial x^5y^6 the result is (5*4 * 6*5*4) * x^3y^3
     */
    std::pair<Degree, Monomial> derivative_binom(const std::map<VarId, Degree> &vars) const {
      Degree multiplicity = 1;
      auto tmp_variables = variables_;

      for(const auto var : vars) {
        auto var_degree_iter = tmp_variables.find(var.first);
        /* If the variable does not appear in the monomial, the derivative
         * must be 0. */
        if (var_degree_iter == tmp_variables.end()) {
          return { 0, Monomial{} };
        }
        else if (var.second > var_degree_iter->second){
          /* If we derive more often wrt. some var than this variable is present, we get 0.
           */
          return { 0, Monomial{} };
        }
        /*
         * this would give the usual derivative:
         *  multiply the coefficient by (K!) where K is the number of times it is derived
         *  multiplicity *= (Degree) boost::math::factorial<float>(var_degree_iter->second);
         */

        // devide automatically by i! if we compute d/dx^i
        multiplicity *= (Degree) boost::math::binomial_coefficient<double>(var_degree_iter->second, var.second);

        //Reduce the multiplicity of the variable by the number of times K it is derived
        tmp_variables.Erase(var.first, var.second);
      }

      /* Remove one of these by removing the first of them and then "multiply"
       * the coefficient with degree_before. */
      return {multiplicity, Monomial{std::move(tmp_variables)} };
    }

    //TODO: avoid creation of temporary objects
    template <typename SR>
    SR derivative_binom_at(const std::map<VarId, Degree> &vars,
                           const std::map<VarId, SR> &valuation) const {




      auto tmp_deriv = derivative_binom(vars);
      SR res = tmp_deriv.second.eval(valuation);
      res *= tmp_deriv.first;
      return res;
    }

    /* Evaluate the monomial given the map from variables to values,
     * If a variable is not interpreted the assertion fails !*/
    template <typename SR>
    SR eval(const std::map<VarId, SR> &values) const {
      auto result = SR::one();

      for (auto var_degree : variables_) {
        auto value_iter = values.find(var_degree.first);
        /* All variables should be in the values map. */
        assert(value_iter != values.end());
        for (Degree i = 0; i < var_degree.second; ++i) {
          result *= value_iter->second;
        }
      }

      return result;
    }

    /* Partially evaluate the monomial. */
    template <typename SR>
    std::pair<SR, Monomial> partial_eval(
        const std::map<VarId, SR> &values) const {

      SR result_value = SR::one();
      Monomial result_monomial;

      for (auto var_degree : variables_) {
        auto value_iter = values.find(var_degree.first);
        if (value_iter == values.end()) {
          /* Variable not found in the mapping, so keep it. */
          result_monomial.variables_.Insert(var_degree.first, var_degree.second);
        } else {
          /* Variable found, use it for evaluation. */
          for (Degree i = 0; i < var_degree.second; ++i) {
            result_value *= value_iter->second;
          }
        }
      }

      return { std::move(result_value), std::move(result_monomial) };
    }

    /* Variable substitution. */
    Monomial subst(const std::map<VarId, VarId> &mapping) const {
      VarDegreeMap tmp_variables;

      for (const auto &var_degree : variables_) {
        auto old_new_iter = mapping.find(var_degree.first);
        if (old_new_iter != mapping.end()) {
          tmp_variables.Insert(old_new_iter->second, var_degree.second);
        } else {
          tmp_variables.Insert(var_degree.first, var_degree.second);
        }
      }

      return Monomial{std::move(tmp_variables)};
    }

    /* Convert this monomial to an element of the free semiring. */
    FreeSemiring make_free() const {
      FreeSemiring result = FreeSemiring::one();
      for (auto var_degree : variables_) {
        FreeSemiring tmp{var_degree.first};
        for (Degree i = 0; i < var_degree.second; ++i) {
          result = result * tmp;
        }
      }
      return result;
    }

    bool operator<(const Monomial &rhs) const {
      return variables_ < rhs.variables_;
    }

    bool operator==(const Monomial &rhs) const {
      return variables_ == rhs.variables_;
    }

    iterator begin() { return variables_.begin(); }
    iterator end() { return variables_.end(); }

    const_iterator begin() const { return variables_.begin(); }
    const_iterator end() const { return variables_.end(); }

    Degree get_degree() const {
      Degree degree = 0;
      for (auto var_degree : variables_) {
        degree += var_degree.second;
      }
      return degree;
    }

    Degree GetDegreeOf(const VarId var) const {
      return variables_.GetDegreeOf(var);
    }


    // FIXME: modify or remove
    std::set<VarId> get_variables() const {
      std::set<VarId> set;
      for (auto var_degree : variables_) {
        set.insert(var_degree.first);
      }
      return set;
    }

    std::string string() const {
      std::stringstream ss;
      ss << variables_;
      return std::move(ss.str());
    }
};

std::ostream& operator<<(std::ostream &out, const Monomial &monomial);
