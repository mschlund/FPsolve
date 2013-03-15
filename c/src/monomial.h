#pragma once

#include <iosfwd>

#include "free-semiring.h"
#include "var_degree_map.h"


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

    /* Evaluate the monomial given the map from variables to values. */
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

    Degree get_degree() const {
      Degree degree = 0;
      for (auto var_degree : variables_) {
        degree += var_degree.second;
      }
      return degree;
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
