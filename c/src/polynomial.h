#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <list>
#include <map>
#include <string>
#include <unordered_map>

#include "free-semiring.h"
#include "matrix.h"
#include "monomial.h"
#include "semiring.h"
#include "var.h"
#include "var_degree_map.h"


/*
 * TODO:
 * - Add SanityCheck that checks all the invariants.
 */


/* FIXME: Polynomials are no semiring in our definition (not starable). */

template <typename SR>
class Polynomial : public Semiring< Polynomial<SR> > {
  private:
    /* Invariant:  The map is never empty.  In particular the 0 element of
     * Polynomial is represented as singleton map with empty monomial pointing
     * to the 0 element of the semiring. */
    std::map<Monomial, SR> monomials_;

    /* The maximum degree of each of the variables that appear in the
     * polynomial. */
    VarDegreeMap variables_;

    static void InsertMonomial(std::map<Monomial, SR> &map, const Monomial &m,
        const SR &c) {
      if (c == SR::null()) {
        return;
      }
      // FIXME: GCC 4.7 is missing emplace
      // map.emplace(m, c);
      map.insert(std::make_pair(m, c));
    }

    void InsertMonomial(std::map<Monomial, SR> &map, Monomial &&m, SR &&c) {
      if (c == SR::null()) {
        return;
      }
      // FIXME: GCC 4.7 is missing emplace
      // map.emplace(std::move(m), std::move(c));
      map.insert(std::make_pair(std::move(m), std::move(c)));
    }

    void InsertMonomial(const Monomial &m, const SR &c) {
      InsertMonomial(monomials_, m, c);
    }

    void InsertMonomial(Monomial &&m, SR &&c) {
      InsertMonomial(monomials_, std::move(m), std::move(c));
    }

    Polynomial(std::map<Monomial, SR> &&ms, VarDegreeMap &&vs)
        : monomials_(std::move(ms)), variables_(std::move(vs)) {
      assert(SanityCheck());
    }

    /* Private constructor to hide the internal data structure. */
    Polynomial(const std::map<Monomial, SR> &ms) : monomials_(ms) {
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
    Polynomial() = default;

    Polynomial(SR &&c, Monomial &&m) {
      variables_ = m.variables_;
      InsertMonomial(std::move(m), std::move(c));

      assert(SanityCheck());
    }

    Polynomial(std::initializer_list< std::pair<SR, Monomial> > init_list) {
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

    Polynomial(const Polynomial &p) = default;
    Polynomial(Polynomial &&p) = default;

    /* Create a 'constant' polynomial. */
    Polynomial(const SR &elem) { InsertMonomial(Monomial{}, elem); }
    Polynomial(SR &&elem) { InsertMonomial(Monomial{}, std::move(elem)); }

    /* Create a polynomial which consists only of one variable. */
    Polynomial(const VarPtr var) {
      InsertMonomial(Monomial{var}, SR::one());
      variables_.Insert(var);
      assert(SanityCheck());
    }

    Polynomial<SR>& operator=(const Polynomial<SR> &p) = default;
    Polynomial<SR>& operator=(Polynomial<SR> &&p) = default;

    Polynomial<SR>& operator+=(const Polynomial<SR> &polynomial) {
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

    Polynomial<SR>& operator*=(const Polynomial<SR> &rhs) {
      if (monomials_.empty()) {
        return *this;
      } else if (rhs.monomials_.empty()) {
        monomials_.clear();
        variables_.clear();
        return *this;
      }

      std::map<Monomial, SR> tmp_monomials;
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
    Polynomial<SR> operator*(const VarPtr &var) const {
      std::map<Monomial, SR> tmp_monomials;
      VarDegreeMap tmp_variables;
      for (const auto &monomial_coeff : monomials_) {
        auto tmp_monomial = monomial_coeff.first * var;
        tmp_variables.Merge(tmp_monomial.variables_);
        InsertMonomial(tmp_monomials, std::move(tmp_monomial),
                                      monomial_coeff.second);
      }

      return Polynomial{std::move(tmp_monomials)};
    }

    friend Polynomial<SR> operator*(const SR &elem,
        const Polynomial<SR> &polynomial) {

      std::map<Monomial, SR> tmp_monomials;

      for (const auto &monomial_coeff : polynomial.monomials_) {
        InsertMonomial(tmp_monomials, monomial_coeff.first,
                       elem * monomial_coeff.second);
      }

      return Polynomial{std::move(tmp_monomials)};
    }

    bool operator==(const Polynomial<SR> &polynomial) const {

      return variables_ == polynomial.variables_ &&
             monomials_ == polynomial.monomials_;

    }

    Polynomial<SR> derivative(const VarPtr& var) const {
      std::map<Monomial, SR> tmp_monomials;
      VarDegreeMap tmp_variables;

      for (const auto &monomial_coeff : monomials_) {
        /* Take the derivative of every monomial and add it to the result. */
        auto count_derivative = monomial_coeff.first.derivative(var);
        auto iter = tmp_monomials.find(count_derivative.second);
        SR tmp_coeff = SR::null();
        for (Degree i = 0; i < count_derivative.first; ++i) {
          tmp_coeff += monomial_coeff.second;
        }
        if (iter == tmp_monomials.end()) {
          tmp_variables.Merge(count_derivative.second.variables_);
          InsertMonomial(tmp_monomials, count_derivative.second,
                         std::move(tmp_coeff));
        } else {
          iter->second += std::move(tmp_coeff);
        }
      }

      return Polynomial{std::move(tmp_monomials), std::move(tmp_variables)};
    }

    Polynomial<SR> derivative(const std::vector<VarPtr> &vars) const {
      Polynomial<SR> tmp_polynomial = *this;
      for (auto var : vars) {
        tmp_polynomial = tmp_polynomial.derivative(var);
      }
      return tmp_polynomial;
    }

    static Matrix< Polynomial<SR> > jacobian(
        const std::vector< Polynomial<SR> > &polynomials,
        const std::vector<VarPtr> &variables) {
      std::vector< Polynomial<SR> > result_vector;
      for (const auto &polynomial : polynomials) {
        for (const auto variable : variables) {
          result_vector.push_back(polynomial.derivative(variable));
        }
      }
      // FIXME: Clean up Matrix and then remove the casts...
      return Matrix< Polynomial<SR> >{ static_cast<int>(variables.size()),
                                       static_cast<int>(polynomials.size()),
                                       std::move(result_vector) };
    };

    SR eval(const std::map<VarPtr, SR> &values) const {
      SR result = SR::null();
      for (const auto &monomial_coeff : monomials_) {
        result += monomial_coeff.second * monomial_coeff.first.eval(values);
      }
      return result;
    }

    /* Evaluate as much as possible (depending on the provided values) and
     * return both the concrete result and the remaining (i.e., unevaluated)
     * monomial. */
    Polynomial<SR> partial_eval(const std::map<VarPtr, SR> &values) const {
      Polynomial<SR> result = Polynomial<SR>::null();
      for (const auto &monomial_coeff : monomials_) {
        auto tmp_coeff_monomial = monomial_coeff.first.partial_eval(values);
        result += Polynomial(monomial_coeff.second * tmp_coeff_monomial.first,
                             std::move(tmp_coeff_monomial.second));
      }
      return result;
    }

    /* Variable substitution. */
    Polynomial<SR> subst(const std::map<VarPtr, VarPtr> &mapping) const {
      std::map<Monomial, SR> tmp_monomials;

      for (const auto &monomial_coeff : monomials_) {
        InsertMonomial(tmp_monomials, monomial_coeff.first.subst(mapping),
                       monomial_coeff.second);
      }

      return Polynomial<SR>{std::move(tmp_monomials)};
    }

    static Matrix<SR> eval(const Matrix< Polynomial<SR> > &poly_matrix,
        const std::map<VarPtr, SR> &values) {
      const std::vector< Polynomial<SR> > &tmp_polynomials = poly_matrix.getElements();
      std::vector<SR> result;
      for (const auto &polynomial : tmp_polynomials) {
        result.emplace_back(polynomial.eval(values));
      }
      return Matrix<SR>{poly_matrix.getColumns(), poly_matrix.getRows(),
                        std::move(result)};
    }

    // FIXME: Why values is passed by value???
    static Matrix<Polynomial<SR> > eval(Matrix<Polynomial<SR> > poly_matrix,
        std::map<VarPtr,Polynomial<SR> > values) {
      const std::vector< Polynomial<SR> > &tmp_polynomials = poly_matrix.getElements();
      std::vector< Polynomial<SR> > result;
      for (const auto &polynomial : tmp_polynomials) {
        result.emplace_back(polynomial.eval(values));
      }
      return Matrix< Polynomial<SR> >{poly_matrix.getColumns(),
                                      poly_matrix.getRows(), result};
    }

    /* Convert this polynomial to an element of the free semiring.  Note that
     * the provided valuation might be modified with new elements. */
    FreeSemiring make_free(std::unordered_map<SR, VarPtr, SR> *valuation) const {
      assert(valuation);

      /* FIXME: Do we need that?  The above assertion should let us know...
      if (!valuation) {
        valuation = new std::unordered_map<SR, VarPtr, SR>();
      }
      */

      auto result = FreeSemiring::null();
      // convert this polynomial by adding all converted monomials
      for (const auto &monomial_coeff : monomials_) {

        if (monomial_coeff.second == SR::null()) {
          result += FreeSemiring::null();
        } else if (monomial_coeff.second == SR::one()) {
          result += FreeSemiring::one();
        } else {
          auto value_iter = valuation->find(monomial_coeff.second);
          if (value_iter == valuation->end()) {
            /* Use a fresh constant - the constructor of Var::getVar() will take
             * care of this. */
            VarPtr tmp_var = Var::getVar();
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
        const Matrix< Polynomial<SR> > &poly_matrix,
        std::unordered_map<SR, VarPtr, SR> *valuation) {

      assert(valuation);
      // FIXME: Again, do we need that?
      // if (!valuation)
      //   valuation = new std::unordered_map<SR, VarPtr, SR>();


      const std::vector< Polynomial<SR> > &tmp_polynomials = poly_matrix.getElements();
      std::vector<FreeSemiring> result;

      for (const auto &polynomial : tmp_polynomials) {
        result.emplace_back(polynomial.make_free(valuation));
      }
      return Matrix<FreeSemiring>{poly_matrix.getColumns(),
                                  poly_matrix.getRows(), std::move(result)};
    }

    Degree get_degree() {
      Degree degree = 0;
      for (auto &monomial_coeff : monomials_) {
        degree = std::max(degree, monomial_coeff.first.get_degree());
      }
      return degree;
    }

    Degree GetMaxDegreeOf(const VarPtr var) {
      auto var_degree_iter = variables_.find(var);
      if (var_degree_iter == variables_.end()) {
        return 0;
      }
      return var_degree_iter->second;
    }

    /* FIXME: Get rid of this. */
    std::set<VarPtr> get_variables() const {
      std::set<VarPtr> vars;
      for (auto var_degree : variables_) {
        vars.insert(var_degree.first);
      }
      return vars;
    }

    // TODO: we cannot star polynomials!
    Polynomial<SR> star() const {
      assert(false);
      return (*this);
    }

    static Polynomial<SR> null() {
      return Polynomial<SR>{};
    }

    static Polynomial<SR> one() {
      return Polynomial<SR>{SR::one()};
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
      ss << " degree info: ";
      for (auto &var_degree : variables_) {
        ss << var_degree.first << " |-> " << var_degree.second;
      }

      return ss.str();
    }
};

template <typename SR> bool Polynomial<SR>::is_commutative = false;
template <typename SR> bool Polynomial<SR>::is_idempotent = false;

template <typename SR>
std::ostream& operator<<(std::ostream& os, const std::map<VarPtr, SR>& values) {
  for (auto value = values.begin(); value != values.end(); ++value) {
    os << value->first << "→" << value->second << ";";
  }
  return os;
}


template <typename SR>
std::ostream& operator<<(std::ostream& os,
    const std::map<VarPtr, Polynomial<SR> >& values) {
  for (const auto &key_value : values) {
    os << key_value->first << "→" << key_value->second << ";";
  }
  return os;
}


#endif
