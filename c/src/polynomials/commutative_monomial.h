#pragma once

#include <iosfwd>

//#include "commutative_polynomial.h"
#include "../semirings/free-semiring.h"
#include "../datastructs/var_degree_map.h"
#include "../utils/string_util.h"

#include <boost/math/special_functions/binomial.hpp>

template <typename SR>
class CommutativePolynomial;

class CommutativeMonomial {
  private:
    /* Maps each variable to its degree. */
    VarDegreeMap variables_;

    template <typename SR>
    friend class CommutativePolynomial;

    /* Private constructor to not leak the internal data structure. */
    CommutativeMonomial(VarDegreeMap &&vs) : variables_(std::move(vs)) {}

  public:
    typedef typename VarDegreeMap::iterator iterator;
    typedef typename VarDegreeMap::const_iterator const_iterator;


    /* Since we don't manage any resources, we can simply use the default
     * constructors and assignment operators. */
    CommutativeMonomial() = default;
    CommutativeMonomial(const CommutativeMonomial &m) = default;
    CommutativeMonomial(CommutativeMonomial &&m) = default;
    CommutativeMonomial& operator=(const CommutativeMonomial &m) = default;
    CommutativeMonomial& operator=(CommutativeMonomial &&m) = default;

    CommutativeMonomial(std::initializer_list<VarId> vs) {
      for (auto var : vs) {
        variables_.Insert(var);
      }
    }

    /* std::vector seems to be a neutral data type and does not leak internal
     * data structure. */
    CommutativeMonomial(std::vector< std::pair<VarId, Degree> > vs) {
      for (auto var_degree : vs) {
        variables_.Insert(var_degree.first, var_degree.second);
      }
    }

    /* Multiply two monomials. */
    CommutativeMonomial operator*(const CommutativeMonomial &monomial) const {
      auto tmp_variables = variables_;
      for (auto var_degree : monomial.variables_) {
        tmp_variables.Insert(var_degree.first, var_degree.second);
      }
      return CommutativeMonomial(std::move(tmp_variables));
    }

    /* Multiply a monomial with a variable. */
    CommutativeMonomial operator*(const VarId &var) const {
      auto tmp_variables = variables_;
      tmp_variables.Insert(var);
      return CommutativeMonomial{std::move(tmp_variables)};
    }

    /* Commutative version of derivative. */
    std::pair<Degree, CommutativeMonomial> derivative(const VarId &var) const {

      auto var_degree_iter = variables_.find(var);

      /* If the variable does not appear in the monomial, the derivative
       * must be 0. */
      if (var_degree_iter == variables_.end()) {
        return { 0, CommutativeMonomial{} };
      }

      /* Remove one of these by removing the first of them and then "multiply"
       * the coefficient with degree_before. */
      auto tmp_variables = variables_;
      tmp_variables.Erase(var);

      return { var_degree_iter->second, CommutativeMonomial{std::move(tmp_variables)} };
    }

    /*
     * take the derivative w.r.t. a "multiindex" of variables (e.g. d/dx^i = d/dx^2y^3 for i=(2,3))
     *  and divide by i! (= 2!*3! in the example),
     *  effectively this deletes as many variables as indicated by i and multiplies
     *  with binom{D}{i} where D is the vector of variable degrees
     *  i.e. if above operator d/dx^2y^3 is applied to the monomial x^5y^6 the result is (5*4 * 6*5*4) * x^3y^3
     */
    std::pair<Degree, CommutativeMonomial> derivative_binom(const VarDegreeMap &vars) const {
      Degree multiplicity = 1;
      auto tmp_variables = variables_;

      for(const auto var : vars) {
        auto var_degree_iter = tmp_variables.find(var.first);
        /* If the variable does not appear in the monomial, the derivative
         * must be 0. */
        if (var_degree_iter == tmp_variables.end()) {
          return { 0, CommutativeMonomial{} };
        }
        else if (var.second > var_degree_iter->second){
          /* If we derive more often wrt. some var than this variable is present, we get 0.
           */
          return { 0, CommutativeMonomial{} };
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
      return {multiplicity, CommutativeMonomial{std::move(tmp_variables)} };
    }

    /* Evaluate the monomial given the map from variables to values,
     * If a variable is not interpreted the assertion fails !*/
    template <typename SR>
    SR eval(const ValuationMap<SR> &values) const {
      auto result = SR::one();

      for (auto var_degree : variables_) {
        auto value_iter = values.find(var_degree.first);
        /* All variables should be in the values map. */
        assert(value_iter != values.end());
        // exponentiation is more efficient than iterated multiplication (binary exp.)
        result *= (value_iter->second ^ var_degree.second);

        /*for (Degree i = 0; i < var_degree.second; ++i) {
          result *= value_iter->second;
        }*/
      }

      return result;
    }


    /*
     * for a monomial \sum_i X_i^{d_i}, the height unfolding is:
     * sum_{i=1}^n
     * \prod_{k=1}^{i-1} (X_k^{<h+1})^d_k *
     * (\sum_{k=0}^{d_i -1} (X_i^{<h+1})^{d_i-k-1}) * (X_i^{=h}) * previous_values(X_i)^k )
     * \prod_{k=i+1}^{n} previous_values(X_k)^{d_k}
     *
     * Important to note: the valuation of X^{<h+1} is computed one round _before_ X^{=h+1}
     *
     * Explanation:
     * example rule X -> X Y Z, becomes X^{=h+1} -> X^{=h} Y^{<h} Z^{<h} + X^{<h+1} Y^{=h} Z^{<h} + X^{<h+1} Y^{<h+1} Z^{=h}
     * (in general we can have variables with exponents >1... this requires the ugly sum-expression in the second line above :))
     * To derive an X-tree of height =h+1, every variable on the rhs can derive a tree of height =h (that is the outer sum)
     * The variables "to the left" of this =h-variable generate trees of height <h,
     * the variables "to the right" of it generate trees of height <(h-1).
     * This way we enumerate every possibility to generate a tree of height =h+1 exactly once
     * (first we choose the "last" variable i to become the =h, then in this block of variables we have to sum over all
     * (d_i-1) partition-points k of the remaining variables into those of type "<h+1" and "<h".
     */


    /* Return the unfolding of the monomial (i.e. a polynomial) using the two variable-renaming-maps given
     * the first for the variables of type <h-1
     * the second one for variables of type <h, all remaining variables will be of type =h
     */
    template <typename SR>
    CommutativePolynomial<SR> HeightUnfolding(const SR& coeff, const std::unordered_map<VarId,VarId>& prev_var_map,
                                                          const std::unordered_map<VarId,VarId>& var_map)  const{
      //std::cout << "(mon) X^{<h}: "<< prev_var_map << std::endl;
      //std::cout << "(mon) X^{<h+1}: "<< var_map << std::endl;

      CommutativePolynomial<SR> result = CommutativePolynomial<SR>::null();

      // holds \prod_{k=1}^{i-1} (X_k^{<h+1})^d_k, initially empty, then add factor in every run of the loop
      VarDegreeMap left_prod = VarDegreeMap();

      // holds \prod_{k=i+1}^{n} previous_values(X_k)^{d_k},
      // initially prod of _all_ vars then remove variables in the loop
      VarDegreeMap right_prod = VarDegreeMap();
      for (auto var_degree : variables_) {
        right_prod.Insert(prev_var_map.at(var_degree.first), var_degree.second);
      }
      SR coeff_copy = coeff;

      for (auto var_degree : variables_) {
        right_prod.EraseAll(prev_var_map.at(var_degree.first));

        for (Degree k=0; k < var_degree.second; ++k) {
          VarDegreeMap center_prod = VarDegreeMap();

          center_prod.Insert(var_map.at(var_degree.first), k);
          center_prod.Insert(var_degree.first, 1);
          center_prod.Insert(prev_var_map.at(var_degree.first), var_degree.second - k -1);

          center_prod.Merge(left_prod);
          center_prod.Merge(right_prod);

          result += CommutativePolynomial<SR>(std::move(coeff_copy), std::move(CommutativeMonomial(std::move(center_prod))));
        }

        left_prod.Insert(var_map.at(var_degree.first), var_degree.second);

      }

      return result;
    }


    /* Partially evaluate the monomial. */
    template <typename SR>
    std::pair<SR, CommutativeMonomial> partial_eval(
        const ValuationMap<SR> &values) const {

      SR result_value = SR::one();
      CommutativeMonomial result_monomial;


      for (auto var_degree : variables_) {
        auto value_iter = values.find(var_degree.first);
        if (value_iter == values.end()) {
          /* Variable not found in the mapping, so keep it. */
          result_monomial.variables_.Insert(var_degree.first, var_degree.second);
        } else {
          /* Variable found, use it for evaluation. */
          result_value *= (value_iter->second ^ var_degree.second);
        }
      }

      // the coefficient is 0 if and only if the monomial is 0---so make sure that the variable-map is empty
      if(SR::null() == result_value) {
    	  result_monomial.variables_.clear();
      }

      return { std::move(result_value), std::move(result_monomial) };
    }

    /* Variable substitution. */
    CommutativeMonomial subst(const std::unordered_map<VarId, VarId> &mapping) const {
      VarDegreeMap tmp_variables;

      for (const auto &var_degree : variables_) {
        auto old_new_iter = mapping.find(var_degree.first);
        if (old_new_iter != mapping.end()) {
          tmp_variables.Insert(old_new_iter->second, var_degree.second);
        } else {
          tmp_variables.Insert(var_degree.first, var_degree.second);
        }
      }

      return CommutativeMonomial{std::move(tmp_variables)};
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

    bool operator<(const CommutativeMonomial &rhs) const {
      return variables_ < rhs.variables_;
    }

    bool operator==(const CommutativeMonomial &rhs) const {
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

    std::size_t Hash() const {
      std::size_t h = 0;
      for (auto &x : variables_) {
        HashCombine(h, x);
      }
      return h;
    }
};

namespace std {

template<>
struct hash<CommutativeMonomial> {
  inline std::size_t operator()(const CommutativeMonomial &m) const {
    return m.Hash();
  }
};

}  /* namespace std */


std::ostream& operator<<(std::ostream &out, const CommutativeMonomial &monomial);
