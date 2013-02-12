#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>

#include "matrix.h"
#include "semiring.h"
#include "var.h"

#ifndef OLD_FREESEMIRING
#include "free-semiring.h"
#else
#include "free-semiring-old.h"
#endif  /* OLD_FREESEMIRING */


/*
 * TODO:
 * - Write move constructors.
 * - Use move semantics to move all the temporary sets/maps, etc.
 */


/*
 * The std::map<VarPtr, Degree> should be a separate class that contains
 * the following functionality.
 */

typedef std::uint16_t Degree;

Degree GetDegreeOf(const std::map<VarPtr, Degree> &map,
    const VarPtr var);

void EraseAll(std::map<VarPtr, Degree> &map, const VarPtr var);

void Insert(std::map<VarPtr, Degree> &map, const VarPtr var, Degree deg = 1);

void Erase(std::map<VarPtr, Degree> &map, const VarPtr var, Degree deg = 1);

void Merge(std::map<VarPtr, Degree> &to_modify,
    const std::map<VarPtr, Degree> &to_merge);

//FIXME: Polynomials are no semiring in our definition (not starable)

template <typename SR>
class Polynomial;

template <typename SR>
class Monomial {
  private:
    SR coefficient_;
    /* Maps each variable to its degree. */
    std::map<VarPtr, Degree> variables_;

    friend class Polynomial<SR>;

    // private constructor to not leak the internal data structure
    Monomial(SR c, std::map<VarPtr, Degree> vs)
        : coefficient_(c), variables_(vs) {}

  public:
    /* Constant. */
    Monomial(SR c) : coefficient_(c) {}

    Monomial(SR c, std::initializer_list<VarPtr> vs)
        : coefficient_(c) {
      for (auto var : vs) {
        Insert(variables_, var);
      }
    }

    // Monomial(SR c, std::initializer_list< std::pair<VarPtr, Degree> > vs)
    //     : coefficient_(c) {
    //   for (auto var_degree : vs) {
    //     Insert(variables_, var_degree.first, var_degree.second);
    //   }
    // }

    /* std::vector seems to be a neutral data type and does not leak internal
     * data structure. */
    Monomial(SR c, std::vector< std::pair<VarPtr, Degree> > vs)
        : coefficient_(c) {
      for (auto var_degree : vs) {
        Insert(variables_, var_degree.first, var_degree.second);
      }
    }

    /* Add the coefficients of two monomials if their variables are equal. */
    Monomial operator+(const Monomial &monomial) const {
      assert(variables_ == monomial.variables_);
      return Monomial{coefficient_ + monomial.coefficient_, variables_};
    }

    // multiply two monomials
    Monomial operator*(const Monomial &monomial) const {
      if (is_null() || monomial.is_null()) {
        return Monomial{SR::null()};
      }
      auto tmp_variables = variables_;
      for (auto var_degree : monomial.variables_) {
        Insert(tmp_variables, var_degree.first, var_degree.second);
      }
      // FIXME: move instead of copying
      return Monomial(coefficient_ * monomial.coefficient_, tmp_variables);
    }

    // multiply a monomial with a variable
    Monomial operator*(const VarPtr &var) const {
      auto tmp_variables = variables_;
      // "add" the variables from one to the other monomial
      Insert(tmp_variables, var);
      // FIXME: move instead of copying
      return Monomial{coefficient_, tmp_variables};
    }

    // commutative version of derivative
    Monomial derivative(const VarPtr &var) const {

      auto var_iter = variables_.find(var);

      /* If the variable does not appear in the monomial, the derivative
       * must be 0. */
      if (var_iter == variables_.end()) {
        return Monomial(SR::null());
      }

      auto degree_before = var_iter->second;

      /* Remove one of these by removing the first of them and then "multiply"
       * the coefficient with degree_before. */
      auto tmp_variables = variables_;
      Erase(tmp_variables, var);

      SR tmp_coefficient = coefficient_;
      for(int i = 0; i < degree_before; ++i) {
        tmp_coefficient = tmp_coefficient + coefficient_;
      }
      // FIXME: move instead of copying
      return Monomial{tmp_coefficient, tmp_variables};
    }

    /* Evaluate the monomial given the map from variables to values. */
    SR eval(const std::map<VarPtr, SR> &values) const {
      SR result{coefficient_};

      for (auto var_degree : variables_) {
        auto value_iter = values.find(var_degree.first);
        /* All variables should be in the values map. */
        assert(value_iter != values.end());
        for (Degree i = 0; i < var_degree.second; ++i) {
          result = result * value_iter->second;
        }
      }

      return result;
    }

    /* Partially evaluate the monomial. */
    Monomial<SR> partial_eval(const std::map<VarPtr, SR> &values) const {
      Monomial<SR> result{coefficient_};

      for (auto var_degree : variables_) {
        auto e = values.find(var_degree.first);
        if (e == values.end()) {
          /* Variable not found in the mapping, so keep it. */
          Insert(result.variables_, var_degree.first, var_degree.second);
        } else {
          /* Variable found, use it for evaluation. */
          result.coefficient_ = result.coefficient_ * e->second;
        }
      }

      return result;
    }

    /* Variable substitution. */
    Monomial<SR> subst(const std::map<VarPtr, VarPtr> &mapping) const {
      auto tmp_variables = variables_;

      for(auto from_to : mapping) {
        auto degree = GetDegreeOf(tmp_variables, from_to.first);
        EraseAll(tmp_variables, from_to.first);
        Insert(tmp_variables, from_to.second, degree);
      }

      // FIXME: move instead of copying
      return Monomial{coefficient_, tmp_variables};
    }

    // convert this monomial to an element of the free semiring
    FreeSemiring make_free(std::unordered_map<SR, VarPtr, SR> *valuation) const {
      FreeSemiring result = FreeSemiring::one();
      for(auto var_degree : variables_) {
        FreeSemiring tmp{var_degree.first};
        for (Degree i = 0; i < var_degree.second; ++i) {
          result = result * tmp;
        }
      }

      // change the SR element to a constant in the free semiring
      auto elem = valuation->find(coefficient_);
      if(elem == valuation->end()) { // this is a new SR element
        // map 'zero' and 'one' element to respective free semiring element
        if (coefficient_ == SR::null()) {
          result = FreeSemiring::null() * result;
        } else if (coefficient_ == SR::one()) {
          result = FreeSemiring::one() * result;
        } else {
          // use a fresh constant - the constructor of Var::getVar() will do this
          VarPtr tmp_var = Var::getVar();
          FreeSemiring tmp{tmp_var};
          valuation->emplace(coefficient_, tmp_var);
          result = tmp * result;
        }
      } else {
        // this is an already known element
        result = (elem->second) * result;
      }
      return result;
    }

    /*
     * FIXME: both operator< and operator== are wrong!
     * Not taking into account the coefficient might lead to incorrect results.
     */

    // a monomial is smaller than another monomial if the variables are smaller
    bool operator<(const Monomial& monomial) const {
      return variables_ < monomial.variables_;
    }

    // a monomial is equal to another monomial if the variables are equal
    // warning: the coefficient will not be regarded --- FIXME:this is extremly dangerous (regarding set-implementation of polynomial)!
    bool operator==(const Monomial& monomial) const {
      return variables_ == monomial.variables_;
    }

    // both monomials are equal if they have the same variables and the same coefficient
    bool equal(const Monomial& monomial) const {
      return variables_ == monomial.variables_ &&
             coefficient_ == monomial.coefficient_;
    }

    Degree get_degree() const {
      Degree degree = 0;
      for (auto var_degree : variables_) {
        degree += var_degree.second;
      }
      return degree;
    }

    SR get_coeff() const {
      return coefficient_;
    }

    bool is_null() const {
      return coefficient_ == SR::null();
    }

    void set_coeff(SR coeff) {
      coefficient_ = coeff;
    }

    void add_to_coeff(const SR coeff) {
      coefficient_ = coefficient_ + coeff;
    }

    std::set<VarPtr> get_variables() const {
      std::set<VarPtr> set;
      for (auto var_degree : variables_) {
        set.insert(var_degree.first);
      }
      return set;
    }

    std::string string() const {
      std::stringstream ss;
      ss << coefficient_;
      if(!is_null()) {
        ss << "*" << variables_;
      }
      return ss.str();
    }
};

template <typename SR>
class Polynomial : public Semiring< Polynomial<SR> > {
  private:
    std::set<Monomial<SR> > monomials_;
    std::map<VarPtr, Degree> variables_;

    Polynomial(std::set< Monomial<SR> > &&ms, std::map<VarPtr, Degree> &&vs)
        : monomials_(std::move(ms)), variables_(std::move(vs)) {}

    // private constructor to hide the internal data structure
    Polynomial(const std::set< Monomial<SR> > &ms) : monomials_(ms) {
      for (auto &monomial : monomials_) {
        Merge(variables_, monomial.variables_);
      }

      // If there is a null-monomial, it should be at the front.
      // If the polynomial has more than one element, delete the null element
      if (monomials_.size() > 1 && monomials_.begin()->is_null()) {
        monomials_.erase(monomials_.begin());
      }
    }

    static std::shared_ptr<Polynomial<SR>> elem_null;
    static std::shared_ptr<Polynomial<SR>> elem_one;
  public:
    Polynomial() = default;

    Polynomial(std::initializer_list< Monomial<SR> > monomial_list) {
      /* FIXME: Is this check necessary? */
      if (monomial_list.size() == 0) {
        monomials_ = { Monomial<SR>{SR::null()} };
      } else {
        for(const auto &monomial : monomial_list) {
          auto iter = monomials_.find(monomial);
          if (iter == monomials_.end()) {
            monomials_.insert(monomial);
          } else {
          /* FIXME: can we use here std::move? */
            auto tmp_monomial = *iter;
            monomials_.erase(iter);
            monomials_.emplace(tmp_monomial + monomial);
          }
          Merge(variables_, monomial.variables_);
        }
      }
    }

    Polynomial(const Polynomial &p)
        : monomials_(p.monomials_), variables_(p.variables_) {}

    // create a 'constant' polynomial
    Polynomial(const SR &elem) {
      monomials_.emplace(Monomial<SR>{elem});
    }

    // create a polynomial which consists only of one variable
    Polynomial(VarPtr var) {
      monomials_.emplace(Monomial<SR>{SR::one(), {var}});
      Insert(variables_, var);
    }

    Polynomial<SR>& operator=(const Polynomial<SR> &p) {
      monomials_ = p.monomials_;
      variables_ = p.variables_;
      return *this;
    }

    Polynomial<SR>& operator+=(const Polynomial<SR> &polynomial) {
      for (const auto &monomial : polynomial.monomials_) {
        auto iter = monomials_.find(monomial);
        if (iter == monomials_.end()) {
          monomials_.insert(monomial);
        } else {
          /* FIXME: can we use here std::move? */
          Monomial<SR> tmp = *iter;
          monomials_.erase(iter);
          monomials_.emplace(tmp + monomial);
        }
      }

      return *this;
    }

    Polynomial<SR> operator*=(const Polynomial<SR> &rhs) {
      std::set< Monomial<SR> > tmp_monomials;
      std::map<VarPtr, Degree> tmp_variables;
      // iterate over both this and the poly polynomial
      for (const auto &lhs_monomial : monomials_) {
        for (const auto &rhs_monomial : rhs.monomials_) {
          Monomial<SR> tmp_monomial = lhs_monomial * rhs_monomial;
          auto iter = tmp_monomials.find(tmp_monomial);
          if (iter == tmp_monomials.end()) {
            tmp_monomials.emplace(std::move(tmp_monomial));
            /* Updating tmp_variables is only necessary if the monomial is not
             * already present in the tmp_monomials.  If it's already there all
             * the degrees should match. */
            Merge(tmp_variables, tmp_monomial.variables_);
          } else { // the monomial was already in the list. Add them.
            tmp_monomial = tmp_monomial + *iter;
            tmp_monomials.erase(iter);
            tmp_monomials.emplace(std::move(tmp_monomial));
          }
        }
      }

      monomials_ = std::move(tmp_monomials);
      variables_ = std::move(tmp_variables);
      return *this;
    }

    // multiplying a polynomial with a variable
    Polynomial<SR> operator*(const VarPtr &var) const {
      std::set< Monomial<SR> > tmp_monomials;
      for(const auto &monomial : monomials_) {
        tmp_monomials.emplace(monomial * var);
      }

      return Polynomial{tmp_monomials};
    }

    friend Polynomial<SR> operator*(const SR &elem,
        const Polynomial<SR> &polynomial) {

      std::set<Monomial<SR> > tmp_monomials;

      for (const auto &monomial : polynomial.monomials_) {
        tmp_monomials.emplace(elem * monomial);
      }

      return Polynomial{tmp_monomials};
    }

    bool operator==(const Polynomial<SR> &polynomial) const {

      return variables_ == polynomial.variables_ &&
             monomials_ == polynomial.monomials_;

      /*
      if (monomials_.size() != polynomial.monomials_.size()) {
        return false;
      }

      for (const auto &monomial : polynomial.monomials_) {
        auto iter = monomials_.find(monomial);
        if (iter == monomials_.end()) {
          return false; // could not find this monomial
        } else {
          if (!monomial.equal(*iter)) // check if the monomial has the same coefficient
            return false;
        }

      }
      return true;
      */
    }

    // convert the given matrix to a matrix containing polynomials
    // TODO: needed?
    /*
    static Matrix<Polynomial<SR> > convert(const Matrix<SR>& mat) {
      std::vector<Polynomial<SR> > ret;
      for(int i=0; i<mat.getColumns() * mat.getRows(); ++i) {
        // create constant polynomials
        ret.push_back(Polynomial(mat.getElements().at(i)));
      }

      return Matrix<Polynomial<SR> >(mat.getColumns(), mat.getRows(), ret);
    }
    */

    Polynomial<SR> derivative(const VarPtr& var) const {
      std::set< Monomial<SR> > tmp_monomials;
      std::map<VarPtr, Degree> tmp_variables;

      for (const auto &monomial : monomials_) {
        if (true) {
        // if (SR::is_commutative) {
          // take the derivative of m_it and add it to the result set
          Monomial<SR> derivative = monomial.derivative(var);
          auto iter = tmp_monomials.find(derivative);
          // TODO: think about this and remove if not needed
          if (iter == tmp_monomials.end()) {
            tmp_monomials.insert(derivative);
            Merge(tmp_variables, monomial.variables_);
          } else {
            Monomial<SR> tmp_monomial = *iter + derivative;
            tmp_monomials.erase(iter);
            tmp_monomials.emplace(std::move(tmp_monomial));
          }
        } else { // non-commutative case
          assert(false);  // TODO: not implemented yet
        }
      }

      return Polynomial{std::move(tmp_monomials), std::move(tmp_variables)};
      /*
      if (monomials.empty()) // TODO: save variables explicit in this class and check if var is in vars
        return Polynomial();
      else
        return Polynomial(monomials);
        */
    }

    Polynomial<SR> derivative(const std::vector<VarPtr> &vars) const {
      Polynomial<SR> tmp_polynomial = *this; // copy the polynomial
      for(auto var : vars) {
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
      return Matrix< Polynomial<SR> >{ variables.size(), polynomials.size(),
                                       std::move(result_vector) };
    };

    SR eval(const std::map<VarPtr, SR> &values) const {
      SR result = SR::null();
      for(const auto &monomial : monomials_) {
        result += monomial.eval(values);
      }
      return result;
    }

    // evaluate the polynomial with the given mapping and return the new polynomial
    Polynomial<SR> partial_eval(const std::map<VarPtr, SR> &values) const {
      Polynomial<SR> result = Polynomial<SR>::null();
      for(const auto &monomial : monomials_) {
        result += Polynomial({monomial.partial_eval(values)});
      }
      return result;
    }

    // substitute variables with other variables
    Polynomial<SR> subst(const std::map<VarPtr, VarPtr> &mapping) const {
      std::set< Monomial<SR> > tmp_monomials;

      for (const auto &monomial : monomials_) {
        tmp_monomials.emplace(monomial.subst(mapping));
      }

      return Polynomial<SR>{std::move(tmp_monomials)};
    }

    static Matrix<SR> eval(const Matrix< Polynomial<SR> > &poly_matrix,
        const std::map<VarPtr, SR> &values) {
      std::vector<Polynomial<SR> > tmp_polynomials = poly_matrix.getElements();
      std::vector<SR> result;
      for(const auto &polynomial : tmp_polynomials) {
        result.emplace_back(polynomial.eval(values));
      }
      return Matrix<SR>{poly_matrix.getColumns(), poly_matrix.getRows(), result};
    }

    // FIXME: Why values is passed by value???
    static Matrix<Polynomial<SR> > eval(Matrix<Polynomial<SR> > poly_matrix,
        std::map<VarPtr,Polynomial<SR> > values) {
      std::vector< Polynomial<SR> > tmp_polynomials = poly_matrix.getElements();
      std::vector< Polynomial<SR> > result;
      for (const auto &polynomial : tmp_polynomials) {
        result.emplace_back(polynomial.eval(values));
      }
      return Matrix< Polynomial<SR> >{poly_matrix.getColumns(),
                                      poly_matrix.getRows(), result};
    }

    // convert this polynomial to an element of the free semiring. regard the valuation map
    // which can already define a conversion from the SR element to a free SR constant
    // the valuation map is extended in this function
    FreeSemiring make_free(std::unordered_map<SR, VarPtr, SR> *valuation) const {
      assert(valuation);
      /*
      if (!valuation) {
        valuation = new std::unordered_map<SR, VarPtr, SR>();
      }
      */

      auto result = FreeSemiring::null();
      // convert this polynomial by adding all converted monomials
      for (const auto &monomial : monomials_) {
        result += monomial.make_free(valuation);
      }

      return result;
    }

    // convert this matrix of polynomials to a matrix with elements of the free semiring
    static Matrix<FreeSemiring> make_free(
        const Matrix< Polynomial<SR> > &poly_matrix,
        std::unordered_map<SR, VarPtr, SR> *valuation) {

      assert(valuation);

      std::vector< Polynomial<SR> > tmp_polynomials = poly_matrix.getElements();
      std::vector<FreeSemiring> result;

      // if(!valuation)
      //   valuation = new std::unordered_map<SR, VarPtr, SR>();

      for(const auto &polynomial : tmp_polynomials) {
        result.emplace_back(polynomial.make_free(valuation));
      }
      return Matrix<FreeSemiring>{poly_matrix.getColumns(),
                                  poly_matrix.getRows(), std::move(result)};
    }

    Degree get_degree() {
      Degree degree = 0;
      for (auto &monomial : monomials_) {
        degree = std::max(degree, monomial.get_degree());
      }
      return degree;
    }

    /* FIXME: Get rid of this. */
    std::set<VarPtr> get_variables() const {
      std::set<VarPtr> vars;
      for (auto var_degree : variables_) {
        vars.insert(var_degree.first);
      }
      return vars;
    }

    // some semiring functions
    Polynomial<SR> star() const {
      // TODO: we cannot star polynomials!
      assert(false);
      return (*this);
    }

    static Polynomial<SR> const null() {
      if (!Polynomial::elem_null)
        Polynomial::elem_null = std::make_shared< Polynomial<SR> >(SR::null());
      return *Polynomial::elem_null;
    }

    static Polynomial<SR> const one() {
      if (!Polynomial::elem_one)
        Polynomial::elem_one = std::make_shared< Polynomial<SR> >(SR::one());
      return *Polynomial::elem_one;
    }

    static bool is_idempotent;
    static bool is_commutative;

    std::string string() const {
      std::stringstream ss;
      for (auto iter = monomials_.begin(); iter != monomials_.end(); ++iter) {
        if(iter != monomials_.begin())
          ss << " + ";
        ss << *iter;
      }

      return ss.str();
    }
};

template <typename SR> bool Polynomial<SR>::is_commutative = false;
template <typename SR> bool Polynomial<SR>::is_idempotent = false;
// initialize pointers
template <typename SR> std::shared_ptr<Polynomial<SR>> Polynomial<SR>::elem_null;
template <typename SR> std::shared_ptr<Polynomial<SR>> Polynomial<SR>::elem_one;

  template <typename SR>
std::ostream& operator<<(std::ostream& os, const Monomial<SR>& monomial) {
  return  os << monomial.string();
}

  template <typename SR>
std::ostream& operator<<(std::ostream& os, const std::map<VarPtr, SR>& values) {
  for(auto value = values.begin(); value != values.end(); ++value) {
    os << value->first << "→" << value->second << ";";
  }
  return os;
}

template <typename SR>
std::ostream& operator<<(std::ostream& os, const std::map<VarPtr, Polynomial<SR> >& values) {
  for(auto value = values.begin(); value != values.end(); ++value) {
    os << value->first << "→" << value->second << ";";
  }
  return os;
}

#endif
