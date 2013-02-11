#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

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


//FIXME: Polynomials are no semiring in our definition (not starable)

template <typename SR>
class Monomial {
  private:
    SR coeff;
    std::multiset<VarPtr,VarPtrSort> variables;

    // private constructor to not leak the internal data structure
    Monomial(SR coeff, std::multiset<VarPtr,VarPtrSort> variables)
        : coeff(coeff), variables(variables) {}

  public:
    // constant monomial coeff
    Monomial(SR coeff) : coeff(coeff) {}

    Monomial(SR coeff, std::initializer_list<VarPtr> variables)
        : coeff(coeff), variables(variables) {}

    // std::vector seems to be a neutral data type and does not leak internal data structure
    Monomial(SR c, std::vector<VarPtr> vs) : coeff(c) {
      variables.insert(vs.begin(), vs.end());
    }

    // add the coefficients of two monomials if there variables are equal
    Monomial operator+(const Monomial& monomial) const {
      assert(this->variables == monomial.variables);
      return Monomial(this->coeff + monomial.coeff, this->variables);
    }

    // multiply two monomials
    Monomial operator*(const Monomial& monomial) const {
      if (is_null() || monomial.is_null()) {
        return Monomial{SR::null()};
      }
      std::multiset<VarPtr,VarPtrSort> variables = this->variables;
      // "add" the variables from one to the other monomial
      variables.insert(monomial.variables.begin(), monomial.variables.end());
      return Monomial(this->coeff * monomial.coeff, variables);
    }

    // multiply a monomial with a variable
    Monomial operator*(const VarPtr& var) const {
      std::multiset<VarPtr,VarPtrSort> variables = this->variables;
      // "add" the variables from one to the other monomial
      variables.insert(var);
      return Monomial(this->coeff, variables);
    }

    // commutative version of derivative
    Monomial derivative(const VarPtr& var) const {
      // count number of occurences of var in variables
      int count = this->variables.count(var);

      // variable is not in variables, derivative is null
      if(count == 0)
        return Monomial(SR::null());

      // remove one of these by removing the first of them and then "multiply"
      // the coefficient with count
      std::vector<VarPtr> variables(this->variables.begin(), this->variables.end());
      SR coeff = this->coeff;
      for(unsigned int i=0; i<variables.size(); ++i) {
        if(variables[i] == var) {
          variables.erase(variables.begin()+i);
          for(int j = 0; j<count-1; ++j)
            coeff = coeff + this->coeff;
          break;
        }
      }
      std::multiset<VarPtr,VarPtrSort> result(variables.begin(), variables.end());
      return Monomial(coeff, result);
    }

    // evaluate the monomial at the position values
    SR eval(const std::map<VarPtr, SR>& values) const {
      SR elem = this->coeff;

      for(auto v_it = this->variables.begin(); v_it != this->variables.end(); ++v_it) {
        auto e = values.find(*v_it);
        assert(e != values.end()); // all variables should be in values
        SR foo = e->second;
        elem = elem * foo;
      }

      return elem;
    }

    // partially evaluate the monomial at the position values
    Monomial<SR> partial_eval(const std::map<VarPtr, SR>& values) const {
      // keep the coefficient
      Monomial<SR> elem(this->coeff);

      // then loop over all variables and try to evaluate them
      for(auto v_it = this->variables.begin(); v_it != this->variables.end(); ++v_it) {
        auto e = values.find(*v_it);
        if(e == values.end()) { // variable not found in the mapping, so keep it
          elem.variables.insert(*v_it);
        }
        else { // variable was found, use it for evaluation
          elem.coeff = elem.coeff * (e->second);
        }
      }

      return elem;
    }

    // substitute variables with other variables
    Monomial<SR> subst(const std::map<VarPtr, VarPtr>& mapping) const {
      SR coeff = this->coeff;
      std::multiset<VarPtr,VarPtrSort> variables = this->variables;

      for(auto m_it = mapping.begin(); m_it != mapping.end(); ++m_it) {
        int count = variables.count(m_it->first);
        variables.erase( (*m_it).first ); // erase all occurences
        for(int i = 0; i < count; ++i)
          variables.insert( (*m_it).second ); // and "replace" them
      }

      return Monomial(coeff, variables);
    }

    // convert this monomial to an element of the free semiring
    FreeSemiring make_free(std::unordered_map<SR, VarPtr, SR>* valuation) const {
      FreeSemiring result = FreeSemiring::one();
      for(auto v_it = this->variables.begin(); v_it != this->variables.end(); ++v_it) {
        result = result * FreeSemiring(*v_it);
      }

      // change the SR element to a constant in the free semiring
      auto elem = valuation->find(this->coeff);
      if(elem == valuation->end()) { // this is a new SR element
        // map 'zero' and 'one' element to respective free semiring element
        if(this->coeff == SR::null()) {
          // valuation->insert(valuation->begin(), std::pair<SR,FreeSemiring>(this->coeff,FreeSemiring::null()));
          result = FreeSemiring::null() * result;
        }
        else if(this->coeff == SR::one()) {
          // valuation->insert(valuation->begin(), std::pair<SR,FreeSemiring>(this->coeff,FreeSemiring::one()));
          result = FreeSemiring::one() * result;
        } else {
          // use a fresh constant - the constructor of Var::getVar() will do this
          VarPtr tmp_var = Var::getVar();
          FreeSemiring tmp(tmp_var);
          valuation->insert(valuation->begin(), std::pair<SR, VarPtr>(this->coeff,tmp_var));
          result = tmp * result;
        }
      }
      else { // this is an already known element
        result = (elem->second) * result;
      }
      return result;
    }

    // a monomial is smaller than another monomial if the variables are smaller
    bool operator<(const Monomial& monomial) const {
      return this->variables < monomial.variables;
    }

    // a monomial is equal to another monomial if the variables are equal
    // warning: the coefficient will not be regarded --- FIXME:this is extremly dangerous (regarding set-implementation of polynomial)!
    bool operator==(const Monomial& monomial) const {
      return this->variables == monomial.variables;
    }

    // both monomials are equal if they have the same variables and the same coefficient
    bool equal(const Monomial& monomial) const {
      return this->variables == monomial.variables && this->coeff == monomial.coeff;
    }

    int get_degree() const {
      return this->variables.size();
    }

    SR get_coeff() const {
      return this->coeff;
    }

    bool is_null() const {
      return coeff == SR::null();
    }

    void set_coeff(SR coeff) {
      this->coeff = coeff;
    }

    void add_to_coeff(const SR coeff) {
      this->coeff = this->coeff + coeff;
    }

    std::set<VarPtr> get_variables() const {
      return std::set<VarPtr>(this->variables.begin(), this->variables.end());
    }

    std::string string() const {
      std::stringstream ss;
      ss << this->coeff;
      if(!(this->coeff == SR::null() || this->variables.empty())) {
        ss << "*" << this->variables;
      }
      return ss.str();
    }
};

template <typename SR>
class Polynomial : public Semiring<Polynomial<SR> > {
  private:
    int degree;
    std::set<Monomial<SR> > monomials;
    std::set<VarPtr> variables;

    // private constructor to hide the internal data structure
    Polynomial(const std::set<Monomial<SR> >& monomials) {
      this->monomials = monomials;
      this->degree = 0;
      for(auto m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it) {
        this->degree = (*m_it).get_degree() > this->degree ? (*m_it).get_degree() : this->degree;
        auto vars = m_it->get_variables();
        this->variables.insert(vars.begin(), vars.end()); // collect all used variables
      }

      // If there is a null-monomial, it should be at the front.
      // If the polynomial has more than one element, delete the null element
      if(this->monomials.size() > 1 && this->monomials.begin()->is_null()) {
        this->monomials.erase(this->monomials.begin());
      }
    }

    static std::shared_ptr<Polynomial<SR>> elem_null;
    static std::shared_ptr<Polynomial<SR>> elem_one;
  public:
    // empty polynomial
    Polynomial() {
      this->degree = 0;
    };

    Polynomial(std::initializer_list<Monomial<SR> > monomials) {
      if(monomials.size() == 0) {
        this->monomials = {Monomial<SR>(SR::null(),{})};
        this->degree = 0;
      }
      else {
        this->monomials = std::set<Monomial<SR> >();
        this->degree = 0;

        for(auto m_it = monomials.begin(); m_it != monomials.end(); ++m_it) {
          auto mon = this->monomials.find(*m_it);

          if(mon == this->monomials.end()) // not yet present as a monomial
            this->monomials.insert(*m_it); // just insert it
          else {
            //monomial already present in this polynomial---add the coefficients (note that c++ containers cannot be modified in-place!)
            Monomial<SR> tmp = *mon;
            this->monomials.erase(mon);
            this->monomials.insert( tmp + (*m_it) );
          }
          this->degree = m_it->get_degree() > this->degree ? m_it->get_degree() : this->degree;
          auto vars = m_it->get_variables();
          this->variables.insert(vars.begin(), vars.end()); // collect all used variables
        }
      }
    }

    Polynomial(const Polynomial& polynomial) {
      this->monomials = polynomial.monomials;
      this->degree = polynomial.degree;
      this->variables = polynomial.variables;
    }

    // create a 'constant' polynomial
    Polynomial(const SR& elem) {
      this->monomials = {Monomial<SR>(elem,{})};
      this->degree = 0;
    }

    // create a polynomial which consists only of one variable
    Polynomial(VarPtr var) {
      this->monomials = {Monomial<SR>(SR::one(),{var})};
      this->degree = 1;
    }

    Polynomial& operator=(const Polynomial& polynomial) {
      this->monomials = polynomial.monomials;
      this->degree = polynomial.degree;
      this->variables = polynomial.variables;
      return (*this);
    }

    Polynomial<SR> operator+=(const Polynomial<SR>& poly) {
      std::set<Monomial<SR> > monomials = this->monomials;
      for(auto m_it = poly.monomials.begin(); m_it != poly.monomials.end(); ++m_it) {
        // check if "same" monomial is already in the set
        auto mon = monomials.find(*m_it);
        if(mon == monomials.end()) // this is not the case
          monomials.insert(*m_it); // just insert it
        else { // monomial with the same variables found
          Monomial<SR> tmp = *mon;
          monomials.erase(mon);
          monomials.insert( tmp + (*m_it) ); // then add both of them and overwrite the old one
        }
      }

      *this = Polynomial(monomials);
      return *this;
    }

    Polynomial<SR> operator*=(const Polynomial<SR>& poly) {
      std::set<Monomial<SR> > monomials;
      // iterate over both this and the poly polynomial
      for(auto m_it1 = this->monomials.begin(); m_it1 != this->monomials.end(); ++m_it1) {
        for(auto m_it2 = poly.monomials.begin(); m_it2 != poly.monomials.end(); ++m_it2) {
          Monomial<SR> tmp = (*m_it1) * (*m_it2);
          auto tmp2 = monomials.find(tmp);
          if(tmp2 == monomials.end())
            monomials.insert(tmp); // multiply them and insert them to the result set
          else { // the monomial was already in the list. Add them.
            tmp = tmp + *tmp2;
            monomials.erase(*tmp2);
            monomials.insert(tmp);
          }
        }
      }

      *this = Polynomial(monomials);
      return *this;
    }

    // multiplying a polynomial with a variable
    Polynomial<SR> operator*(const VarPtr& var) const {
      std::set<Monomial<SR> > monomials;
      for(auto m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it) {
        monomials.insert( (*m_it) * var );
      }

      return Polynomial(monomials);
    }

    friend Polynomial<SR> operator*(const SR& elem, const Polynomial<SR>& polynomial) {
      std::set<Monomial<SR> > monomials;

      for(auto m_it = polynomial.monomials.begin(); m_it != polynomial.monomials.end(); ++m_it) {
        monomials.insert(elem * (*m_it));
      }

      return Polynomial(monomials);
    }

    bool operator==(const Polynomial<SR>& polynomial) const {
      if(this->monomials.size() != polynomial.monomials.size())
        return false;

      for(auto m_it = polynomial.monomials.begin(); m_it != polynomial.monomials.end(); ++m_it) {
        auto monomial = this->monomials.find(*m_it); // search with variables
        if(monomial == this->monomials.end())
          return false; // could not find this monomial
        else {
          if(!monomial->equal(*m_it)) // check if the monomial has the same coefficient
            return false;
        }

      }
      return true;
    }

    // convert the given matrix to a matrix containing polynomials
    // TODO: needed?
    static Matrix<Polynomial<SR> > convert(const Matrix<SR>& mat) {
      std::vector<Polynomial<SR> > ret;
      for(int i=0; i<mat.getColumns() * mat.getRows(); ++i) {
        // create constant polynomials
        ret.push_back(Polynomial(mat.getElements().at(i)));
      }

      return Matrix<Polynomial<SR> >(mat.getColumns(), mat.getRows(), ret);
    }

    Polynomial<SR> derivative(const VarPtr& var) const {
      std::set<Monomial<SR> > monomials;

      for(auto m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it) {
        //if(SR::is_commutative) // TODO: check if compiler is optimizing this out
        if(true) {
          // take the derivative of m_it and add it to the result set
          Monomial<SR> derivative = (*m_it).derivative(var);
          auto monomial = monomials.find(derivative);
          if(monomial == monomials.end()) { // TODO: think about this and remove if not needed
            monomials.insert(derivative);
          } else {
            Monomial<SR> tmp = (*monomial) + derivative;
            monomials.erase(monomial); // remove
            monomials.insert(tmp); // and insert the updated version
          }
        }
        else { // non-commutative case
          assert(false);  // TODO: not implemented yet
        }
      }

      if(monomials.empty()) // TODO: save variables explicit in this class and check if var is in vars
        return Polynomial();
      else
        return Polynomial(monomials);
    }

    Polynomial<SR> derivative(const std::vector<VarPtr>& vars) const {
      Polynomial<SR> polynomial = *this; // copy the polynomial
      for(auto var = vars.begin(); var != vars.end(); ++var) {
        polynomial = polynomial.derivative(*var);
      }
      return polynomial;
    }

    static Matrix<Polynomial<SR> > jacobian(const std::vector<Polynomial<SR> >& polynomials, const std::vector<VarPtr>& variables) {
      std::vector<Polynomial<SR> > ret;
      for(auto poly = polynomials.begin(); poly != polynomials.end(); ++poly) {
        for(auto var = variables.begin(); var != variables.end(); ++var) {
          ret.push_back(poly->derivative(*var));
        }
      }
      return Matrix<Polynomial<SR> >(variables.size(), polynomials.size(), ret);
    };

    // TODO: i need the variables in this function!
    Matrix<Polynomial<SR> > hessian() const {
      std::vector<Polynomial<SR> > ret;
      for(auto var2 = this->variables.begin(); var2 != this->variables.end(); ++var2) {
        Polynomial<SR> tmp = this->derivative(*var2);
        for(auto var1 = this->variables.begin(); var1 != this->variables.end(); ++var1) {
          ret.push_back(tmp.derivative(*var1));
        }
      }
      return Matrix<Polynomial<SR> >(this->variables.size(), this->variables.size(), ret);
    }

    SR eval(const std::map<VarPtr,SR>& values) const {
      SR result = SR::null();
      for(auto m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it) {
        SR elem = m_it->eval(values);
        result = result + elem;
      }
      return result;
    }

    // evaluate the polynomial with the given mapping and return the new polynomial
    Polynomial<SR> partial_eval(const std::map<VarPtr,SR>& values) const {
      Polynomial<SR> result = Polynomial<SR>::null();
      for(auto m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it) {
        Monomial<SR> elem = m_it->partial_eval(values);
        result = result + Polynomial({elem});
      }
      return result;
    }

    // substitute variables with other variables
    Polynomial<SR> subst(const std::map<VarPtr, VarPtr>& mapping) const {
      std::set<Monomial<SR> > monomials;

      for(auto m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it)
        monomials.insert((*m_it).subst(mapping));

      return Polynomial<SR>(monomials);
    }

    static Matrix<SR> eval(const Matrix<Polynomial<SR> >& polys, const std::map<VarPtr,SR>& values) {
      std::vector<Polynomial<SR> > polynomials = polys.getElements();
      std::vector<SR> ret;
      for(int i = 0; i < polys.getRows()*polys.getColumns(); i++) {
        ret.push_back(polynomials[i].eval(values));
      }
      return Matrix<SR>(polys.getColumns(), polys.getRows(), ret);
    }

    static Matrix<Polynomial<SR> > eval(Matrix<Polynomial<SR> > polys, std::map<VarPtr,Polynomial<SR> > values) {
      std::vector<Polynomial<SR> > polynomials = polys.getElements();
      std::vector<Polynomial<SR> > ret;
      for(int i = 0; i < polys.getRows()*polys.getColumns(); i++) {
        ret.push_back(polynomials[i].eval(values));
      }
      return Matrix<Polynomial<SR> >(polys.getColumns(), polys.getRows(), ret);
    }

    // convert this polynomial to an element of the free semiring. regard the valuation map
    // which can already define a conversion from the SR element to a free SR constant
    // the valuation map is extended in this function
    FreeSemiring make_free(std::unordered_map<SR,VarPtr, SR>* valuation) {
      if(!valuation)
        valuation = new std::unordered_map<SR, VarPtr, SR>();

      FreeSemiring result = FreeSemiring::null();
      // convert this polynomial by adding all converted monomials
      for(auto m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it) {
        result = result + m_it->make_free(valuation);
      }

      return result;
    }

    // convert this matrix of polynomials to a matrix with elements of the free semiring
    static Matrix<FreeSemiring> make_free(const Matrix<Polynomial<SR> >& polys, std::unordered_map<SR, VarPtr, SR>* valuation)
    {
      std::vector<Polynomial<SR> > polynomials = polys.getElements();
      std::vector<FreeSemiring> ret;
      if(!valuation)
        valuation = new std::unordered_map<SR, VarPtr, SR>();

      for(int i = 0; i < polys.getRows()*polys.getColumns(); i++) {
        ret.push_back(polynomials[i].make_free(valuation));
      }
      return Matrix<FreeSemiring>(polys.getColumns(), polys.getRows(), ret);
    }

    int get_degree() {
      return this->degree;
    }

    std::set<VarPtr> get_variables() const {
      return this->variables;
    }

    // some semiring functions
    Polynomial<SR> star() const {
      // TODO: we cannot star polynomials!
      assert(false);
      return (*this);
    }

    static Polynomial<SR> const null() {
      if(!Polynomial::elem_null)
        Polynomial::elem_null = std::shared_ptr<Polynomial<SR>>(new Polynomial(SR::null()));
      return *Polynomial::elem_null;
    }

    static Polynomial<SR> const one() {
      if(!Polynomial::elem_one)
        Polynomial::elem_one = std::shared_ptr<Polynomial<SR>>(new Polynomial(SR::one()));
      return *Polynomial::elem_one;
    }

    static bool is_idempotent;
    static bool is_commutative;

    std::string string() const {
      std::stringstream ss;
      for (auto m_it = this->monomials.begin(); m_it != this->monomials.end(); ++m_it) {
        if(m_it != this->monomials.begin())
          ss << " + ";
        ss << (*m_it);
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
