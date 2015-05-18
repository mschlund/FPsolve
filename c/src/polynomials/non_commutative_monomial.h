#pragma once

#include <iosfwd>
#include <map>

#include "../datastructs/var_degree_map.h"
#include "../semirings/free-semiring.h"
#include "commutative_polynomial.h"


template <typename SR> class NonCommutativePolynomialBase;
template <typename SR> class NonCommutativeMonomial;

enum elemType {Variable, SemiringType};


template <typename SR>
class NonCommutativeMonomialBase {

protected:
  friend class NonCommutativeMonomial<SR>;
  friend class NonCommutativePolynomialBase<SR>;

  /* a monomial is represented by three vectors
   * idx_ holds pairs of (type, index), type is either variable or sr
   *   the index gives absolute position inside the chosen type vector
   * variables_ holds all variables
   * srs_ holds all semiring elements */
  std::vector<std::pair<elemType,int>> idx_;
  std::vector<VarId> variables_;
  std::vector<SR> srs_;


  /* Private constructor to not leak the internal data structure. */
  //NonCommutativeMonomialBase(VarDegreeMap &&vs) : variables_(std::move(vs)) {}
  NonCommutativeMonomialBase(std::vector<std::pair<elemType,int>> &idx,
      std::vector<VarId> &variables,
      std::vector<SR> &srs) :
        idx_(idx), variables_(variables), srs_(srs) {}


  /*
   * evaluate all variables except the one at the exceptional position
   * returns a monomial of the form aXb where a,b are in SR.
   * TODO: check, that all variables are present in the valuation!
   */
  NonCommutativeMonomialBase<SR> eval_all_except(const int except_pos, const ValuationMap<SR>& valuation) const {
    SR a = SR::one();
    SR b = SR::one();

    assert((unsigned int)except_pos < this->variables_.size());

    /*
     * go through the idx_-vector,
     * - evaluate all variables encountered as long as its position is less than except_pos
     * - on the way: multiply everything (evaluated vars and sr-elems) onto a
     * - after except_pos has been found: multiply everything encountered onto b
     */

    SR* prod = &a;
    for(auto p : this->idx_) {
      if(p.first == elemType::Variable && except_pos == p.second) {
        prod = &b;
        continue;
      }
      if(p.first == elemType::Variable)
        (*prod) *= valuation.at(this->variables_.at(p.second));
      else
        (*prod) *= this->srs_.at(p.second);
    }

    std::vector<std::pair<elemType,int> > info{std::make_pair(elemType::SemiringType,0),
      std::make_pair(elemType::Variable,0),
      std::make_pair(elemType::SemiringType,1),};
    std::vector<VarId> vars{this->variables_.at(except_pos)};
    std::vector<SR> sr{a,b};

    return NonCommutativeMonomialBase(info, vars, sr);
  }

public:

  /* Since we don't manage any resources, we can simply use the default
   * constructors and assignment operators. */
  NonCommutativeMonomialBase() = default;
  NonCommutativeMonomialBase(const NonCommutativeMonomialBase &m) = default;
  NonCommutativeMonomialBase(NonCommutativeMonomialBase &&m) = default;

  NonCommutativeMonomialBase(std::initializer_list<std::pair<elemType,int>> idx,
      std::initializer_list<VarId> variables,
      std::initializer_list<SR> srs) :
        idx_(idx), variables_(variables), srs_(srs_) {}


  NonCommutativeMonomialBase& operator=(const NonCommutativeMonomialBase &m) = default;
  NonCommutativeMonomialBase& operator=(NonCommutativeMonomialBase &&m) = default;


  /* Multiply two monomials and return a result in normalform if both monomials
   * already are in normal form */
  NonCommutativeMonomialBase operator*(const NonCommutativeMonomialBase &monomial) const {
    auto tmp_idx = this->idx_;
    auto tmp_variables = this->variables_;
    auto tmp_srs = this->srs_;
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
        tmp_idx.push_back({Variable, p.second + offset_variables});
      else if (p.first == SemiringType)
        tmp_idx.push_back({SemiringType, p.second + offset_srs});
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

    //return NonCommutativeMonomialBase(std::move(tmp_idx), std::move(tmp_variables), std::move(tmp_srs));
    return NonCommutativeMonomialBase(tmp_idx, tmp_variables, tmp_srs);
  }


  /* Multiply a monomial with a variable. */
  NonCommutativeMonomialBase operator*(const VarId &var) const {
    auto tmp_idx = this->idx_;
    auto tmp_variables = this->variables_;
    auto tmp_srs = this->srs_;
    tmp_idx.push_back(std::pair<elemType, int>(Variable, tmp_variables.size()));
    tmp_variables.push_back(var);

    return NonCommutativeMonomialBase(tmp_idx, tmp_variables, tmp_srs);
  }

  /* Multiply a monomial with a semiring element and return it in normal form. */
  NonCommutativeMonomialBase operator*(const SR &sr) const {
    auto tmp_idx = this->idx_;
    auto tmp_variables = this->variables_;
    auto tmp_srs = this->srs_;

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

    return NonCommutativeMonomialBase(tmp_idx, tmp_variables, tmp_srs);
  }

  /* convert this monomial into an commutative one
   * return a Polynomial, because the semiring element is not saved
   * in the commutative version of the monomial but in the polynomial */
  CommutativePolynomial<SR> make_commutative() const {
    CommutativePolynomial<SR> result_polynomial = CommutativePolynomial<SR>::one();

    for(auto const &sr : this->srs_) {
      result_polynomial *= sr;
    }

    for(auto const &var : this->variables_) {
      result_polynomial *= var;
    }

    return result_polynomial;
  }

  /* derivative function which is used in the polynomial derivative function.
   * for the variables for the d-1-th iteration we use the given map 'substitution' */
  NonCommutativePolynomialBase<SR> derivative(const SubstitutionMap &substitution) const {
    NonCommutativePolynomialBase<SR> result; // empty polynomial
    auto subst_monomial = this->subst(substitution); // substitute all variables
    for(unsigned int position = 0; position < this->variables_.size(); position++) {
      // the variable at the position is the variable which should not be touched
      auto tmp = subst_monomial;
      tmp.variables_.at(position) = this->variables_.at(position); // therefore restore this one...
      result += tmp;
    }
    return result;
  }

  /*
   * compute the linearization of the polynomial at the given point "valuation"
   * to this end, we linearize every monomial and sum them up
   */
  NonCommutativePolynomialBase<SR> differential_at(const ValuationMap<SR> &valuation) const {
    NonCommutativePolynomialBase<SR> result; // empty polynomial = 0 (constant monomials have differential f(x)=0)

    for(unsigned int position = 0; position < this->variables_.size(); position++) {
      // the variable at the position is the variable which should not be touched, all others will be evaluated
      result += this->eval_all_except(position, valuation);
    }
    return result;
  }

  SR calculate_delta_helper(
      const std::vector<bool> &permutation,
      const ValuationMap<SR> &de2, // [d-2], true
      const ValuationMap<SR> &dl1  // (d-1), false
  ) const {
    SR tmp = SR::one();
    for(auto const &p : this->idx_) {
      if(p.first == Variable) {
        if(permutation.at(p.second) == true) { // use [d-2]
          auto value_iter = de2.find(this->variables_.at(p.second));
          assert(value_iter != de2.end());
          tmp *= value_iter->second;
        } else { // if (permutation.at(p.second) == false) // use (d-1)
          auto value_iter = dl1.find(this->variables_.at(p.second));
          assert(value_iter != dl1.end());
          tmp *= value_iter->second;
        }
      } else if (p.first == SemiringType)
        tmp *= this->srs_.at(p.second);
    }
    return tmp;
  }

  SR calculate_delta(
      const ValuationMap<SR> &de2, // [d-2], true
      const ValuationMap<SR> &dl1  // (d-1), false
  ) const {
    SR result = SR::null();

    /* outer loop handles the different cases (trees with exactly n-times dim == d-1 )
     * start with n = 2, which means, exactly 2 children have dimensions exactly d-1 */
    for(unsigned int n = 2; n <= this->variables_.size(); n++)
    {
      /* order of a vector of bools is [false, true] < [true, false] */
      std::vector<bool> permutation(n, false); // these are the (d-1) elements
      std::vector<bool> permutation2(this->variables_.size()-n, true); // these are the [d-2] elements
      permutation.insert(permutation.end(), permutation2.begin(), permutation2.end());
      do {
        result += calculate_delta_helper(permutation, de2, dl1);
      } while(std::next_permutation(permutation.begin(), permutation.end()));
    }

    return result;
  }

  /* Evaluate the monomial given the map from variables to values. */
  SR eval(const ValuationMap<SR> &values) const {
    auto result = SR::one();

    for (const auto &p : this->idx_)
    {
      if (p.first == Variable)
      {
        auto value_iter = values.find(this->variables_.at(p.second));
        /* All variables should be in the values map. */
        assert(value_iter != values.end());
        result *= value_iter->second;
      }
      else if (p.first == SemiringType)
        result *= this->srs_.at(p.second);
    }

    return result;
  }

  /* Partially evaluate the monomial. */
  NonCommutativeMonomialBase partial_eval(
      const ValuationMap<SR> &values) const {

    NonCommutativeMonomialBase result_monomial;

    for(auto p : this->idx_)
    {
      if (p.first == Variable)
      {
        auto value_iter = values.find(this->variables_.at(p.second));
        if (value_iter == values.end())
        { /* Variable not found in the mapping, so keep it. */
          result_monomial.idx_.push_back(std::pair<elemType, int>(Variable, result_monomial.variables_.size()));
          result_monomial.variables_.push_back(this->variables_.at(p.second));
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
          result_monomial.srs_.push_back(elem * this->srs_.at(p.second));

        } else {
          result_monomial.idx_.push_back(std::pair<elemType, int>(SemiringType, result_monomial.srs_.size()));
          result_monomial.srs_.push_back(this->srs_.at(p.second));
        }
      }
    }

    return result_monomial;
  }

  /* Variable substitution. */
  NonCommutativeMonomialBase subst(const SubstitutionMap &mapping) const {
    VarDegreeMap tmp_variables;

    auto result_monomial = *this; // copy it to work with it

    for(auto p : result_monomial.idx_)
    {
      if(p.first == Variable)
      {
        auto old_new_iter = mapping.find(this->variables_.at(p.second));
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

    for (auto p : this->idx_) {
      if (p.first == Variable) {
        result *= FreeSemiring(this->variables_.at(p.second));
      } else if (p.first == SemiringType) {
        {
          auto tmp_sr = this->srs_.at(p.second);
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
              FreeSemiring tmp_var_free{tmp_var};
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

  bool operator<(const NonCommutativeMonomialBase &rhs) const {
    // lexicographic ordering
    if(this->idx_ != rhs.idx_) return this->idx_ < rhs.idx_;
    if(this->variables_ != rhs.variables_) return this->variables_ < rhs.variables_;
    if(this->srs_ != rhs.srs_) return this->srs_ < rhs.srs_;

    // they are equal
    return false;
  }

  bool operator==(const NonCommutativeMonomialBase &rhs) const {
    return
        this->idx_ == rhs.idx_ &&
        this->variables_ == rhs.variables_ &&
        this->srs_ == rhs.srs_;
  }

  Degree get_degree() const {
    return this->variables_.size();
  }

  // FIXME: modify or remove
  std::set<VarId> get_variables() const {
    std::set<VarId> set;
    for (auto var : this->variables_) {
      set.insert(var);
    }
    return set;
  }

  /*
   * If the monomial has form xYz where x is an element of the semiring,
   * Y is a monomial over the semiring, and z is an element of the semiring,
   * this method gives back x wrapped in a monomial.
   */
  SR getLeadingSR() const {
    SR leadingSR = SR::one();

    // give me a copy of the semiring factor on the leading side of the monomial
    for(unsigned int i = 0; i < this->idx_.size(); i++) {

      // multiply as long as there was no variable encountered
      if(this->idx_.at(i).first == SemiringType) {
        leadingSR = leadingSR  * this->srs_.at(this->idx_.at(i).second);
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
    SR trailingSR = SR::one();

    // give me a copy of the semiring factor on the trailing side of the monomial
    for(int i = this->idx_.size() - 1; i >= 0; i--) {

      // multiply as long as there was no variable encountered
      if(this->idx_.at(i).first == SemiringType) {
        trailingSR = this->srs_.at(this->idx_.at(i).second) * trailingSR;
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

  std::string string() const {
    std::stringstream ss;
    for(auto p = this->idx_.begin(); p != this->idx_.end(); p++) {
      if(p != this->idx_.begin())
        ss << " * ";

      if(p->first == Variable)
        ss << this->variables_.at(p->second);
      else
        ss << this->srs_.at(p->second);
    }
    return std::move(ss.str());
  }
};

template <typename SR>
std::ostream& operator<<(std::ostream &out, const NonCommutativeMonomialBase<SR> &monomial) {
  return out << monomial.string();
}


template <typename SR>
class NonCommutativeMonomial : public NonCommutativeMonomialBase<SR> {
public:
  using NonCommutativeMonomialBase<SR>::NonCommutativeMonomialBase;

  // constructor for implicit conversions
  NonCommutativeMonomial(const NonCommutativeMonomialBase<SR>& p) {
    this->idx_ = p.idx_;
    this->srs_ = p.srs_;
    this->variables_ = p.variables_;
  }
};





















