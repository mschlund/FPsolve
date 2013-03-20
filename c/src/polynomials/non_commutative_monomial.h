#pragma once

#include <iosfwd>

#include "../semirings/free-semiring.h"

template <typename SR>
class NonCommutativePolynomial;
enum elemType {Variable, SemiringType};

template <typename SR>
class NonCommutativeMonomial {
  private:
    friend class NonCommutativePolynomial<SR>;
    
    /* a monomial is represented by three vectors
     * idx_ holds pairs of (type, index), type is either variable or sr
     *   the index gives absolute position inside the chosen type vector
     * variables_ holds all variables
     * srs_ holds all semiring elements */
    std::vector<std::pair<elemType,int>> idx_;
    std::vector<VarId> variables_;
    std::vector<SR> srs_;


    /* Private constructor to not leak the internal data structure. */
    //NonCommutativeMonomial(VarDegreeMap &&vs) : variables_(std::move(vs)) {}
    NonCommutativeMonomial(std::vector<std::pair<elemType,int>> &&idx,
        std::vector<VarId> &&variables,
        std::vector<SR> &&srs) :
      idx_(std::move(idx)), variables_(std::move(variables)), srs_(std::move(srs)) {}

  public:

    /* Since we don't manage any resources, we can simply use the default
     * constructors and assignment operators. */
    NonCommutativeMonomial() = default;
    NonCommutativeMonomial(const NonCommutativeMonomial &m) = default;
    NonCommutativeMonomial(NonCommutativeMonomial &&m) = default;
    NonCommutativeMonomial& operator=(const NonCommutativeMonomial &m) = default;
    NonCommutativeMonomial& operator=(NonCommutativeMonomial &&m) = default;

    NonCommutativeMonomial(std::initializer_list<std::pair<elemType,int>> idx,
        std::initializer_list<VarId> variables,
        std::initializer_list<SR> srs) :
      idx_(idx), variables_(variables), srs_(srs_) {}

    /* std::vector seems to be a neutral data type and does not leak internal
     * data structure. */
/*    NonCommutativeMonomial(std::vector< std::pair<elemType, int> > vs) {
      for (auto var_degree : vs) {
        variables_.Insert(var_degree.first, var_degree.second);
      }
    }
*/

    /* Multiply two monomials. */
    NonCommutativeMonomial operator*(const NonCommutativeMonomial &monomial) const {
      auto tmp_idx = idx_;
      auto tmp_variables = variables_;
      auto tmp_srs = srs_;
      unsigned int offset_variables = tmp_variables.size();
      unsigned int offset_srs = tmp_srs.size();

      for (auto &p : monomial.idx_) {
        if(p.first == Variable)
          p.second += offset_variables;
        else if (p.first == SemiringType)
          p.second += offset_srs;
        tmp_idx.push_back(p);
      }

      for (auto v : monomial.variables_)
        tmp_variables.push_back(v);

      for (auto s : monomial.srs_)
        tmp_srs.push_back(s);

      return NonCommutativeMonomial(std::move(tmp_idx), std::move(tmp_variables), std::move(tmp_srs));
    }


    /* Multiply a monomial with a variable. */
    NonCommutativeMonomial operator*(const VarId &var) const {
      auto tmp_idx = idx_;
      auto tmp_variables = variables_;
      auto tmp_srs = srs_;
      tmp_variables.push_back(var);
      tmp_idx.push_back(std::pair<elemType, int>(Variable, tmp_variables.size()));

      return NonCommutativeMonomial(std::move(tmp_idx), std::move(tmp_variables), std::move(tmp_srs));
    }

    /* Multiply a monomial with a semiring element. */
    // TODO: should we optimize two consecutive semiring elements?
    NonCommutativeMonomial operator*(const SR &sr) const {
      auto tmp_idx = idx_;
      auto tmp_variables = variables_;
      auto tmp_srs = srs_;
      tmp_srs.push_back(sr);
      tmp_idx.push_back(std::pair<elemType, int>(SemiringType, tmp_srs.size()));

      return NonCommutativeMonomial(std::move(tmp_idx), std::move(tmp_variables), std::move(tmp_srs));
    }

    /* Commutative version of derivative. */
/*    std::pair<Degree, NonCommutativeMonomial> derivative(const VarId &var) const {

      auto var_degree_iter = variables_.find(var);
*/
      /* If the variable does not appear in the monomial, the derivative
       * must be 0. */
/*      if (var_degree_iter == variables_.end()) {
        return { 0, NonCommutativeMonomial{} };
      }
*/
      /* Remove one of these by removing the first of them and then "multiply"
       * the coefficient with degree_before. */
/*      auto tmp_variables = variables_;
      tmp_variables.Erase(var);

      return { var_degree_iter->second, NonCommutativeMonomial{std::move(tmp_variables)} };
    }
*/

    /* Evaluate the monomial given the map from variables to values. */
    SR eval(const std::map<VarId, SR> &values) const {
      auto result = SR::one();

      for (const auto &p : idx_)
      {
        if (p.first == Variable)
        {
          auto value_iter = values.find(variables_.at(p.second));
          /* All variables should be in the values map. */
          assert(value_iter != values.end());
          result *= value_iter->second;
        }
        else if (p.first == SemiringType)
          result *= srs_.at(p.second);

      }

      return result;
    }

    /* Partially evaluate the monomial. */
    NonCommutativeMonomial partial_eval(
        const std::map<VarId, SR> &values) const {

      NonCommutativeMonomial result_monomial;

      for(auto p : idx_)
      {
        if (p.first == Variable)
        {
          auto value_iter = values.find(variables_.at(p.second));
          if (value_iter == values.end())
          { /* Variable not found in the mapping, so keep it. */
            result_monomial.idx_.push_back(std::pair<elemType, int>(Variable, result_monomial.variables_.size()));
            result_monomial.variables_.push_back(variables_.at(p.second));
          } else
          { /* Variable found, use it for evaluation */
            // FIXME: multiply SR elements immediatly if possible
            result_monomial.idx_.push_back(std::pair<elemType, int>(SemiringType, result_monomial.srs_.size()));
            result_monomial.srs_.push_back(value_iter->second);
          }
          
        } else if (p.first == SemiringType)
        {
          // FIXME: probably possible to multiply with new semiring element in front of this element
          result_monomial.idx_.push_back(std::pair<elemType, int>(SemiringType, result_monomial.srs_.size()));
          result_monomial.srs_.push_back(srs_.at(p.second));
        }
      }

      return result_monomial;
    }

    /* Variable substitution. */
    NonCommutativeMonomial subst(const std::map<VarId, VarId> &mapping) const {
      VarDegreeMap tmp_variables;

      auto result_monomial = *this; // copy it to work with it

      for(auto p : result_monomial.idx_)
      {
        if(p.first == Variable)
        {
          auto old_new_iter = mapping.find(variables_.at(p.second));
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

      for (auto p : idx_) {
        if (p.first == Variable) {
          result *= FreeSemiring(variables_.at(p.second));
        } else if (p.first == SemiringType) {
          {
            auto tmp_sr = srs_.at(p.second);
            if(tmp_sr == SR::null()) {
              result *= FreeSemiring::null();
            } else if (tmp_sr == SR::one()) {
              result *= FreeSemiring::one();
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
                result *= *value_iter;
              }
            }
          }
        }
      }

      return result;
    }

    bool operator<(const NonCommutativeMonomial &rhs) const {
      // FIXME
      return variables_ < rhs.variables_;
    }

    bool operator==(const NonCommutativeMonomial &rhs) const {
      // FIXME
      return variables_ == rhs.variables_;
    }

    Degree get_degree() const {
      return variables_.size();
    }

    // FIXME: modify or remove
    /*std::set<VarId> get_variables() const {
      std::set<VarId> set;
      for (auto var_degree : variables_) {
        set.insert(var_degree.first);
      }
      return set;
    }*/

    std::string string() const {
      std::stringstream ss;
      // FIXME
      ss << variables_;
      return std::move(ss.str());
    }
};

template <typename SR>
std::ostream& operator<<(std::ostream &out, const NonCommutativeMonomial<SR> &monomial);
