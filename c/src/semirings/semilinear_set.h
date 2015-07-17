#pragma once

#include <initializer_list>
#include <set>

#include "linear_set.h"
#include "semiring.h"
#include "../datastructs/sparse_vec.h"
#include "../utils/string_util.h"
#include "../utils/profiling-macros.h"
#ifdef USE_GENEPI
#include "semilinSetNdd.h"
#endif

#include <boost/algorithm/string.hpp>

class VarId;

template <
  typename VarType = VarId,
  typename Value = Counter,
  DIVIDER_TEMPLATE_TYPE VecDivider = DummyDivider,
  VEC_SIMPL_TEMPLATE_TYPE VecSimpl = DummyVecSimplifier,
  LIN_SIMPL_TEMPLATE_TYPE LinSimpl = DummyLinSimplifier>
class SemilinearSet;

/* Compatibility with the old implementation. */
typedef SemilinearSet<> SemilinSetExp;

/* Three typedefs for easy comparison between the impact of various
 * simplifications performs no simplification at all. */

typedef SemilinearSet<VarId, Counter, DummyDivider,
                      SparseVecSimplifier, DummyLinSimplifier
                      > SemilinearSetV;

typedef SemilinearSet<VarId, Counter, DummyDivider,
                      SparseVecSimplifier, LinearSetSimplifier
                      > SemilinearSetL;

/* DivSemilinearSet additionally divides the SparseVec by its gcd.  NOTE: this
 * is an over-approximation, the result might no longer be precise. */
typedef SemilinearSet< VarId, Counter, GcdDivider, SparseVecSimplifier, LinearSetSimplifier> DivSemilinearSet;


template <typename VarType,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl,
          LIN_SIMPL_TEMPLATE_TYPE LinSimpl>
class SemilinearSet : public StarableSemiring< SemilinearSet<VarType, Value, VecDivider,
                                                     VecSimpl, LinSimpl>,
                                       Commutativity::Commutative,
                                       Idempotence::Idempotent > {
  public:
    typedef SparseVec<VarType, Value, DummyDivider> OffsetType;
    typedef SparseVec<VarType, Value, VecDivider> GeneratorType;
    typedef LinearSet<VarType, Value, VecDivider, VecSimpl> LinearSetType;
    typedef LinSimpl<VarType, Value, VecDivider, VecSimpl> LinSimplType;

    SemilinearSet() = default;
    SemilinearSet(std::initializer_list<LinearSetType> list)
        : set_(list) {}
    SemilinearSet(const SemilinearSet &slset) = default;
    SemilinearSet(SemilinearSet &&slset) = default;

    SemilinearSet(const LinearSetType &lset) : set_({lset}) {}
    SemilinearSet(LinearSetType &&lset) : set_({std::move(lset)}) {}

    SemilinearSet(const VarType &v, Counter c)
      : set_({ LinearSetType{ OffsetType{v, c} } }) {}
    SemilinearSet(const VarType &v) : SemilinearSet(v, 1) {}

    template <DIVIDER_TEMPLATE_TYPE OldVecDivider,
              VEC_SIMPL_TEMPLATE_TYPE OldVecSimpl,
              LIN_SIMPL_TEMPLATE_TYPE OldLinSimpl>
    SemilinearSet(const SemilinearSet<VarType, Value, OldVecDivider,
                                      OldVecSimpl, OldLinSimpl> &slset) {
      for (const auto &lset : slset) {
        /* Elements of slset are already sorted, so it's safe to use
         * emplace_back. */
        set_.emplace_back(LinearSetType{lset});
      }
    }

    // parse a description of the form e.g. "<a:3,b:2>" or of the form "a"
    SemilinearSet(const std::string &str_val) {
      assert(str_val.length() > 0);

      if(str_val.compare("<>") == 0 || str_val.compare("()") == 0) {
        set_.emplace_back(LinearSetType{ OffsetType() } );
      }
      else if(str_val.front() == '<') {
        //split by "," , then every component by ":"
        std::vector<std::string> elems;
        std::vector<std::pair<VarType,Value> > sparse_values;
        boost::split(elems, str_val, boost::is_any_of(",<>"), boost::algorithm::token_compress_on);

        for (std::string s : elems) {
          if(s.empty()) {
            continue;
          }
          //std::cout << "Token: " << s << std::endl;
          std::vector<std::string> var_val;
          boost::split(var_val, s, boost::is_any_of(":"), boost::algorithm::token_compress_on);

          if(var_val.size() != 2) {
            std::cerr << "Bad Input (semilinear set): \"" << str_val << "\"" << std::endl;
          }
          else {
            std::istringstream i(var_val[1]);
            Counter c;
            if (!(i >> c))
            {
              std::cerr << "ERROR: Bad string value (" << var_val[1] << ") for semilinear-set constructor (element count) --defaulting to 0!"<< std::endl;
              c = 0;
            }
            VarType v = Var::GetVarId(var_val[0]);
            if(c != 0) {
              sparse_values.push_back(std::make_pair(v,c));
            }
          }
        }
        set_.emplace_back(LinearSetType{ OffsetType(std::move(sparse_values)) } );
      }
      else {
        set_.emplace_back(LinearSetType{ OffsetType{Var::GetVarId(str_val), 1} } );
      }
    }

    ~SemilinearSet() = default;

// if we use genepi then we can check equivalence using NDDs (sound+complete)
#ifdef USE_GENEPI
    bool operator==(const SemilinearSet &rhs) const {
      if (this->IsZero() && rhs.IsZero() || this->IsOne() && rhs.IsOne())
        return true;

      if (this->IsZero() && rhs.IsOne() || this->IsOne() && rhs.IsZero())
        return false;

      auto V1 = this->getVariables();
      auto V2 = rhs.getVariables();
      if (V1 != V2)
        return false;

      int k = V1.size();

      //std::cout << "numvars: " << k << std::endl;

      SemilinSetNdd::solver_init(k);
      //std::cout << *this << " ?== " << rhs << std::endl;
      //std::cout << "testing for eq in slsetNDD" << std::endl;
      bool eq = (SemilinSetNdd(this->set_) == SemilinSetNdd(rhs.set_));
      //std::cout << "eq solved: " << *this << " ?==" << rhs << ": " << eq << std::endl;
      SemilinSetNdd::solver_dealloc();
      return eq;
    }
#else
    bool operator==(const SemilinearSet &rhs) const {
      return set_ == rhs.set_;
    }
#endif

    static SemilinearSet null() {
      return SemilinearSet{};
    }

    static SemilinearSet one() {
      return SemilinearSet{LinearSetType{}};
    }

    SemilinearSet& operator=(const SemilinearSet &slset) = default;
    SemilinearSet& operator=(SemilinearSet &&slset) = default;

    SemilinearSet& operator+=(const SemilinearSet &rhs) {
      OPADD;
      if (IsZero()) {
        set_ = rhs.set_;
      } else if (!rhs.IsZero()) {
        set_ = VecSetUnion(set_, rhs.set_);
        SimplifySet<LinSimplType>(set_);
      }

      return *this;
    }

    SemilinearSet& operator*=(const SemilinearSet &rhs) {
      OPMULT;
      if (IsZero() || rhs.IsZero()) {
        *this = null();
      } else if (IsOne()) {
        set_ = rhs.set_;
      } else if (!rhs.IsOne()) {
        std::set<LinearSetType> tmp_result;
        for(auto &lin_set_rhs : rhs.set_) {
          for(auto &lin_set_lhs : set_) {
            tmp_result.insert(lin_set_lhs + lin_set_rhs);
          }
        }
        set_.clear();
        for (auto &lset : tmp_result) {
          set_.emplace_back(std::move(lset));
        }
        SimplifySet<LinSimplType>(set_);
      }

      return *this;
    }

    SemilinearSet star(const LinearSetType &lset) const {
      /* Remember to check if an offset is a zero vector and don't put that to
       * generators... */

      /* If we do not have any generators, i.e.,
       *   ls = w  (for some word w)
       * just return
       *   w*
       * instead of 1 + ww*. */
      if (lset.GetGenerators().empty()) {
        VecSet<GeneratorType> result_gens;
        /* If w is not the one-element, move w to the generators. */
        if (lset.GetOffset() != OffsetType{} && !lset.GetOffset().IsZero()) {
          result_gens.emplace_back(GeneratorType{lset.GetOffset()});
        }
        return SemilinearSet{ LinearSetType{
                                OffsetType{}, std::move(result_gens)} };
      }

      /* Star of a linear set is a semilinear set:
       * (w_0.w_1*.w_2*...w_n*)* = 1 + (w_0.w_0*.w_1*.w_2*...w_n*) */

      VecSet<GeneratorType> result_gens =
        VecSetUnionWith(
            lset.GetGenerators(),
            VecSet<GeneratorType>{ GeneratorType{lset.GetOffset()} },
            [](const GeneratorType &o) { return !o.IsZero(); });

      /* We need to add one(), but to avoid additional simplification step, we
       * just "inline" it... */
      SemilinearSet result;
      LinearSetType tmp_lin{lset.GetOffset(), std::move(result_gens)};
      LinearSetType tmp_one;
      if (tmp_lin < tmp_one) {
        result.set_.emplace_back(std::move(tmp_lin));
        result.set_.emplace_back(std::move(tmp_one));
      } else {
        result.set_.emplace_back(std::move(tmp_one));
        result.set_.emplace_back(std::move(tmp_lin));
      }

      return result;
    }

    SemilinearSet star() const {
      OPSTAR;
      SemilinearSet result = one();
      for (auto &ls : set_) {
        result *= star(ls);
      }

      return result;
    }

    bool IsZero() const {
      return set_.size() == 0;
    }

    bool IsOne() const {
      return set_.size() == 1 && set_.begin()->IsZero();
    }

    std::string string() const {
      std::stringstream sout;
      sout << "{ " << std::endl
           << ToStringSorted(set_, "\n")
           << "}" << std::endl;
      return std::move(sout.str());
    }

    // note that this is potentially expensive ..
    std::set<VarType> getVariables() const {
      std::set<VarType> res;
      for (const auto &linset : set_) {
        auto ovars = linset.GetOffset().getVariables();
        res.insert(ovars.begin(), ovars.end());
        for (const auto &g : linset.GetGenerators()) {
          auto gvars = g.getVariables();
          res.insert(gvars.begin(),gvars.end());
        }
      }
      return res;
    }

    typedef typename VecSet<LinearSetType>::const_iterator const_iterator;

    const_iterator begin() const { return set_.begin(); }
    const_iterator end() const { return set_.end(); }


  private:
    SemilinearSet(VecSet<LinearSetType> &&s) : set_(s) {}
    VecSet<LinearSetType> set_;
};
