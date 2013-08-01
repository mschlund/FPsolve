#pragma once

#include <initializer_list>
#include <set>

#include "linear_set.h"
#include "semiring.h"
#include "../datastructs/sparse_vec.h"
#include "../string_util.h"
#include "../utils/profiling-macros.h"

class VarId;

template <
  typename Var = VarId,
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
typedef SemilinearSet< VarId, Counter, GcdDivider> DivSemilinearSet;

template <typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl,
          LIN_SIMPL_TEMPLATE_TYPE LinSimpl>
class SemilinearSet : public StarableSemiring< SemilinearSet<Var, Value, VecDivider,
                                                     VecSimpl, LinSimpl>,
                                       Commutativity::Commutative,
                                       Idempotence::Idempotent > {
  public:
    typedef SparseVec<Var, Value, DummyDivider> OffsetType;
    typedef SparseVec<Var, Value, VecDivider> GeneratorType;
    typedef LinearSet<Var, Value, VecDivider, VecSimpl> LinearSetType;
    typedef LinSimpl<Var, Value, VecDivider, VecSimpl> LinSimplType;

    SemilinearSet() = default;
    SemilinearSet(std::initializer_list<LinearSetType> list)
        : set_(list) {}
    SemilinearSet(const SemilinearSet &slset) = default;
    SemilinearSet(SemilinearSet &&slset) = default;

    SemilinearSet(const LinearSetType &lset) : set_({lset}) {}
    SemilinearSet(LinearSetType &&lset) : set_({std::move(lset)}) {}

    SemilinearSet(const Var &v, Counter c)
      : set_({ LinearSetType{ OffsetType{v, c} } }) {}
    SemilinearSet(const Var &v) : SemilinearSet(v, 1) {}

    template <DIVIDER_TEMPLATE_TYPE OldVecDivider,
              VEC_SIMPL_TEMPLATE_TYPE OldVecSimpl,
              LIN_SIMPL_TEMPLATE_TYPE OldLinSimpl>
    SemilinearSet(const SemilinearSet<Var, Value, OldVecDivider,
                                      OldVecSimpl, OldLinSimpl> &slset) {
      for (const auto &lset : slset) {
        /* Elements of slset are already sorted, so it's safe to use
         * emplace_back. */
        set_.emplace_back(LinearSetType{lset});
      }
    }

    ~SemilinearSet() = default;

    bool operator==(const SemilinearSet &rhs) const {
      return set_ == rhs.set_;
    }

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
            [](const OffsetType &o) { return !o.IsZero(); });

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

    typedef typename VecSet<LinearSetType>::const_iterator const_iterator;

    const_iterator begin() const { return set_.begin(); }
    const_iterator end() const { return set_.end(); }

  private:
    SemilinearSet(VecSet<LinearSetType> &&s) : set_(s) {}

    VecSet<LinearSetType> set_;
};
