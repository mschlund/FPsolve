#pragma once

#include <initializer_list>
#include <set>

#include "linear_set.h"
#include "semiring.h"
#include "../datastructs/sparse_vec.h"
#include "../string_util.h"

class VarId;

template <
  typename Var = VarId,
  typename Value = Counter,
  DIVIDER_TEMPLATE_TYPE VecDivider = DummyDivider,
  VEC_SIMPL_TEMPLATE_TYPE VecSimpl = DummyVecSimplifier2,
  LIN_SIMPL_TEMPLATE_TYPE LinSimpl = DummyLinSimplifier>
class SemilinearSet;

/* Compatibility with the old implementation. */
typedef SemilinearSet<> SemilinSetExp;

/* Three typedefs for easy comparison between the impact of various
 * simplifications performs no simplification at all. */

typedef SemilinearSet<VarId, Counter, DummyDivider,
                      SparseVecSimplifier2, DummyLinSimplifier
                      > SemilinearSetV;

typedef SemilinearSet<VarId, Counter, DummyDivider,
                      SparseVecSimplifier2, LinearSetSimplifier
                      > SemilinearSetL;

/* DivSemilinearSet additionally divides the SparseVec by its gcd.  NOTE: this
 * is an over-approximation, the result might no longer be precise. */
typedef SemilinearSet< VarId, Counter, GcdDivider> DivSemilinearSet;


template <typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl,
          LIN_SIMPL_TEMPLATE_TYPE LinSimpl>
class SemilinearSet : public Semiring< SemilinearSet<Var, Value, VecDivider,
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
      if (IsZero()) {
        set_ = rhs.set_;
      } else if (!rhs.IsZero()) {
        // VecSet<LinearSetType> result;
        set_ = VecSetUnion(set_, rhs.set_);
        // std::set_union(set_.begin(), set_.end(),
        //                rhs.set_.begin(), rhs.set_.end(),
        //                std::inserter(result, result.begin()));

        // set_ = std::move(result);

        SimplifySet<LinSimplType>(set_);
      }

      return *this;
    }

    SemilinearSet& operator*=(const SemilinearSet &rhs) {
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
        // set_ = std::move(result);
        SimplifySet<LinSimplType>(set_);
      }

      return *this;
    }

    SemilinearSet star(const LinearSetType &lset) const {

      /* If we do not have any generators, i.e.,
       *   ls = w  (for some word w)
       * just return
       *   w*
       * instead of 1 + ww*. */
      if (lset.GetGenerators().empty()) {
        VecSet<GeneratorType> result_gens;
        /* If w is not the one-element, move w to the generators. */
        if (lset.GetOffset() != OffsetType{}) {
          result_gens.emplace_back(GeneratorType{lset.GetOffset()});
        }
        return SemilinearSet{ LinearSetType{
                                OffsetType{}, std::move(result_gens)} };
      }

      /* Star of a linear set is a semilinear set:
       * (w_0.w_1*.w_2*...w_n*)* = 1 + (w_0.w_0*.w_1*.w_2*...w_n*) */

      VecSet<GeneratorType> result_gens =
        VecSetUnion(lset.GetGenerators(),
                    VecSet<GeneratorType>{ GeneratorType{lset.GetOffset()} });

      SemilinearSet result{
        LinearSetType{lset.GetOffset(), std::move(result_gens)} };

      /* Insert one.  We're inlining the definition for efficiency. */
      // result.set_.insert(LinearSetType{});

      return one() + result;
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
