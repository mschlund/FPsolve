#pragma once

#include <initializer_list>
#include <set>

#include "linear_set.h"
#include "semiring.h"
#include "sparse_vec.h"

class VarId;

template <typename Var = VarId,
          typename Value = Counter,
          typename VecDivider = DummyDivider,
          typename VecSimpl = SparseVecSimplifier<Var, Value, VecDivider>,
          typename LinSimpl = LinearSetSimplifier< Var, Value, VecDivider, VecSimpl> >
class SemilinearSet;

/* Compatibility with the old implementation. */
typedef SemilinearSet<> SemilinSetExp;

/* SimpleLinearSet performs no simplification at all. */
typedef SemilinearSet<VarId, Counter, DummyDivider,
                      DummySimplifier, DummySimplifier
                      > SimpleSemilinearSet;

/* DivSemilinearSet additionally divides the SparseVec by its gcd.  NOTE: this
 * is an over-approximation, the result might no longer be precise. */
typedef SemilinearSet< VarId, Counter,
                       GcdDivider<VarId, Counter> > DivSemilinearSet;

template <typename Var,
          typename Value,
          typename VecDivider,
          typename VecSimpl,
          typename LinSimpl>
class SemilinearSet : public Semiring< SemilinearSet<Var, Value, VecDivider,
                                                     VecSimpl, LinSimpl> > {
  public:
    typedef SparseVec<Var, Value, DummyDivider> OffsetType;
    typedef SparseVec<Var, Value, VecDivider> GeneratorType;
    typedef LinearSet<Var, Value, VecDivider, VecSimpl> LinearSetType;

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

    template <typename OldVecDivider, typename OldVecSimpl, typename OldLinSimpl>
    SemilinearSet(const SemilinearSet<Var, Value, OldVecDivider,
                                      OldVecSimpl, OldLinSimpl> &slset) {
      for (const auto &lset : slset) {
        set_.insert(LinearSetType{lset});
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
        std::set<LinearSetType> result;
        std::set_union(set_.begin(), set_.end(),
                       rhs.set_.begin(), rhs.set_.end(),
                       std::inserter(result, result.begin()));

        set_ = std::move(result);

        SimplifySet<LinSimpl>(set_);
      }

      return *this;
    }

    SemilinearSet& operator*=(const SemilinearSet &rhs) {
      if (IsZero() || rhs.IsZero()) {
        *this = null();
      } else if (IsOne()) {
        set_ = rhs.set_;
      } else if (!rhs.IsOne()) {
        std::set< LinearSetType > result;
        for(auto &lin_set_rhs : rhs.set_) {
          for(auto &lin_set_lhs : set_) {
            result.insert(lin_set_lhs + lin_set_rhs);
          }
        }
        set_ = std::move(result);
        SimplifySet<LinSimpl>(set_);
      }

      return *this;
    }

    // for any n in Nat: n*S = S if n>0 and 0 otherwise (the counting-SR is idempotent)
    SemilinearSet& operator*=(const std::uint_fast16_t &rhs) {
      if(rhs == 0){
        *this = null();
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
        std::set<GeneratorType> result_gens;
        /* If w is not the one-element, move w to the generators. */
        if (lset.GetOffset() != OffsetType{}) {
          result_gens.insert(GeneratorType{lset.GetOffset()});
        }
        return SemilinearSet{ LinearSetType{
                                OffsetType{}, std::move(result_gens)} };
      }

      /* Star of a linear set is a semilinear set:
       * (w_0.w_1*.w_2*...w_n*)* = 1 + (w_0.w_0*.w_1*.w_2*...w_n*) */

      std::set<GeneratorType> result_gens = lset.GetGenerators();
      result_gens.insert(GeneratorType{lset.GetOffset()});

      SemilinearSet result{ LinearSetType{
                              lset.GetOffset(), std::move(result_gens)} };

      /* Insert one.  We're inlining the definition for efficiency. */
      result.set_.insert(LinearSetType{});

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
      sout << "{ " << std::endl;
      for (const auto &ls : set_) {
        sout << ls << std::endl;
      }
      sout << "}" << std::endl;
      return std::move(sout.str());
    }

    typedef typename std::set<LinearSetType>::const_iterator const_iterator;

    const_iterator begin() const { return set_.begin(); }
    const_iterator end() const { return set_.end(); }

    static const bool is_idempotent = true;
    static const bool is_commutative = true;

  private:
    SemilinearSet(std::set<LinearSetType> &&s) : set_(s) {}

    std::set<LinearSetType> set_;
};
