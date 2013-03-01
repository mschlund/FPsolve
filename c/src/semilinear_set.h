#pragma once

#include <initializer_list>
#include <iostream>
#include <set>

#include "linear_set.h"
#include "semiring.h"
#include "sparse_vec.h"
#include "var.h"

template <typename SemiSimpl, typename LinSimpl, typename V>
class SemilinearSet : public Semiring<
                               SemilinearSet<SemiSimpl, LinSimpl, V> > {
  public:
    SemilinearSet() = default;
    SemilinearSet(std::initializer_list< LinearSet<LinSimpl, V> > list)
        : set_(list) {}
    SemilinearSet(const SemilinearSet &slset) = default;
    SemilinearSet(SemilinearSet &&slset) = default;

    SemilinearSet(const LinearSet<LinSimpl, V> &lset) : set_({lset}) {}
    SemilinearSet(LinearSet<LinSimpl, V> &&lset) : set_({std::move(lset)}) {}

    SemilinearSet(const V &v, Counter c)
      : set_({ LinearSet<LinSimpl, V>{ SparseVec<V>{v, c} } }) {}
    SemilinearSet(const V &v) : SemilinearSet(v, 1) {}

    ~SemilinearSet() = default;

    bool operator==(const SemilinearSet &rhs) const {
      return set_ == rhs.set_;
    }

    static SemilinearSet null() {
      return SemilinearSet{};
    }

    static SemilinearSet one() {
      return SemilinearSet{LinearSet<LinSimpl, V>{}};
    }

    SemilinearSet& operator=(const SemilinearSet &slset) = default;
    SemilinearSet& operator=(SemilinearSet &&slset) = default;

    SemilinearSet& operator+=(const SemilinearSet &rhs) {
      std::set< LinearSet<LinSimpl, V> > result;
      std::set_union(set_.begin(), set_.end(),
                     rhs.set_.begin(), rhs.set_.end(),
                     std::inserter(result, result.begin()));

      set_ = std::move(result);

      Simplify();

      return *this;
    }

    SemilinearSet& operator*=(const SemilinearSet &rhs) {
      // std::cout << "-> SemilinearSet::operator*" << std::endl;
      // std::cout << *this << std::endl;
      // std::cout << rhs << std::endl;
      std::set< LinearSet<LinSimpl, V> > result;
      for(auto &lin_set_rhs : rhs.set_) {
        for(auto &lin_set_lhs : set_) {
          result.insert(lin_set_lhs + lin_set_rhs);
        }
      }

      set_ = std::move(result);

      Simplify();

      // std::cout << *this << std::endl;
      // std::cout << "<- SemilinearSet::operator*" << std::endl;
      return *this;
    }

    SemilinearSet star(const LinearSet<LinSimpl, V> &lset) const {

      /* If we do not have any generators, i.e.,
       *   ls = w  (for some word w)
       * just return
       *   w*
       * instead of 1 + ww*. */
      if (lset.GetGenerators().empty()) {
        std::set< SparseVec<V> > result_gens;
        /* If w is not the one-element, move w to the generators. */
        if (lset.GetOffset() != SparseVec<V>{}) {
          result_gens.insert(lset.GetOffset());
        }
        return SemilinearSet{ LinearSet<LinSimpl, V>{
                                SparseVec<V>{}, std::move(result_gens)} };
      }

      /* Star of a linear set is a semilinear set:
       * (w_0.w_1*.w_2*...w_n*)* = 1 + (w_0.w_0*.w_1*.w_2*...w_n*) */

      std::set< SparseVec<V> > result_gens = lset.GetGenerators();
      result_gens.insert(lset.GetOffset());

      SemilinearSet result{ LinearSet<LinSimpl, V>{
                              lset.GetOffset(), std::move(result_gens)} };

      /* Insert one.  We're inlining the definition for efficiency. */
      result.set_.insert(LinearSet<LinSimpl, V>{});

      return result;
    }

    SemilinearSet star() const {
      SemilinearSet result = one();
      for (auto &ls : set_) {
        result *= star(ls);
      }

      return result;
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

    template <typename F>
    void Iterate(F fun) const {
      for (const auto &lset : set_) {
        fun(lset);
      }
    }

    static const bool is_idempotent = true;
    static const bool is_commutative = true;

  private:
    SemilinearSet(std::set< LinearSet<LinSimpl, V> > &&s) : set_(s) {}

    void Simplify() {
      SimplifySet(simplifier_, set_);
    }

    std::set< LinearSet<LinSimpl, V> > set_;
    SemiSimpl simplifier_;

    template <typename S21, typename S22, typename S11, typename S12, typename VV>
    friend SemilinearSet<S21, S22, VV> ChangeSimplifiers(
      const SemilinearSet<S11, S12, VV> &slset);
};

template <typename S21, typename S22, typename S11, typename S12, typename V>
SemilinearSet<S21, S22, V> ChangeSimplifiers(
    const SemilinearSet<S11, S12, V> &slset) {
  std::set< LinearSet<S22, V> > result_set;
  for (const auto &x : slset.set_) {
    // FIXME: GCC 4.7 doesn't have emplace.
    result_set.insert(ChangeLinearSimplifier<S22>(x));
  }
  return SemilinearSet<S21, S22, V>{std::move(result_set)};
}

template <typename LinearSimpl, typename Var>
class SemilinearSimplifier {
  public:
    bool IsActive() const { return true; }
    bool IsCovered(LinearSet<LinearSimpl, Var> &lset,
                   const std::set< LinearSet<LinearSimpl, Var> > &rhs_lsets) {
      for (auto &rhs_lset : rhs_lsets) {
        auto new_offset = lset.GetOffset() - rhs_lset.GetOffset();
        if (!new_offset.IsValid() ||
            !linear_simpl_.IsCovered(new_offset, rhs_lset.GetGenerators())) {
          continue;
        }
        for (auto &lhs_gen : lset.GetGenerators()) {
          if (!linear_simpl_.IsCovered(lhs_gen, rhs_lset.GetGenerators())) {
            continue;
          }
        }
        return true;
      }
      return false;
    }
  private:
    LinearSimpl linear_simpl_;
};

template <typename LinearSimpl, typename Var>
class SemilinearSubsetSimplifier {
  public:
    bool IsActive() const { return true; }
    bool IsCovered(LinearSet<LinearSimpl, Var> &lset,
                   const std::set< LinearSet<LinearSimpl, Var> > &rhs_lsets) {
      for (auto &rhs_lset : rhs_lsets) {
        auto new_offset = lset.GetOffset() - rhs_lset.GetOffset();
        if (new_offset.IsValid() &&
            linear_simpl_.IsCovered(new_offset, rhs_lset.GetGenerators()) &&
            std::includes(rhs_lset.GetGenerators().begin(),
                          rhs_lset.GetGenerators().end(),
                          lset.GetGenerators().begin(),
                          lset.GetGenerators().end())) {
          return true;
        }
      }
      return false;
    }
  private:
    LinearSimpl linear_simpl_;
};

/* Compatibility with old implementation. */
typedef SemilinearSet<SemilinearSimplifier< SmartSimplifier<VarPtr>, VarPtr >,
                      SmartSimplifier<VarPtr>,
                      VarPtr> SemilinSetExp;
// typedef SemilinearSet<DummySimplifier, SmartSimplifier<VarPtr>, VarPtr> SemilinSetExp;
// typedef SemilinearSet<DummySimplifier, NaiveSimplifier<VarPtr>, VarPtr> SemilinSetExp;
// typedef SemilinearSet<DummySimplifier, DummySimplifier, VarPtr> SemilinSetExp;
