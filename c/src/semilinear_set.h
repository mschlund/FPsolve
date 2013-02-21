#pragma once

#include <initializer_list>
#include <set>

#include "linear_set.h"
#include "semiring.h"
#include "sparse_vec.h"
#include "var.h"

template <typename Simplifier1, typename Simplifier2, typename V>
class SemilinearSet : public Semiring<
                               SemilinearSet<Simplifier1, Simplifier2, V> > {
  public:
    SemilinearSet() = default;
    SemilinearSet(std::initializer_list< LinearSet<Simplifier2, V> > list)
        : sets_(list) {}
    SemilinearSet(const SemilinearSet &slset) = default;
    SemilinearSet(SemilinearSet &&slset) = default;

    SemilinearSet(const LinearSet<Simplifier2, V> &lset) : sets_({lset}) {}
    SemilinearSet(LinearSet<Simplifier2, V> &&lset) : sets_({std::move(lset)}) {}

    SemilinearSet(const V &v, Counter c)
      : sets_({ LinearSet<Simplifier2, V>{ SparseVec<V>{v, c} } }) {}
    SemilinearSet(const V &v) : SemilinearSet(v, 1) {}

    ~SemilinearSet() = default;

    bool operator==(const SemilinearSet &rhs) const {
      return sets_ == rhs.sets_;
    }

    static SemilinearSet null() {
      return SemilinearSet{};
    }

    static SemilinearSet one() {
      return SemilinearSet{LinearSet<Simplifier2, V>{}};
    }

    SemilinearSet& operator=(const SemilinearSet &slset) = default;
    SemilinearSet& operator=(SemilinearSet &&slset) = default;

    SemilinearSet& operator+=(const SemilinearSet &rhs) {
      std::set< LinearSet<Simplifier2, V> > result;
      std::set_union(sets_.begin(), sets_.end(),
                     rhs.sets_.begin(), rhs.sets_.end(),
                     std::inserter(result, result.begin()));
      sets_ = std::move(result);
      // FIXME: simplification
      return *this;
    }

    SemilinearSet& operator*=(const SemilinearSet &rhs) {
      std::set< LinearSet<Simplifier2, V> > result;
      for(auto &lin_set_rhs : rhs.sets_) {
        for(auto &lin_set_lhs : sets_) {
          result.insert(lin_set_lhs + lin_set_rhs);
        }
      }
      sets_ = std::move(result);

      return *this;
    }

    SemilinearSet star(const LinearSet<Simplifier2, V> &lset) const {

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
        return SemilinearSet{ LinearSet<Simplifier2, V>{
                                SparseVec<V>{}, std::move(result_gens)} };
      }

      /* Star of a linear set is a semilinear set:
       * (w_0.w_1*.w_2*...w_n*)* = 1 + (w_0.w_0*.w_1*.w_2*...w_n*) */

      std::set< SparseVec<V> > result_gens = lset.GetGenerators();
      result_gens.insert(lset.GetOffset());

      SemilinearSet result{ LinearSet<Simplifier2, V>{
                              lset.GetOffset(), std::move(result_gens)} };

      /* Insert one.  We're inlining the definition for efficiency. */
      result.sets_.insert(LinearSet<Simplifier2, V>{});

      return result;
    }

    SemilinearSet star() const {
      SemilinearSet result = one();
      for (auto &ls : sets_) {
        result *= star(ls);
      }
      return result;
    }

    std::string string() const {
      std::stringstream sout;
      sout << "{ ";
      for (const auto &ls : sets_) {
        sout << ls << " ";
      }
      sout << "}";
      return std::move(sout.str());
    }

  private:
    std::set< LinearSet<Simplifier2, V> > sets_;
};

typedef SemilinearSet<DummySimplifier, DummySimplifier, VarPtr>
        SemilinearSetSimple;
