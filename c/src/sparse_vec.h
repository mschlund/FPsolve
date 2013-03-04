#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

#include "hash.h"
#include "unique_vector_map.h"

typedef std::uint_fast32_t Counter;

/*
 * Sparse vector representing the mapping from variables to counters.  We never
 * create a vector that has already be created, therefore operations such as
 * equality are really just pointer equality (i.e., they are very efficient).
 * Moreover we never modify a vector, so copy constructor, assignment operator,
 * etc. are really only copying a pointer (which is again very efficient).
 */
template <typename Var, typename Value = Counter>
class SparseVec {
  typedef UniqueVMap<Var, Value> UniqueVMap_;
  typedef UniqueVMapPtr<Var, Value> UniqueVMapPtr_;
  typedef UniqueVMapBuilder<Var, Value> UniqueVMapBuilder_;

  public:
    SparseVec() {
      vmap_ = builder_.New({});
    }

    SparseVec(const SparseVec &v) = default;
    SparseVec(SparseVec &&v) = default;

    SparseVec& operator=(const SparseVec &v) = default;
    SparseVec& operator=(SparseVec &&v) = default;

    SparseVec(std::vector< std::pair<Var, Value> > &&vector) {
      vmap_ = builder_.New(std::move(vector));
    }

    SparseVec(std::initializer_list< std::pair<Var, Value> > list)
        : SparseVec(std::vector< std::pair<Var, Value> >{list}) {}

    SparseVec(const Var &v, Value c)
      : SparseVec({ std::make_pair(v, c) }) {}

    bool operator==(const SparseVec &rhs) const {
      assert(Sanity() && rhs.Sanity());
      if (vmap_ == nullptr || rhs.vmap_ == nullptr) {
        return false;
      }
      return vmap_ == rhs.vmap_;
    }

    bool operator!=(const SparseVec &rhs) const {
      return !(*this == rhs);
    }

    /* Note that this is _not_ the lexicographical ordering of the contents of
     * the vectors! */
    bool operator<(const SparseVec &rhs) const {
      assert(Sanity() && rhs.Sanity());
      if (vmap_ == nullptr || rhs.vmap_ == nullptr) {
        return false;
      }
      return vmap_ < rhs.vmap_;
    }

    SparseVec operator+(const SparseVec &rhs) const {
      assert(Sanity() && rhs.Sanity());
      if (vmap_ == nullptr || rhs.vmap_ == nullptr) {
        return SparseVec{nullptr};
      }
      return builder_.NewSum(*vmap_, *rhs.vmap_);
    }

    /*
     * IMPORTANT: This might create an _invalid_ SparseVec.  Whenever you use
     * operator- you should check that SparseVec.IsValid()!
     */
    SparseVec operator-(const SparseVec &rhs) const {
      assert(Sanity() && rhs.Sanity());
      if (vmap_ == nullptr || rhs.vmap_ == nullptr) {
        return SparseVec{nullptr};
      }
      return builder_.NewDiff(*vmap_, *rhs.vmap_);
    }

    bool IsZero() const {
      assert(Sanity());
      return vmap_ == nullptr ? false : vmap_->empty();
    }

    bool IsValid() const {
      return vmap_ != nullptr;
    }

    bool Divides(const SparseVec &rhs) const {
      assert(Sanity() && rhs.Sanity());

      /* Propagate invalid SparseVec. */
      if (vmap_ == nullptr || rhs.vmap_ == nullptr) {
        return false;
      }

      if (vmap_->size() != rhs.vmap_->size()) {
        return false;
      }

      Counter k = 0;

      for (auto lhs_iter = vmap_->begin(), rhs_iter = rhs.vmap_->begin();
           lhs_iter != vmap_->end();
           ++lhs_iter, ++rhs_iter) {

        assert(rhs_iter != rhs.vmap_->end());

        if (/* If the domains are different. */
            lhs_iter->first != rhs_iter->first ||
            /* Or we have a multiple and it doesn't work here. */
            (k != 0 && rhs_iter->second != k * lhs_iter->second) ||
            /* Or we can't divide. */
            rhs_iter->second % lhs_iter->second != 0) {
          return false;
        }

        /* Otherwise we have a candidate for the multiple. */
        k = rhs_iter->second / lhs_iter->second;

      }
      return true;
    }


    friend std::ostream& operator<<(std::ostream &out, const SparseVec &svector) {
      assert(svector.Sanity());
      out << *svector.vmap_;
      return out;
    }

    std::size_t Hash() const {
      std::hash<UniqueVMapPtr_> h;
      return h(vmap_);
    }

  private:
    SparseVec(UniqueVMapPtr_ v) : vmap_(v) {}

    bool Sanity() const {
      if (vmap_ == nullptr) {
        return false;
      }
      return true;
    }

    UniqueVMapPtr_ vmap_;
    static UniqueVMapBuilder_ builder_;
};

template <typename Var, typename Value>
UniqueVMapBuilder<Var, Value> SparseVec<Var, Value>::builder_;

namespace std {

template<typename Var, typename Value>
struct hash< SparseVec<Var, Value> > {
  inline std::size_t operator()(const SparseVec<Var, Value> &vec) const {
    return vec.Hash();
  }

};

}  /* namespace std */


class DummySimplifier {
  public:
    bool IsActive() const { return false; }

    template <typename Elem>
    bool IsCovered(const Elem &lhs, const std::set<Elem> &rhs_set) {
      return false;
    }
};


template <typename Var, typename Value = Counter>
class NaiveSimplifier {
  public:
    bool IsActive() const { return true; }

    bool IsCovered(const SparseVec<Var, Value> &lhs,
        const std::set< SparseVec<Var, Value> > &rhs_set) {
      if (rhs_set.count(lhs) > 0) {
        return true;
      }
      for (const auto &rhs_gen : rhs_set) {
        if (rhs_gen.Divides(lhs)) {
          return true;
        }
      }
      return false;
    }
};

/*
 * TODO: This simplifier assumes that we can use any of the elements of the set
 * arbitrary many times.  We should also have one that uses every element of the
 * set only once...
 */

template <typename Var, typename Value = Counter>
class SparseVecSimplifier : public NaiveSimplifier<Var, Value> {
  public:
    bool IsActive() const { return true; }

    bool IsCovered(const SparseVec<Var, Value> &lhs,
        const std::set< SparseVec<Var, Value> > &rhs_set) {
      /* Check the cheap and naive simplifier. */
      if (NaiveSimplifier<Var, Value>::IsCovered(lhs, rhs_set)) {
        return true;
      }
      std::unordered_set< SparseVec<Var, Value> > failed;
      return IsCovered_(lhs, rhs_set, failed);
    }
  private:
    /* Dynamic programming/memoization */
    bool IsCovered_(const SparseVec<Var, Value> &lhs,
        const std::set< SparseVec<Var, Value> > &rhs_set,
        std::unordered_set< SparseVec<Var, Value> > &failed) {
      if (0 < failed.count(lhs)) {
        return false;
      }
      /* Should we go from the back? */
      for (auto &rhs : rhs_set) {
        if (rhs.IsZero()) {
          assert(false);
          continue;
        }
        auto new_lhs = lhs - rhs;
        if (new_lhs.IsValid() && (new_lhs.IsZero() ||
                                  IsCovered_(new_lhs, rhs_set, failed))) {
          return true;
        }
      }
      failed.insert(lhs);
      return false;
    }
};
