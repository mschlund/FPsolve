#pragma once

#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "key_wrapper.h"
#include "semilinear_util.h"
#include "sparse_vec.h"
#include "unique_linear_set.h"

template <typename Simplifier, typename V>
class LinearSet {
  public:
    LinearSet() : set_ptr_(builder_.New(SparseVec<V>{}, std::set< SparseVec<V> >{})) {}

    LinearSet(const LinearSet &lset) = default;
    LinearSet(LinearSet &&lset) = default;

    LinearSet(const SparseVec<V> &v) : set_ptr_(builder_.New(v, {})) {}
    LinearSet(SparseVec<V> &&v) : set_ptr_(builder_.New(std::move(v), {})) {}

    LinearSet(const SparseVec<V> &o, const std::set< SparseVec<V> > &vs)
        : set_ptr_(builder_.New(o, vs)) {}

    LinearSet(SparseVec<V> &&o, std::set< SparseVec<V> > &&vs)
        : set_ptr_(builder_.New(std::move(o), std::move(vs))) {}

    // FIXME: do we need this?
    // LinearSet& operator=(const LinearSet &s) = default;
    // LinearSet& operator=(LinearSet &&s) = default;

    ~LinearSet() = default;

    bool operator==(const LinearSet &rhs) const {
      return set_ptr_ == rhs.set_ptr_;
    }

    bool operator<(const LinearSet &rhs) const {
      return set_ptr_ < rhs.set_ptr_;
    }


    LinearSet operator+(const LinearSet &rhs) const {
      auto result_offset = GetOffset() + rhs.GetOffset();
      std::set< SparseVec<V> > result_generators;

      std::set_union(GetGenerators().begin(), GetGenerators().end(),
                     rhs.GetGenerators().begin(), rhs.GetGenerators().end(),
                     inserter(result_generators, result_generators.begin()));

      SimplifySet(simplifier_, result_generators);

      return LinearSet{builder_.New(std::move(result_offset),
                                    std::move(result_generators))};
    }

    std::size_t Hash() const {
      std::hash< UniqueLinearSetPtr<V> > h;
      return h(set_ptr_);
    }

    friend std::ostream& operator<<(std::ostream &out, const LinearSet lset) {
      out << *lset.set_ptr_;
      return out;
    }

    const SparseVec<V>& GetOffset() const { return set_ptr_->GetOffset(); }
    const std::set< SparseVec<V> >& GetGenerators() const {
      return set_ptr_->GetGenerators();
    }

  private:
    LinearSet(UniqueLinearSetPtr<V> s) : set_ptr_(s) {}

    UniqueLinearSetPtr<V> set_ptr_;

    /* TODO: Try to get rid of static... */
    static UniqueLinearSetBuilder<V> builder_;
    static Simplifier simplifier_;

    template <typename S2, typename S1, typename VV>
    friend LinearSet<S2, VV> ChangeLinearSimplifier(const LinearSet<S1, VV> &lset);

};

namespace std {

template<typename S, typename V>
struct hash< LinearSet<S, V> > {
  inline std::size_t operator()(const LinearSet<S, V> &set) const {
    return set.Hash();
  }

};

}  /* namespace std */

template <typename Simplifier, typename V>
UniqueLinearSetBuilder<V> LinearSet<Simplifier, V>::builder_;

template <typename Simplifier, typename V>
Simplifier LinearSet<Simplifier, V>::simplifier_;

template <typename S2, typename S1, typename V>
LinearSet<S2, V> ChangeLinearSimplifier(const LinearSet<S1, V> &lset) {
  return LinearSet<S2, V>{lset.set_ptr_};
}
