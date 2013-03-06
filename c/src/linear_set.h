#pragma once

#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "key_wrapper.h"
#include "semilinear_util.h"
#include "sparse_vec.h"
#include "unique_set.h"

template <typename Simplifier, typename V>
class LinearSet {
  public:
    LinearSet() : offset_(), generators_(builder_.New({})) {}

    LinearSet(const LinearSet &lset) = default;
    LinearSet(LinearSet &&lset) = default;

    LinearSet(const SparseVec<V> &v)
        : offset_(builder_.New(v, {})), generators_(builder_.New({})) {}

    LinearSet(SparseVec<V> &&v)
        : offset_(std::move(v)), generators_(builder_.New({})) {}

    LinearSet(const SparseVec<V> &o, const std::set< SparseVec<V> > &vs)
        : offset_(o), generators_(builder_.New(vs)) {}

    LinearSet(SparseVec<V> &&o, std::set< SparseVec<V> > &&vs)
        : offset_(std::move(o)), generators_(builder_.New(std::move(vs))) {}

    LinearSet& operator=(const LinearSet &s) = default;
    LinearSet& operator=(LinearSet &&s) = default;

    ~LinearSet() = default;

    bool operator==(const LinearSet &rhs) const {
      return offset_ == rhs.offset_ && generators_ == rhs.generators_;
    }

    bool operator<(const LinearSet &rhs) const {
      if (offset_ < rhs.offset_) {
        return true;
      } else if (offset_ == rhs.offset_) {
        return generators_ < rhs.generators_;
      }
      return false;
    }


    LinearSet operator+(const LinearSet &rhs) const {
      auto result_offset = offset_ + rhs.offset_;
      std::set< SparseVec<V> > result_generators;

      std::set_union(GetGenerators().begin(), GetGenerators().end(),
                     rhs.GetGenerators().begin(), rhs.GetGenerators().end(),
                     inserter(result_generators, result_generators.begin()));

      SimplifySet(simplifier_, result_generators);

      return LinearSet{std::move(result_offset),
                       builder_.New(std::move(result_generators))};
    }

    std::size_t Hash() const {
      std::size_t hash = 0;
      HashCombine(hash, offset_);
      HashCombine(hash, generators_);
      return hash;
    }

    friend std::ostream& operator<<(std::ostream &out, const LinearSet lset) {
      out << "<" << lset.offset_ << " : " << *lset.generators_ << ">";
      return out;
    }

    const SparseVec<V>& GetOffset() const { return offset_; }
    const std::set< SparseVec<V> >& GetGenerators() const {
      return generators_->GetSet();
    }

  private:
    LinearSet(SparseVec<V> &&o, UniqueSetPtr< SparseVec<V> > s)
        : offset_(std::move(o)), generators_(s) {}

    SparseVec<V> offset_;
    UniqueSetPtr< SparseVec<V> > generators_;

    /* TODO: Try to get rid of static... */
    static UniqueSetBuilder< SparseVec<V> > builder_;
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
UniqueSetBuilder< SparseVec<V> > LinearSet<Simplifier, V>::builder_;

template <typename Simplifier, typename V>
Simplifier LinearSet<Simplifier, V>::simplifier_;

template <typename S2, typename S1, typename V>
LinearSet<S2, V> ChangeLinearSimplifier(const LinearSet<S1, V> &lset) {
  return LinearSet<S2, V>{lset.set_ptr_};
}
