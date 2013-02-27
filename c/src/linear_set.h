#pragma once

#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "key_wrapper.h"
#include "sparse_vec.h"
#include "semilinear_util.h"

// FIXME: this should be a class
template <typename V>
using OffsetGenerators = std::pair< SparseVec<V>, std::set< SparseVec<V> > >;

template <typename V>
using OffsetGeneratorsPtr = OffsetGenerators<V>*;

namespace {

template <typename V>
class OffsetGeneratorsFactory;

}  /* Anonymous namespace. */

template <typename Simplifier, typename V>
class LinearSet {
  public:
    LinearSet() {
      off_gens_ = factory_.NewOffsetGenerators(new OffsetGenerators<V>{
          SparseVec<V>{}, std::set< SparseVec<V> >{}});
    }

    LinearSet(const LinearSet &lset) = default;
    LinearSet(LinearSet &&lset) = default;

    LinearSet(const SparseVec<V> &o, const std::set< SparseVec<V> > &vs) {
      off_gens_ = factory_.NewOffsetGenerators(new OffsetGenerators<V>{o, vs});
    }

    LinearSet(SparseVec<V> &&o, std::set< SparseVec<V> > &&vs) {
      off_gens_ = factory_.NewOffsetGenerators(
          new OffsetGenerators<V>{std::move(o), std::move(vs)});
    }

    LinearSet(const SparseVec<V> &v) {
      off_gens_ = factory_.NewOffsetGenerators(new OffsetGenerators<V>{v, {}});
    }
    LinearSet(SparseVec<V> &&v) {
      off_gens_ = factory_.NewOffsetGenerators(
          new OffsetGenerators<V>{std::move(v), {}});
    }

    // FIXME: do we need this?
    // LinearSet& operator=(const LinearSet &s) = default;
    // LinearSet& operator=(LinearSet &&s) = default;

    ~LinearSet() = default;

    bool operator==(const LinearSet &rhs) const {
      return off_gens_ == rhs.off_gens_;
    }

    bool operator<(const LinearSet &rhs) const {
      return off_gens_ < rhs.off_gens_;
    }


    LinearSet operator+(const LinearSet &rhs) const {
      auto result = std::unique_ptr< OffsetGenerators<V> >{
        new OffsetGenerators<V>{GetOffset() + rhs.GetOffset(),
                                std::set< SparseVec<V> >{}}};

      std::set_union(GetGenerators().begin(), GetGenerators().end(),
                     rhs.GetGenerators().begin(), rhs.GetGenerators().end(),
                     inserter(result->second, result->second.begin()));

      SimplifySet(simplifier_, result->second);

      return LinearSet{factory_.NewOffsetGenerators(result.release())};
    }

    std::size_t Hash() const {
      std::hash< OffsetGeneratorsPtr<V> > h;
      return h(off_gens_);
    }

    friend std::ostream& operator<<(std::ostream &out, const LinearSet lset) {
      out << "<";
      out << lset.GetOffset();
      out << ":";
      out << "{";
      for (const auto &v : lset.GetGenerators()) {
        out << v;
      }
      out << "}";
      out << ">";
      return out;
    }

    const SparseVec<V>& GetOffset() const { return off_gens_->first; }
    const std::set< SparseVec<V> >& GetGenerators() const {
      return off_gens_->second;
    }

  private:
    LinearSet(OffsetGeneratorsPtr<V> ogs) : off_gens_(ogs) {}

    OffsetGeneratorsPtr<V> off_gens_;

    /* TODO: Try to get rid of static... */
    static OffsetGeneratorsFactory<V> factory_;
    static Simplifier simplifier_;

    template <typename S2, typename S1, typename VV>
    friend LinearSet<S2, VV> ChangeLinearSimplifier(const LinearSet<S1, VV> &lset);

};

template <typename Simplifier, typename V>
OffsetGeneratorsFactory<V> LinearSet<Simplifier, V>::factory_;

template <typename Simplifier, typename V>
Simplifier LinearSet<Simplifier, V>::simplifier_;

template <typename S2, typename S1, typename V>
LinearSet<S2, V> ChangeLinearSimplifier(const LinearSet<S1, V> &lset) {
  return LinearSet<S2, V>{lset.off_gens_};
}

namespace {

/*
 * FIXME: Currently we don't handle removing unused VarVectors.  This should be
 * pretty easy:
 * - make OffsetGenerators a class with ref counter and pointer to the factory,
 * - use intrusive pointer in the LinearSet, which when counter gets to 0 calls
 *   the factory to delete the mapping and then deletes the object.
 */
template <typename V>
class OffsetGeneratorsFactory {
  public:
    OffsetGeneratorsFactory() = default;

    ~OffsetGeneratorsFactory() {
      // std::cout << "Number of OffsetGenerator objects: " << map_.size() << std::endl;
      for (auto &key_value : map_) { delete key_value.second; }
    }

    OffsetGeneratorsFactory(const OffsetGeneratorsFactory &f) = delete;
    OffsetGeneratorsFactory(OffsetGeneratorsFactory &&f) = delete;

    OffsetGeneratorsFactory& operator=(const OffsetGeneratorsFactory &f) = delete;
    OffsetGeneratorsFactory& operator=(OffsetGeneratorsFactory &&f) = delete;

    OffsetGeneratorsPtr<V> NewOffsetGenerators(const OffsetGeneratorsPtr<V> off_gens) {
      assert(off_gens);
      auto iter_inserted =
        map_.emplace(KeyWrapper< OffsetGenerators<V> >{off_gens}, off_gens);
      /* Sanity check -- the actual vectors must be the same. */
      assert(*iter_inserted.first->second == *off_gens);
      return iter_inserted.first->second;
    }

  private:
    std::unordered_map< KeyWrapper< OffsetGenerators<V> >,
                        OffsetGeneratorsPtr<V> > map_;
};

}  /* Anonymous namespace. */



namespace std {

template<typename S, typename V>
struct hash< LinearSet<S, V> > {
  inline std::size_t operator()(const LinearSet<S, V> &set) const {
    return set.Hash();
  }

};

}  /* namespace std */
