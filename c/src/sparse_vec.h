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
template <typename V>
class SparseVec {
  typedef UniqueVMap<V, Counter> UniqueVMap_;
  typedef UniqueVMapPtr<V, Counter> UniqueVMapPtr_;
  typedef UniqueVMapBuilder<V, Counter> UniqueVMapBuilder_;

  public:
    SparseVec() {
      vmap_ = builder_.New({});
    }

    SparseVec(const SparseVec &v) = default;
    SparseVec(SparseVec &&v) = default;

    SparseVec& operator=(const SparseVec &v) = default;
    SparseVec& operator=(SparseVec &&v) = default;

    SparseVec(std::vector< std::pair<V, Counter> > &&vector) {
      vmap_ = builder_.New(std::move(vector));
    }

    SparseVec(std::initializer_list< std::pair<V, Counter> > list)
        : SparseVec(std::vector< std::pair<V, Counter> >{list}) {}

    SparseVec(const V &v, Counter c)
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

      unsigned int k = 0;

      for (auto &pair : *rhs.vmap_) {

        auto iter = vmap_->Find(pair.first);

        if (iter == vmap_->end()) {
          /* Domains are different! => this does not divide rhs. */
          return false;
        }

        // have we already obtained a candidate for a multiple?
        if (0 != k && pair.second != k * iter->second) {
            return false;
        }
        // no candidate for k yet
        else if (pair.second % iter->second != 0) {
          return false;
        }
        // division leaves no remainder -> candidate k found
        else {
          k = pair.second / iter->second;
        }
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

template <typename V>
UniqueVMapBuilder<V, Counter> SparseVec<V>::builder_;

namespace std {

template<typename V>
struct hash< SparseVec<V> > {
  inline std::size_t operator()(const SparseVec<V> &vec) const {
    return vec.Hash();
  }

};

}  /* namespace std */
