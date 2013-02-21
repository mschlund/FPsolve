#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <vector>

#include "hash.h"
#include "key_wrapper.h"

typedef std::uint_fast32_t Counter;

template <typename V>
class SparseVec;

namespace {

template <typename V>
using VarVector = std::vector< std::pair<V, Counter> >;

template <typename V>
using VarVectorPtr = VarVector<V>*;

template <typename V>
class VarVectorFactory;

}  /* Anonymous namespace. */

/*
 * Sparse vector representing the mapping from variables to counters.  We never
 * create a vector that has already be created, therefore operations such as
 * equality are really just pointer equality (i.e., they are very efficient).
 * Moreover we never modify a vector, so copy constructor, assignment operator,
 * etc. are really only copying a pointer (which is again very efficient).
 */
template <typename V>
class SparseVec {
  public:
    SparseVec() {
      vector_ptr_ = factory_.NewVarVector(new VarVector<V>);
    }

    SparseVec(const SparseVec &v) = default;
    SparseVec(SparseVec &&v) = default;

    SparseVec& operator=(const SparseVec &v) = default;
    SparseVec& operator=(SparseVec &&v) = default;

    SparseVec(VarVector<V> &&vector) {
      std::sort(vector.begin(), vector.end());
      std::unique_ptr< VarVector<V> > result{new VarVector<V>};

      for (auto &pair : vector) {
        if (result->empty()) {
          result->emplace_back(pair);
        } else if (result->back().first < pair.first) {
          result->emplace_back(pair);
        } else {
          assert(result->back().first == pair.first);
          result->back().second += pair.second;
        }
      }

      vector_ptr_ = factory_.NewVarVector(result.release());
    }

    SparseVec(std::initializer_list< std::pair<V, Counter> > list)
        : SparseVec(VarVector<V>{list}) {}

    SparseVec(const V &v, Counter c)
      : SparseVec({ std::pair<V, Counter>{v, c} }) {}

    bool operator==(const SparseVec &rhs) const {
      return vector_ptr_ == rhs.vector_ptr_;
    }

    bool operator!=(const SparseVec &rhs) const {
      return vector_ptr_ != rhs.vector_ptr_;
    }

    bool operator<(const SparseVec &rhs) const {
      return vector_ptr_ < rhs.vector_ptr_;
    }

    SparseVec operator+(const SparseVec &rhs) const {
      std::unique_ptr< VarVector<V> > result{new VarVector<V>};

      auto lhs_iter = vector_ptr_->begin();
      auto rhs_iter = rhs.vector_ptr_->begin();
      const auto lhs_iter_end = vector_ptr_->end();
      const auto rhs_iter_end = rhs.vector_ptr_->end();

      while (lhs_iter != lhs_iter_end && rhs_iter != rhs_iter_end) {
        if (lhs_iter->first < rhs_iter->first) {
          result->emplace_back(*lhs_iter);
          ++lhs_iter;
        } else if (lhs_iter->first > rhs_iter->first) {
          result->emplace_back(*rhs_iter);
          ++rhs_iter;
        } else {
          /* lhs_iter->first == rhs_iter->first */
          result->emplace_back(lhs_iter->first, lhs_iter->second + rhs_iter->second);
          ++lhs_iter;
          ++rhs_iter;
        }
      }

      for (; lhs_iter != lhs_iter_end; ++lhs_iter) {
        result->emplace_back(*lhs_iter);
      }

      for (; rhs_iter != rhs_iter_end; ++rhs_iter) {
        result->emplace_back(*rhs_iter);
      }

      return SparseVec{factory_.NewVarVector(result.release())};
    }

    friend std::ostream& operator<<(std::ostream &out, const SparseVec &svector) {
      out << "[";
      for (const auto &pair : *svector.vector_ptr_) {
        out << "(" << pair.first << ", " << pair.second << ")";
      }
      out << "]";
      return out;
    }

    std::size_t Hash() const {
      std::hash< VarVectorPtr<V> > h;
      return h(vector_ptr_);
    }

  private:
    SparseVec(VarVectorPtr<V> v) : vector_ptr_(v) {}

    VarVectorPtr<V> vector_ptr_;
    static VarVectorFactory<V> factory_;
};

template <typename V>
VarVectorFactory<V> SparseVec<V>::factory_;

namespace {

template <typename V>
class VarVectorFactory {
  public:
    VarVectorFactory() = default;

    ~VarVectorFactory() {
      std::cout << "Number of SparseVec objects: " << map_.size() << std::endl;
      for (auto &key_value : map_) { delete key_value.second; }
    }

    VarVectorFactory(const VarVectorFactory &f) = delete;
    VarVectorFactory(VarVectorFactory &&f) = delete;

    VarVectorFactory& operator=(const VarVectorFactory &f) = delete;
    VarVectorFactory& operator=(VarVectorFactory &&f) = delete;

    VarVectorPtr<V> NewVarVector(const VarVectorPtr<V> vector_ptr) {
      assert(vector_ptr);
      auto iter_inserted =
        map_.emplace(KeyWrapper< VarVector<V> >{vector_ptr}, vector_ptr);
      /* Sanity check -- the actual vectors must be the same. */
      assert(*iter_inserted.first->second == *vector_ptr);
      return iter_inserted.first->second;
    }

  private:
    std::unordered_map< KeyWrapper< VarVector<V> >, VarVectorPtr<V> > map_;
};

}  /* Anonymous namespace. */

namespace std {

template<typename V>
struct hash< SparseVec<V> > {
  inline std::size_t operator()(const SparseVec<V> &vec) const {
    return vec.Hash();
  }

};

}  /* namespace std */


