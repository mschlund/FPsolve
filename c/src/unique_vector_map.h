#pragma once

#include <initializer_list>

#include <boost/intrusive_ptr.hpp>

#include "key_wrapper.h"

typedef std::uint_fast32_t RefCounter;

template <typename K, typename V>
class UniqueVMapBuilder;

template <typename K, typename V>
class UniqueVMap {
  typedef std::vector< std::pair<K, V> > Vector_;
  public:
    UniqueVMap() = delete;

    typename Vector_::const_iterator Find(const K &key) const {
      auto iter = std::lower_bound(vector_.begin(), vector_.end(), key, Less{});
      if (iter != vector_.end() && iter->first == key) {
        return iter;
      }
      return vector_.end();
    }

    std::size_t Hash() const {
      std::hash<Vector_> h;
      return h(vector_);
    }

    bool empty() const { return vector_.empty(); }
    bool size() const { return vector_.size(); }

    typename Vector_::const_iterator begin() const { return vector_.begin(); }

    typename Vector_::const_iterator end() const { return vector_.end(); }

    bool operator==(const UniqueVMap &rhs) const {
      return vector_ == rhs.vector_;
    }

    bool operator<(const UniqueVMap &rhs) const {
      return vector_ < rhs.vector_;
    }


    friend inline void intrusive_ptr_add_ref(const UniqueVMap *vec_map) {
      ++vec_map->ref_count_;
    }

    friend inline void intrusive_ptr_release(const UniqueVMap *vec_map) {
      --vec_map->ref_count_;
      if (vec_map->ref_count_ == 0) {
        vec_map->builder_.Delete(vec_map);
      }
    }

    friend std::ostream& operator<<(std::ostream &out,
                                    const UniqueVMap &vec_map) {
      out << "[";
      for (const auto &key_value : vec_map) {
        out << "(" << key_value.first << ", " << key_value.second << ")";
      }
      out << "]";
      return out;
    }

  private:
    UniqueVMap(UniqueVMapBuilder<K, V> &b) : builder_(b) {}
    UniqueVMap(UniqueVMapBuilder<K, V> &b, Vector_ &&v)
        : builder_(b), vector_(std::move(v)) {
      assert(Sanity());
    }


    bool Sanity() const {
      for (auto &var_count : vector_) {
        if (var_count.second == 0) {
          return false;
        }
      }
      return true;
    }

    struct Less {
      bool operator()(const std::pair<K, V> &key_value, const K &key)
        const { return key_value.first < key; }
    };

    friend UniqueVMapBuilder<K, V>;

    UniqueVMapBuilder<K, V> &builder_;
    mutable RefCounter ref_count_ = 0;

    // FIXME: should that be const?
    Vector_ vector_;
};

namespace std {

template<typename K, typename V>
struct hash< UniqueVMap<K, V> > {
  inline std::size_t operator()(const UniqueVMap<K, V> &vec_map) const {
    return vec_map.Hash();
  }

};

}  /* namespace std */


template <typename A>
using IntrPtr = boost::intrusive_ptr<A>;

template <typename K, typename V>
using UniqueVMapPtr = IntrPtr< UniqueVMap<K, V> >;

template <typename K, typename V>
class UniqueVMapBuilder {
  public:
    UniqueVMapBuilder() = default;

    ~UniqueVMapBuilder() {
      std::cout << "Number of UniqueVMap objects: " << map_.size() << std::endl;
      for (auto &key_value : map_) { delete key_value.second; }
    }

    UniqueVMapBuilder(const UniqueVMapBuilder &f) = delete;
    UniqueVMapBuilder(UniqueVMapBuilder &&f) = delete;

    UniqueVMapBuilder& operator=(const UniqueVMapBuilder &f) = delete;
    UniqueVMapBuilder& operator=(UniqueVMapBuilder &&f) = delete;

    UniqueVMapPtr<K, V> TryLookup(UniqueVMap<K, V> *vmap) {
      assert(vmap);
      auto iter_inserted =
        map_.emplace(KeyWrapper< UniqueVMap<K, V> >{vmap}, vmap);
      /* Sanity check -- the actual vectors must be the same. */
      assert(*iter_inserted.first->second == *vmap);
      return UniqueVMapPtr<K, V>{iter_inserted.first->second};
    }

    UniqueVMapPtr<K, V> New(std::vector< std::pair<K, V> > &&input_vector) {
      std::sort(input_vector.begin(), input_vector.end());
      auto result =
        std::unique_ptr< UniqueVMap<K, V> >{new UniqueVMap<K, V>(*this)};

      for (auto &pair : input_vector) {
        if (result->vector_.empty() || result->vector_.back().first < pair.first) {
          result->vector_.emplace_back(std::move(pair));
        } else {
          assert(result->vector_.back().first == pair.first);
          result->vector_.back().second += std::move(pair.second);
        }
      }

      return TryLookup(result.release());
    }

    UniqueVMapPtr<K, V> NewSum(const UniqueVMap<K, V> &lhs,
                               const UniqueVMap<K, V> &rhs) {

      std::unique_ptr< UniqueVMap<K, V> > result{new UniqueVMap<K, V>(*this)};

      auto lhs_iter = lhs.vector_.begin();
      auto rhs_iter = rhs.vector_.begin();
      const auto lhs_iter_end = lhs.vector_.end();
      const auto rhs_iter_end = rhs.vector_.end();

      while (lhs_iter != lhs_iter_end && rhs_iter != rhs_iter_end) {
        if (lhs_iter->first < rhs_iter->first) {
          result->vector_.emplace_back(*lhs_iter);
          ++lhs_iter;
        } else if (lhs_iter->first > rhs_iter->first) {
          result->vector_.emplace_back(*rhs_iter);
          ++rhs_iter;
        } else {
          /* lhs_iter->first == rhs_iter->first */
          result->vector_.emplace_back(lhs_iter->first,
                                       lhs_iter->second + rhs_iter->second);
          ++lhs_iter;
          ++rhs_iter;
        }
      }

      for (; lhs_iter != lhs_iter_end; ++lhs_iter) {
        result->vector_.emplace_back(*lhs_iter);
      }

      for (; rhs_iter != rhs_iter_end; ++rhs_iter) {
        result->vector_.emplace_back(*rhs_iter);
      }

      return TryLookup(result.release());
    }

    UniqueVMapPtr<K, V> NewDiff(const UniqueVMap<K, V> &lhs,
                                const UniqueVMap<K, V> &rhs) {

      if (lhs.vector_.size() < rhs.vector_.size()) {
        return UniqueVMapPtr<K, V>{nullptr};
      }

      std::unique_ptr< UniqueVMap<K, V> > result{new UniqueVMap<K, V>(*this)};

      auto lhs_iter = lhs.vector_.begin();
      auto rhs_iter = rhs.vector_.begin();
      const auto lhs_iter_end = lhs.vector_.end();
      const auto rhs_iter_end = rhs.vector_.end();

      /* We must go through all the elements of rhs (but not neccesarily
       * lhs). */
      while (rhs_iter != rhs_iter_end) {

        /* If there's nothing left in lhs, or there is unmatched variable in
         * rhs, return invalid SparseVec. */
        if (lhs_iter == lhs_iter_end || lhs_iter->first > rhs_iter->first) {
          return UniqueVMapPtr<K, V>{nullptr};
        }

        if (lhs_iter->first < rhs_iter->first) {
          result->vector_.emplace_back(*lhs_iter);
          ++lhs_iter;
          continue;
        }

        /* lhs_iter->first == rhs_iter->first */
        if (lhs_iter->second < rhs_iter->second) {
          return UniqueVMapPtr<K, V>{nullptr};
        }

        /* If the lhs value is greater, subtract the value from rhs, if they're
         * equal, do nothing (i.e., the result is zero so we don't add anything
         * to the vector. */
        if (lhs_iter->second > rhs_iter->second) {
          result->vector_.emplace_back(lhs_iter->first,
                                       lhs_iter->second - rhs_iter->second);
        }

        ++lhs_iter;
        ++rhs_iter;
      }

      /* If we didn't go through the whole lhs, add what remained. */
      for (; lhs_iter != lhs_iter_end; ++lhs_iter) {
        result->vector_.emplace_back(*lhs_iter);
      }

      return TryLookup(result.release());
    }

    void Delete(const UniqueVMap<K, V> *vector_ptr) {
      auto iter = map_.find(KeyWrapper< UniqueVMap<K, V> >{vector_ptr});
      if (iter != map_.end()) {
        auto tmp_ptr = iter->second;
        assert(tmp_ptr == vector_ptr);
        map_.erase(iter);
        delete tmp_ptr;
      } else {
        assert(false);
      }
    }

  private:
    std::unordered_map< KeyWrapper< UniqueVMap<K, V> >, UniqueVMap<K, V>* > map_;
};

namespace std {

template<typename K, typename V>
struct hash< UniqueVMapPtr<K, V> > {
  inline std::size_t operator()(const UniqueVMapPtr<K, V> &ptr) const {
    std::hash<UniqueVMap<K, V>*> h;
    return h(ptr.get());
  }
};

}  /* namespace std */

