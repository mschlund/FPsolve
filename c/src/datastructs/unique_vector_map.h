#pragma once

#include <initializer_list>

#include <boost/intrusive_ptr.hpp>

#include "key_wrapper.h"
#include "../debug_output.h"


#define DIVIDER_TEMPLATE_TYPE template <typename, typename> class

typedef std::uint_fast32_t RefCounter;

template <typename K, typename V>
class UniqueVMapBuilder;

template <typename K, typename V>
class UniqueVMapDivider;

template <typename K, typename V>
class GcdDivider;

template <typename K, typename V>
class UniqueVMap {
  typedef std::vector< std::pair<K, V> > Vector_;
  public:
    UniqueVMap() = delete;
    UniqueVMap(const UniqueVMap &vmap) = delete;
    UniqueVMap(UniqueVMap &&vmap) = delete;

    UniqueVMap& operator=(const UniqueVMap &vmap) = delete;
    UniqueVMap& operator=(UniqueVMap &&vmap) = delete;


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
    typename Vector_::size_type size() const { return vector_.size(); }

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
      assert(vec_map->ref_count_ > 0);
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

    void Clear() {
      vector_.clear();
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
    friend GcdDivider<K, V>;

    UniqueVMapBuilder<K, V> &builder_;
    mutable RefCounter ref_count_ = 0;

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

template <typename K, typename V>
class DummyDivider {
  public:
    inline void operator()(const UniqueVMap<K, V> &vec) const {}
};


template <typename K, typename V>
class GcdDivider {
  public:
    void operator()(UniqueVMap<K, V> &vec) const {
      if (vec.empty()) {
        return;
      }
      V divisor = vec.vector_.front().second;
      if (divisor == 1) {
        return;
      }
      DMSG("Gcd: Before: " << vec);
      // FIXME: When dealing with long vectors, maybe we could do something like
      // logarithmic number of Gcd calls -- first for every pair of numbers,
      // then for every pair of previous results, etc.
      // Make sure that this is sound...
      for (auto iter = vec.begin() + 1; iter != vec.end(); ++iter) {
        divisor = Gcd(divisor, iter->second);
        if (divisor == 1) {
          DMSG("Gcd: Nothing happens...");
          return;
        }
      }
      assert(divisor > 1);
      for (auto &x : vec.vector_) {
        assert(x.second % divisor == 0);
        x.second = x.second / divisor;
      }
      DMSG("Gcd: After: " << vec);
    }

  private:
    V Gcd(V a, V b) const {
      if (a < b) {
        std::swap(a, b);
      }
      V tmp;
      while (b != 0) {
        tmp = b;
        b = a % b;
        a = tmp;
      }
      return a;
    }
};



template <typename A>
using IntrPtr = boost::intrusive_ptr<A>;

template <typename K, typename V>
using UniqueVMapPtr = IntrPtr< const UniqueVMap<K, V> >;

template <typename K, typename V>
class UniqueVMapBuilder {

  struct Deleter {
    void operator()(UniqueVMap<K, V> *vmap) {
      vmap->builder_.Deallocate(vmap);
    }
  };

  public:
    UniqueVMapBuilder() = default;

    ~UniqueVMapBuilder() {
      for (auto &ptr : allocated_) { delete ptr; }
      /* At this point there shouldn't be any outstanding pointers left... */
      assert(map_.empty());
      for (auto &key_value : map_) { DMSG(*key_value.second); delete key_value.second; }
    }

    UniqueVMapBuilder(const UniqueVMapBuilder &f) = delete;
    UniqueVMapBuilder(UniqueVMapBuilder &&f) = delete;

    UniqueVMapBuilder& operator=(const UniqueVMapBuilder &f) = delete;
    UniqueVMapBuilder& operator=(UniqueVMapBuilder &&f) = delete;

    template <DIVIDER_TEMPLATE_TYPE Divider = DummyDivider>
    UniqueVMapPtr<K, V> TryLookup(std::unique_ptr< UniqueVMap<K, V>, Deleter > &&vmap) {
      assert(vmap);
      assert(vmap->ref_count_ == 0);

      Divider<K, V> divider;
      divider(*vmap);

      /* Note that the unique_ptr vmap still owns ptr. */
      auto ptr = vmap.get();

      auto iter_inserted =
        map_.emplace(KeyWrapper< UniqueVMap<K, V> >{ptr}, ptr);
      /* Sanity check -- the actual vectors must be the same. */
      assert(*iter_inserted.first->second == *vmap);

      /* If we actually inserted the ptr, vmap should release its ownership.
       * Otherwise the ~unique_ptr will call Deallocate. */
      if (iter_inserted.second) {
        vmap.release();
      }

      return UniqueVMapPtr<K, V>{iter_inserted.first->second};
    }

    template <DIVIDER_TEMPLATE_TYPE Divider = DummyDivider>
    UniqueVMapPtr<K, V> New(const UniqueVMap<K, V> &vmap) {
      auto result = Allocate();
      /* Copy over the vector, since it might be modified in-place by
       * the Divider. */
      result->vector_ = vmap.vector_;
      return TryLookup<Divider>(std::move(result));
    }

    template <DIVIDER_TEMPLATE_TYPE Divider = DummyDivider>
    UniqueVMapPtr<K, V> New(std::vector< std::pair<K, V> > &&input_vector) {
      std::sort(input_vector.begin(), input_vector.end());
      auto result = Allocate();

      for (auto &pair : input_vector) {
        if (result->vector_.empty() || result->vector_.back().first < pair.first) {
          result->vector_.emplace_back(std::move(pair));
        } else {
          assert(result->vector_.back().first == pair.first);
          result->vector_.back().second += std::move(pair.second);
        }
      }

      return TryLookup<Divider>(std::move(result));
    }

    template <DIVIDER_TEMPLATE_TYPE Divider = DummyDivider>
    UniqueVMapPtr<K, V> NewSum(const UniqueVMap<K, V> &lhs,
                               const UniqueVMap<K, V> &rhs) {

      auto result = Allocate();

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
          assert(lhs_iter->first == rhs_iter->first);
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

      return TryLookup<Divider>(std::move(result));
    }

    template <DIVIDER_TEMPLATE_TYPE Divider = DummyDivider>
    UniqueVMapPtr<K, V> NewDiff(const UniqueVMap<K, V> &lhs,
                                const UniqueVMap<K, V> &rhs) {

      if (lhs.vector_.size() < rhs.vector_.size()) {
        return UniqueVMapPtr<K, V>{nullptr};
      }

      auto result = Allocate();

      auto lhs_iter = lhs.vector_.begin();
      auto rhs_iter = rhs.vector_.begin();
      const auto lhs_iter_end = lhs.vector_.end();
      const auto rhs_iter_end = rhs.vector_.end();

      /* We must go through all the elements of rhs (but not neccesarily of
       * lhs, i.e., lhs might be "larger"). */
      while (rhs_iter != rhs_iter_end) {

        /* If there's nothing left in lhs, or there is unmatched variable in
         * rhs, return null pointer. */
        if (lhs_iter == lhs_iter_end || lhs_iter->first > rhs_iter->first) {
          return UniqueVMapPtr<K, V>{nullptr};
        }

        if (lhs_iter->first < rhs_iter->first) {
          result->vector_.emplace_back(*lhs_iter);
          ++lhs_iter;
          continue;
        }

        assert(lhs_iter->first == rhs_iter->first);
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
        /*
         * If lhs value is equal to rhs value we just advance to the next variable
         */

        ++lhs_iter;
        ++rhs_iter;
      }

      /* If we didn't go through the whole lhs, copy whatever remained. */
      std::copy(lhs_iter, lhs_iter_end, std::back_inserter(result->vector_));

      return TryLookup<Divider>(std::move(result));
    }

    void Delete(const UniqueVMap<K, V> *vector_ptr) {
      auto iter = map_.find(KeyWrapper< UniqueVMap<K, V> >{vector_ptr});
      if (iter != map_.end()) {
        Deallocate(iter->second);
        map_.erase(iter);
      } else {
        assert(false);
      }
    }

  private:
    std::unique_ptr< UniqueVMap<K, V>, Deleter > Allocate() {
      if (allocated_.empty()) {
        return std::unique_ptr< UniqueVMap<K, V>, Deleter >(
                 new UniqueVMap<K, V>(*this));
      }
      std::unique_ptr< UniqueVMap<K, V>, Deleter > ptr{allocated_.back()};
      allocated_.pop_back();
      return ptr;
    }

    void Deallocate(UniqueVMap<K, V> *ptr) {
      assert(ptr->ref_count_ == 0);
      /* Clear the vector so there's no garbage left.  This should _not_
       * deallocate the buffer. */
      ptr->Clear();

      if (allocated_.size() < 1024) {
        allocated_.push_back(ptr);
      } else {
        delete ptr;
      }
    }


    std::unordered_map< KeyWrapper< UniqueVMap<K, V> >, UniqueVMap<K, V>* > map_;
    std::vector<UniqueVMap<K, V>*> allocated_;
};

namespace std {

template<typename K, typename V>
struct hash< UniqueVMapPtr<K, V> > {
  inline std::size_t operator()(const UniqueVMapPtr<K, V> &ptr) const {
    std::hash<const UniqueVMap<K, V>*> h;
    return h(ptr.get());
  }
};




}  /* namespace std */
