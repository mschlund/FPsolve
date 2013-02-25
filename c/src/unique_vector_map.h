#pragma once

#include <initializer_list>

#include "key_wrapper.h"

typedef std::uint_fast32_t RefCounter;

template <typename K, typename V>
class UniqueVectorMapBuilder;

template <typename K, typename V>
class UniqueVectorMap {
  typedef std::vector< std::pair<K, V> > Vector_;
  public:
    UniqueVectorMap() = delete;

    UniqueVectorMap(UniqueVectorMapBuilder<K, V> &b, Vector_ &&input_vector)
        : builder_(b) {

      std::sort(input_vector.begin(), input_vector.end());

      for (auto &pair : input_vector) {
        if (vector_.empty() || vector_.back().first < pair.first) {
          vector_.emplace_back(std::move(pair));
        } else {
          assert(vector_.back().first == pair.first);
          vector_.back().second += std::move(pair.second);
        }
      }
    }

    UniqueVectorMap(UniqueVectorMapBuilder<K, V> &b, const V &v, Counter c)
        : builder_(b), vector_({ std::make_pair(v, c) }) {}

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

    typename Vector_::const_iterator begin() const { return vector_.begin(); }

    typename Vector_::const_iterator end() const { return vector_.end(); }

    bool operator==(const UniqueVectorMap &rhs) const {
      return vector_ == rhs.vector_;
    }

    bool operator<(const UniqueVectorMap &rhs) const {
      return vector_ < rhs.vector_;
    }


    friend inline void intrusive_ptr_add_ref(const UniqueVectorMap *vec_map) {
      ++vec_map->ref_count_;
    }

    friend inline void intrusive_ptr_release(const UniqueVectorMap *vec_map) {
      --vec_map->ref_count_;
      if (vec_map->ref_count_ == 0) {
        vec_map->builder_.Delete(vec_map);
      }
    }

    friend std::ostream& operator<<(std::ostream &out,
                                    const UniqueVectorMap &vec_map) {
      out << "[";
      for (const auto &key_value : vec_map) {
        out << "(" << key_value.first << ", " << key_value.second << ")";
      }
      out << "]";
      return out;
    }

  private:
    struct Less {
      bool operator()(const std::pair<K, V> &key_value, const K &key)
        const { return key_value.first < key; }
    };

    UniqueVectorMapBuilder<K, V> &builder_;
    mutable RefCounter ref_count_;

    // FIXME: should that be const?
    Vector_ vector_;
};

namespace std {

template<typename K, typename V>
struct hash< UniqueVectorMap<K, V> > {
  inline std::size_t operator()(const UniqueVectorMap<K, V> &vec_map) const {
    return vec_map.Hash();
  }

};

}  /* namespace std */


template <typename A>
using IntrPtr = boost::intrusive_ptr<A>;

template <typename K, typename V>
using UniqueVectorMapPtr = IntrPtr< UniqueVectorMap<K, V> >;

template <typename K, typename V>
class UniqueVectorMapBuilder {
  public:
    UniqueVectorMapBuilder() = default;

    ~UniqueVectorMapBuilder() {
      std::cout << "Number of UniqueVectorMap objects: " << map_.size() << std::endl;
      for (auto &key_value : map_) { delete key_value.second; }
    }

    UniqueVectorMapBuilder(const UniqueVectorMapBuilder &f) = delete;
    UniqueVectorMapBuilder(UniqueVectorMapBuilder &&f) = delete;

    UniqueVectorMapBuilder& operator=(const UniqueVectorMapBuilder &f) = delete;
    UniqueVectorMapBuilder& operator=(UniqueVectorMapBuilder &&f) = delete;

    UniqueVectorMapPtr<K, V> New(std::vector< std::pair<K, V> > &&vector) {
      auto vector_ptr = new UniqueVectorMap<K, V>(*this, std::move(vector));
      assert(vector_ptr);
      auto iter_inserted =
        map_.emplace(KeyWrapper< UniqueVectorMap<K, V> >{vector_ptr}, vector_ptr);
      /* Sanity check -- the actual vectors must be the same. */
      assert(*iter_inserted.first->second == *vector_ptr);
      return UniqueVectorMapPtr<K, V>{iter_inserted.first->second};
    }

    UniqueVectorMapPtr<K, V> New(std::initializer_list< std::pair<K, V> > list) {
      return New(std::vector< std::pair<K, V> >{list});
    }

    UniqueVectorMapPtr<K, V> New(const K &key, const V &value) {
      return New({ {key, value} });
    }

    UniqueVectorMapPtr<K, V> NewMerge(const UniqueVectorMap<K, V> &lhs,
                                      const UniqueVectorMap<K, V> &rhs) const {
    }

    UniqueVectorMapPtr<K, V> NewSubtract(const UniqueVectorMap<K, V> &lhs,
                                         const UniqueVectorMap<K, V> &rhs) const {
    }

    void Delete(const UniqueVectorMap<K, V> *vector_ptr) {
      std::cout << "UniqueVectorMap: deleting: " << *vector_ptr << std::endl;
      auto iter = map_.find(KeyWrapper< UniqueVectorMap<K, V> >{vector_ptr});
      assert(iter != map_.end());
      auto tmp_ptr = iter->second;
      map_.erase(iter);
      delete tmp_ptr;
    }

  private:
    std::unordered_map< KeyWrapper< UniqueVectorMap<K, V> >,
                        UniqueVectorMap<K, V>* > map_;
};

