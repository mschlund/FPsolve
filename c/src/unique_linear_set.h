#pragma once

#include "hash.h"
#include "sparse_vec.h"

template <typename V>
class UniqueLinearSetBuilder;

template <typename V>
class UniqueLinearSet {
  public:
    UniqueLinearSet() = delete;
    UniqueLinearSet(const UniqueLinearSet &ulset) = delete;
    UniqueLinearSet(UniqueLinearSet &&ulset) = delete;

    UniqueLinearSet& operator=(const UniqueLinearSet &ulset) = delete;
    UniqueLinearSet& operator=(UniqueLinearSet &&ulset) = delete;

    ~UniqueLinearSet() = default;

    UniqueLinearSet(UniqueLinearSetBuilder<V> &b) : builder_(b) {}

    UniqueLinearSet(UniqueLinearSetBuilder<V> &b, SparseVec<V> &&o,
        std::set< SparseVec<V> > &&gs)
        : builder_(b), offset_(std::move(o)), generators_(std::move(gs)) {}

    UniqueLinearSet(UniqueLinearSetBuilder<V> &b, const SparseVec<V> &o,
        const std::set< SparseVec<V> > &gs)
        : builder_(b), offset_(o), generators_(gs) {}

    const SparseVec<V>& GetOffset() const { return offset_; }
    const std::set< SparseVec<V> >& GetGenerators() const { return generators_; }

    bool operator==(const UniqueLinearSet &rhs) const {
      return GetOffset() == rhs.GetOffset() &&
             GetGenerators() == rhs.GetGenerators();
    }

    // bool operator<(const UniqueLinearSet &rhs) const {
    //   return GetOffset() == rhs.GetOffset() &&
    //          GetGenerators() == rhs.GetGenerators();
    // }

    std::size_t Hash() const {
      std::size_t hash = 0;
      HashCombine(hash, offset_);
      HashCombine(hash, generators_);
      return hash;
    }

    friend inline void intrusive_ptr_add_ref(const UniqueLinearSet *ulset) {
      ++ulset->ref_count_;
    }

    friend inline void intrusive_ptr_release(const UniqueLinearSet *ulset) {
      assert(ulset->ref_count_ > 0);
      --ulset->ref_count_;
      if (ulset->ref_count_ == 0) {
        ulset->builder_.Delete(ulset);
      }
    }

    friend std::ostream& operator<<(std::ostream &out, const UniqueLinearSet &lset) {
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

  private:
    UniqueLinearSetBuilder<V> &builder_;

    SparseVec<V> offset_;
    std::set< SparseVec<V> > generators_;
    mutable RefCounter ref_count_ = 0;

    friend UniqueLinearSetBuilder<V>;
};

namespace std {

template<typename V>
struct hash< UniqueLinearSet<V> > {
  inline std::size_t operator()(const UniqueLinearSet<V> &lset) const {
    return lset.Hash();
  }
};

}  /* namespace std */


template <typename V>
using UniqueLinearSetPtr = IntrPtr< const UniqueLinearSet<V> >;

namespace std {

template<typename V>
struct hash< UniqueLinearSetPtr<V> > {
  inline std::size_t operator()(const UniqueLinearSetPtr<V> &ptr) const {
    std::hash<const UniqueLinearSet<V>*> h;
    return h(ptr.get());
  }
};

}  /* namespace std */


template <typename V>
class UniqueLinearSetBuilder {
  public:
    UniqueLinearSetBuilder() = default;

    ~UniqueLinearSetBuilder() {
      /* At this point there shouldn't be any outstanding pointers left... */
      assert(map_.empty());
      for (auto &key_value : map_) { DMSG(*key_value.second); delete key_value.second; }
    }

    UniqueLinearSetBuilder(const UniqueLinearSetBuilder &f) = delete;
    UniqueLinearSetBuilder(UniqueLinearSetBuilder &&f) = delete;

    UniqueLinearSetBuilder& operator=(const UniqueLinearSetBuilder &f) = delete;
    UniqueLinearSetBuilder& operator=(UniqueLinearSetBuilder &&f) = delete;

    UniqueLinearSetPtr<V> TryLookup(std::unique_ptr< UniqueLinearSet<V> > ulset) {
      assert(ulset);
      assert(ulset->ref_count_ == 0);
      auto iter_inserted =
        map_.emplace(KeyWrapper< UniqueLinearSet<V> >{ulset.get()}, ulset.get());
      /* Sanity check -- the actual vectors must be the same. */
      assert(*iter_inserted.first->second == *ulset);

      if (iter_inserted.second) {
        ulset.release();
      }

      return UniqueLinearSetPtr<V>{iter_inserted.first->second};
    }

    UniqueLinearSetPtr<V> New(SparseVec<V> &&o, std::set< SparseVec<V> > &&gs) {
      return TryLookup(std::unique_ptr< UniqueLinearSet<V> >(
            new UniqueLinearSet<V>(*this, std::move(o), std::move(gs))));
    }

    UniqueLinearSetPtr<V> New(const SparseVec<V> &o,
                              const std::set< SparseVec<V> > &gs) {
      return TryLookup(std::unique_ptr< UniqueLinearSet<V> >(
            new UniqueLinearSet<V>(*this, o, gs)));
    }

    void Delete(const UniqueLinearSet<V> *ulset) {
      auto iter = map_.find(KeyWrapper< UniqueLinearSet<V> >{ulset});
      if (iter != map_.end()) {
        delete iter->second;
        map_.erase(iter);
      } else {
        assert(false);
      }
    }

  private:
    std::unordered_map< KeyWrapper< UniqueLinearSet<V> >,
                        UniqueLinearSet<V>* > map_;
};
