#pragma once

#include "hash.h"
#include "vec_set.h"


template <typename A>
class UniqueSetBuilder;

template <typename A>
class UniqueSet {
  public:
    typedef typename VecSet<A>::iterator iterator;
    typedef typename VecSet<A>::const_iterator const_iterator;

    UniqueSet() = delete;
    UniqueSet(const UniqueSet &set) = delete;
    UniqueSet(UniqueSet &&set) = delete;

    UniqueSet& operator=(const UniqueSet &set) = delete;
    UniqueSet& operator=(UniqueSet &&set) = delete;

    ~UniqueSet() = default;

    UniqueSet(UniqueSetBuilder<A> &b) : builder_(b) {}

    UniqueSet(UniqueSetBuilder<A> &b, VecSet<A> &&s)
        : builder_(b), set_(std::move(s)) {}

    UniqueSet(UniqueSetBuilder<A> &b, const VecSet<A> &s)
        : builder_(b), set_(s) {}

    const VecSet<A>& GetVecSet() const { return set_; }

    bool empty() const { return set_.empty(); }

    bool operator==(const UniqueSet &rhs) const {
      return set_ == rhs.set_;
    }

    bool operator<(const UniqueSet &rhs) const {
      return set_ < rhs.set_;
    }

    std::size_t Hash() const {
      std::hash< VecSet<A> > h;
      return h(set_);
    }

    iterator begin() { return set_.begin(); }
    iterator end() { return set_.end(); }

    const_iterator begin() const { return set_.begin(); }
    const_iterator end() const { return set_.end(); }

    friend inline void intrusive_ptr_add_ref(const UniqueSet *set) {
      ++set->ref_count_;
    }

    friend inline void intrusive_ptr_release(const UniqueSet *set) {
      assert(set->ref_count_ > 0);
      --set->ref_count_;
      if (set->ref_count_ == 0) {
        set->builder_.Delete(set);
      }
    }

    friend std::ostream& operator<<(std::ostream &out, const UniqueSet &set) {
      out << "{";
      for (const auto &v : set) {
        out << v;
      }
      out << "}";
      return out;
    }

  private:
    UniqueSetBuilder<A> &builder_;
    VecSet<A> set_;
    mutable RefCounter ref_count_ = 0;

    friend UniqueSetBuilder<A>;
};


namespace std {

template<typename A>
struct hash< UniqueSet<A> > {
  inline std::size_t operator()(const UniqueSet<A> &set) const {
    return set.Hash();
  }
};

}  /* namespace std */


template <typename A>
using UniqueSetPtr = IntrPtr< const UniqueSet<A> >;

namespace std {

template<typename A>
struct hash< UniqueSetPtr<A> > {
  inline std::size_t operator()(const UniqueSetPtr<A> &ptr) const {
    std::hash<const UniqueSet<A>*> h;
    return h(ptr.get());
  }
};

}  /* namespace std */


template <typename A>
class UniqueSetBuilder {
  public:
    UniqueSetBuilder() = default;

    ~UniqueSetBuilder() {
      /* At this point there shouldn't be any outstanding pointers left... */
      assert(map_.empty());
      for (auto &key_value : map_) { delete key_value.second; }
    }

    UniqueSetBuilder(const UniqueSetBuilder &b) = delete;
    UniqueSetBuilder(UniqueSetBuilder &&b) = delete;

    UniqueSetBuilder& operator=(const UniqueSetBuilder &b) = delete;
    UniqueSetBuilder& operator=(UniqueSetBuilder &&b) = delete;

    UniqueSetPtr<A> TryLookup(std::unique_ptr< UniqueSet<A> > &&set) {
      assert(set);
      assert(set->ref_count_ == 0);
      auto iter_inserted =
        map_.emplace(KeyWrapper< UniqueSet<A> >{set.get()}, set.get());
      /* Sanity check -- the actual vectors must be the same. */
      assert(*iter_inserted.first->second == *set);

      if (iter_inserted.second) {
        set.release();
      }

      return UniqueSetPtr<A>{iter_inserted.first->second};
    }

    UniqueSetPtr<A> New(VecSet<A> &&s) {
      return TryLookup(std::unique_ptr< UniqueSet<A> >(
            new UniqueSet<A>(*this, std::move(s))));
    }

    UniqueSetPtr<A> New(const VecSet<A> &s) {
      return TryLookup(std::unique_ptr< UniqueSet<A> >(
            new UniqueSet<A>(*this, s)));
    }

    void Delete(const UniqueSet<A> *set) {
      auto iter = map_.find(KeyWrapper< UniqueSet<A> >{set});
      if (iter != map_.end()) {
        delete iter->second;
        map_.erase(iter);
      } else {
        assert(false);
      }
    }

  private:
    std::unordered_map< KeyWrapper< UniqueSet<A> >, UniqueSet<A>* > map_;
};
