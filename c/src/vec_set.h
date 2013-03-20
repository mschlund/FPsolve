#pragma once

#include <algorithm>
#include <cassert>

#include <boost/iterator/iterator_facade.hpp>

#include "debug_output.h"
#include "hash.h"

/*
 * The idea behind VecSet is to have a structure that supports efficient inorder
 * traversal and union.
 */

template <typename A>
class VecSet;

namespace {

typedef std::size_t Index;

static const Index nil_index = std::numeric_limits<Index>::max();

template <typename A>
struct ElemInfo {
  ElemInfo() = default;
  ElemInfo(const ElemInfo &e) = default;
  ElemInfo(ElemInfo &&e) = default;

  ElemInfo(const A &v, Index p, Index n) : value(v), prev(p), next(n) {}
  ElemInfo(A &&v, Index p, Index n) : value(std::move(v)), prev(p), next(n) {}

  ElemInfo& operator=(const ElemInfo &e) = default;
  ElemInfo& operator=(ElemInfo &&e) = default;


  A value;
  Index prev;
  Index next;
  bool erased = false;
};


template <bool B, typename A>
struct AddConst;

template <typename A>
struct AddConst<true, A> {
  typedef const A type;
};

template <typename A>
struct AddConst<false, A> {
  typedef A type;
};

}  /* Anonymous namespace */

template <bool Const, typename A>
class VecIter : public boost::iterator_facade<VecIter<Const, A>,
                                              typename AddConst<Const, A>::type,
                                              boost::forward_traversal_tag> {
  typedef typename AddConst< Const, std::vector< ElemInfo<A> > >::type
          Vector_;
  typedef typename AddConst< Const, A >::type A_;
  public:
    VecIter() : vector_(nullptr), index_(nil_index) {}

    explicit VecIter(Vector_ *v, Index i)
        : vector_(v), index_(i) {}

  private:

    inline void increment() {
      assert(Ok());
      index_ = (*vector_)[index_].next;
    }

    inline bool equal(const VecIter &rhs) const {
      assert(vector_ == rhs.vector_);
      return index_ == rhs.index_;
    }

    inline A_& dereference() const {
      assert(vector_ != nullptr && index_ < vector_->size());
      return (*vector_)[index_].value;
    }

    friend class boost::iterator_core_access;
    friend class VecSet<A>;

    bool Ok() const {
      if (index_ >= vector_->size() || (*vector_)[index_].erased) {
        return false;
      }
      return true;
    }



    Vector_ *vector_;
    Index index_;
};

template <typename A>
class VecSet {
  public:
    typedef VecIter<false, A> iterator;
    typedef VecIter<true, A> const_iterator;

    typedef A value_type;

    VecSet() {}

    VecSet(const VecSet &vs) = default;
    VecSet(VecSet &&vs) = default;

    VecSet& operator=(const VecSet &vs) = default;
    VecSet& operator=(VecSet &&vs) = default;

    VecSet(std::vector<A> &&v) {
      std::sort(v.begin(), v.end());
      Index next = 0;
      for (Index i = 0; i < v.size(); ++i) {
        /* Skip duplicates. */
        if (next > 0 && vector_[next - 1].value == v[i]) {
          continue;
        }
        assert(vector_.size() == next);
        emplace_back(std::move(v[i]));
        ++next;
      }
      assert(Ok());
    }

    VecSet(std::initializer_list<A> v) : VecSet(std::vector<A>{v}) {}

    bool Contains(const A &elem) const {
      /* The idea here is to use the fact that the elements are in sorted order
       * (even if some of them are marked "erased").  So we can do standard
       * binary search and then check if it hasn't been erased. */
      struct Less {
        bool operator()(const ElemInfo<A> &lhs, const A &rhs) {
          return lhs.value < rhs;
        }
      };
      Less less;
      auto iter = std::lower_bound(vector_.begin(), vector_.end(), elem, less);
      if (iter != vector_.end() && !iter->erased && iter->value == elem) {
        return true;
      }
      return false;
    }

    iterator MarkErased(iterator iter) {
      assert(marked_erased == nil_index);
      ElemInfo<A> &ei = vector_[iter.index_];
      ei.erased = true;
      if (iter.index_ == begin_) {
        begin_ = ei.next;
      }
      ++num_erased_;
      if (ei.next != vector_.size()) {
        vector_[ei.next].prev = ei.prev;
      }
      if (ei.prev != nil_index) {
        vector_[ei.prev].next = ei.next;
      }
      marked_erased = iter.index_;
      assert(Ok());
      return iterator{&vector_, ei.next};
    }

    void CommitErase() {
      marked_erased = nil_index;
      assert(Ok());
    }

    void AbortErase() {
      assert(marked_erased != nil_index);
      ElemInfo<A> &ei = vector_[marked_erased];
      ei.erased = false;
      if (marked_erased < begin_) {
        begin_ = marked_erased;
      }
      --num_erased_;
      if (ei.next != vector_.size()) {
        vector_[ei.next].prev = marked_erased;
      }
      if (ei.prev != nil_index) {
        vector_[ei.prev].next = marked_erased;
      }
      marked_erased = nil_index;
      assert(Ok());
    }

    inline std::size_t Hash() const {
      std::size_t h = 0;
      for (auto &x : *this) {
        HashCombine(h, x);
      }
      return h;
    }

    Index size() const {
      return vector_.size() - num_erased_;
    }

    bool empty() const {
      return size() == 0;
    }

    bool operator<(const VecSet &rhs) const {
      return std::lexicographical_compare(this->begin(), this->end(),
                                          rhs.begin(), rhs.end());
    }

    bool operator==(const VecSet &rhs) const {
      if (size() != rhs.size()) {
        return false;
      }
      return std::equal(this->begin(), this->end(), rhs.begin());
    }

    iterator begin() {
      return iterator{&vector_, begin_};
    }

    iterator end() {
      return iterator{&vector_, vector_.size()};
    }

    const_iterator begin() const {
      return const_iterator{&vector_, begin_};
    }

    const_iterator end() const {
      return const_iterator{&vector_, vector_.size()};
    }

    void push_back(const A &a) {
      assert(marked_erased == nil_index);
      emplace_back(a);
    }

    void emplace_back(const A &a) {
      assert(marked_erased == nil_index);
      emplace_back(A{a});
    }

    // FIXME: should we track the last element of the list or just require that
    // emplace_back cannot be called when there are erased elements in it..?
    void emplace_back(A &&a) {
      assert(marked_erased == nil_index);
      Index index = vector_.size();
      Index prev_index = index == 0 ? nil_index : index - 1;
      vector_.emplace_back(ElemInfo<A>{std::move(a), prev_index, index + 1});
      assert(Ok());
    }

  private:
    void Debug() const {
      DMSG("empty() == " << empty());
      DMSG("begin_ == " << begin_);
      DMSG("begin().index_ == " << begin().index_);
      DMSG("end().index_ == " << end().index_);
      DMSG("num_erased_ == " << num_erased_);
      DMSG("vector_.size() == " << vector_.size());
      for (auto &ei : vector_) {
        DMSG(ei.prev << " " << ei.next << " " << ei.value);
      }
    }

    bool Ok() const {
      if (empty() && begin_ == vector_.size() && begin() == end()) {
        return true;
      } else if (empty() || begin_ == vector_.size() || begin() == end()) {
        DMSG("VecSet: Ok: empty");
        Debug();
        return false;
      }

      Index prev = begin_;
      for (Index current = begin_;
           current < vector_.size();
           prev = current, current = vector_[current].next) {
        if (current >= vector_.size()) {
          DMSG("VecSet: Ok: currupted index");
          Debug();
          return false;
        }
        if (current == begin_ && vector_[current].prev != nil_index) {
          DMSG("VecSet: Ok: begin.prev is wrong");
          Debug();
          return false;
        }
        if (current > begin_ && vector_[prev].next != current) {
          DMSG("VecSet: Ok: prev.next is wrong");
          Debug();
          return false;
        }
        if (current > begin_ && vector_[current].prev != prev) {
          DMSG("VecSet: Ok: current.prev is wrong");
          Debug();
          return false;
        }
        if (vector_[current].erased) {
          DMSG("VecSet: Ok: erased node reachable");
          Debug();
          return false;
        }
      }
      if (vector_[prev].next != end().index_ ||
          vector_[prev].next != vector_.size()) {
        DMSG("VecSet: Ok: wrong end index");
        Debug();
        return false;
      }
      Index count = 0;
      for (auto &x : vector_) {
        if (x.erased) {
          ++count;
        }
      }
      if (count != num_erased_) {
        DMSG("Wrong num_erased");
          Debug();
        return false;
      }
      return true;
    }

    friend class VecIter<false, A>;
    friend class VecIter<true, A>;


    std::vector< ElemInfo<A> > vector_;
    Index begin_ = 0;
    Index marked_erased = nil_index;
    Index num_erased_ = 0;
};


/* Helper function for creating a new VecSet that is a union of the arguments.
 * The only difference with std::set_union is that we explicitly use A to
 * insert an element of type B. */
template<typename A, typename B>
VecSet<A> VecSetUnion(const VecSet<A> &lhs, const VecSet<B> &rhs) {
  DMSG("VecSetUnion");
  auto lhs_iter = lhs.begin();
  auto lhs_iter_end = lhs.end();
  auto rhs_iter = rhs.begin();
  auto rhs_iter_end = rhs.end();
  VecSet<A> result;
  auto output = std::back_inserter(result);
  while (lhs_iter != lhs_iter_end) {
    if (rhs_iter == rhs_iter_end) {
      std::copy(lhs_iter, lhs_iter_end, output);
      return result;
    }
    auto rhs_elem = A{*rhs_iter};
    if (rhs_elem < *lhs_iter) {
      *output = rhs_elem;
      ++rhs_iter;
    } else {
      *output = *lhs_iter;
      if (!(*lhs_iter < rhs_elem)) {
        ++rhs_iter;
      }
      ++lhs_iter;
    }
  }
  for (; rhs_iter != rhs_iter_end; ++rhs_iter) {
    *output = A{*rhs_iter};
  }
  return result;
}


namespace std {

template<typename A>
struct hash< VecSet<A> > {
  inline std::size_t operator()(const VecSet<A> &set) const {
    std::size_t h = 0;
    for (auto &x : set) {
      HashCombine(h, x);
    }
    return h;
  }
};

}  /* namespace std */

