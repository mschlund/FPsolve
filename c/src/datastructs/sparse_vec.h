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
#include "var.h"
#include "vec_set.h"
#include "../utils/string_util.h"


#define VEC_SIMPL_TEMPLATE_TYPE \
  template <typename, typename, template <typename, typename> class> class

typedef std::uint_fast32_t Counter;


/*
 * Since SparseVec has a template parameter Divider, then we can't simply have a
 * static UniqueVMapBuilder there, because if there are SparseVec's with
 * different Dividers they will have different UniqueVMapBuilders.  We want
 * UniqueVMapBuilders to be global (and be able to easily change the Dividers
 * without modifying the underlying UniqueVMap).  Having this proxy
 * GlobalVMapBuilder solves the problem.
 */
template <typename Var, typename Value>
class GlobalVMapBuilder {
  public:
    static UniqueVMapBuilder<Var, Value>& Get() { return builder_; }
  private:
    static UniqueVMapBuilder<Var, Value> builder_;
};

template <typename Var, typename Value>
UniqueVMapBuilder<Var, Value> GlobalVMapBuilder<Var, Value>::builder_;


/*
 * Sparse vector representing the mapping from variables to counters.  We never
 * create a vector that has already be created, therefore operations such as
 * equality are really just pointer equality (i.e., they are very efficient).
 * Moreover we never modify a vector, so copy constructor, assignment operator,
 * etc. are really only copying a pointer (which is again very efficient).
 */
template <typename Var = VarId,
          typename Value = Counter,
          DIVIDER_TEMPLATE_TYPE Divider = DummyDivider>
class SparseVec {
  typedef UniqueVMap<Var, Value> UniqueVMap_;
  typedef UniqueVMapPtr<Var, Value> UniqueVMapPtr_;
  typedef UniqueVMapBuilder<Var, Value> UniqueVMapBuilder_;
  typedef GlobalVMapBuilder<Var, Value> GlobalVMapBuilder_;

  public:
    SparseVec()
        : vmap_(GlobalVMapBuilder_::Get().template New<Divider>(
            std::vector< std::pair<Var, Value> >{})) {}

    SparseVec(const SparseVec &v) = default;
    SparseVec(SparseVec &&v) = default;

    SparseVec(std::vector< std::pair<Var, Value> > &&vector)
        : vmap_(GlobalVMapBuilder_::Get().template New<Divider>(std::move(vector))) {}

    SparseVec(std::initializer_list< std::pair<Var, Value> > list)
        : SparseVec(std::vector< std::pair<Var, Value> >{list}) {}

    SparseVec(const Var &v, Value c)
      : SparseVec({ std::make_pair(v, c) }) {}

    /* Make it possible to construct a SparseVec from one with a different
     * divider.  For that we need to use builder_ to actually invoke the
     * new Divider.  Make it explicit since we don't want to do those
     * convertions by accident. */
    template <DIVIDER_TEMPLATE_TYPE OldDivider>
    explicit SparseVec(const SparseVec<Var, Value, OldDivider> &v)
        : vmap_(GlobalVMapBuilder_::Get().template New<Divider>(*v.vmap_)) {}

    /* Allow casting a SparseVec to one with a different Divider, but without
     * actually invoking the new Divider.  This is useful for simplifiers that
     * just need to check set membership or subtract vectors... */
    template <DIVIDER_TEMPLATE_TYPE AnyDivider>
    SparseVec<Var, Value, AnyDivider> Cast() const {
      return SparseVec<Var, Value, AnyDivider>{vmap_};
    }

    SparseVec& operator=(const SparseVec &v) = default;
    SparseVec& operator=(SparseVec &&v) = default;

    bool operator==(const SparseVec &rhs) const {
      assert(Sanity() && rhs.Sanity());
      if (vmap_ == nullptr || rhs.vmap_ == nullptr) {
        return false;
      }
      assert((vmap_ == rhs.vmap_ && Equal(rhs)) ||
             (vmap_ != rhs.vmap_ && !Equal(rhs)));
      return vmap_ == rhs.vmap_;
    }

    bool operator!=(const SparseVec &rhs) const { return !(*this == rhs); }

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
      return GlobalVMapBuilder_::Get().template NewSum<Divider>(*vmap_, *rhs.vmap_);
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
      /* We're currently running NewDiff without simplifying the vector with a
       * divider.  This is probably both unnecessary and would probably make the
       * simplification (SparseVecSimplifier) less likely to optimize away
       * unnecessary vector. */
      return GlobalVMapBuilder_::Get().template NewDiff<DummyDivider>(*vmap_, *rhs.vmap_);
    }

    bool IsZero() const {
      assert(Sanity());
      return IsValid() ? vmap_->empty() : false;
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

      Counter k = 0;

      for (auto lhs_iter = vmap_->begin(), rhs_iter = rhs.vmap_->begin();
           lhs_iter != vmap_->end();
           ++lhs_iter, ++rhs_iter) {

        assert(rhs_iter != rhs.vmap_->end());

        if (/* If the domains are different. */
            lhs_iter->first != rhs_iter->first ||
            /* Or we have a multiple and it doesn't work here. */
            (k != 0 && rhs_iter->second != k * lhs_iter->second) ||
            /* Or we can't divide. */
            rhs_iter->second % lhs_iter->second != 0) {
          return false;
        }

        /* Otherwise we have a candidate for the multiple. */
        k = rhs_iter->second / lhs_iter->second;

      }
      return true;
    }


    friend std::ostream& operator<<(std::ostream &out, const SparseVec &svector) {
      assert(svector.Sanity());
      out << "[";
      out << ToStringSorted(*svector.vmap_);
      out << "]";
      return out;
    }

    void RawString(std::ostream &out) const {
      out << *vmap_;
    }

    std::size_t Hash() const {
      std::hash<UniqueVMapPtr_> h;
      return h(vmap_);
    }

    bool Equal(const SparseVec &rhs) const {
      return *vmap_ == *rhs.vmap_;
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

    /* Make other SparseVec a friend so that we can convert between them, e.g.,
     * when changing the divider. */
    template <typename SomeVar, typename SomeValue,
              DIVIDER_TEMPLATE_TYPE SomeDivider>
    friend class SparseVec;
};

namespace std {

template<typename Var, typename Value, DIVIDER_TEMPLATE_TYPE Divider>
struct hash< SparseVec<Var, Value, Divider> > {
  inline std::size_t operator()(const SparseVec<Var, Value, Divider> &vec) const {
    return vec.Hash();
  }

};

}  /* namespace std */

template<typename Var, typename Value, DIVIDER_TEMPLATE_TYPE Divider>
class DummyVecSimplifier {
  public:
    typedef SparseVec<Var, Value, Divider> SparseVecType;
    DummyVecSimplifier(const VecSet<SparseVecType> &s) {}

    static bool IsActive() { return false; }

    bool IsCovered(const SparseVecType &lhs_any) {
      return false;
    }
};

template <typename Var, typename Value, DIVIDER_TEMPLATE_TYPE Divider>
class NaiveSimplifier {
  public:
    typedef SparseVec<Var, Value, Divider> SparseVecType;

    NaiveSimplifier(const VecSet<SparseVecType> &s) : rhs_set_(s) {}

    static bool IsActive() { return true; }

    template <DIVIDER_TEMPLATE_TYPE AnyDivider>
    bool IsCovered(const SparseVec<Var, Value, AnyDivider> &lhs_any) {
      /* Cast the vector without actually invoking the Divider. */
      auto lhs = lhs_any.template Cast<Divider>();
      if (lhs.IsZero() || rhs_set_.Contains(lhs)) {
        return true;
      }
      for (const auto &rhs_gen : rhs_set_) {
        if (rhs_gen.Divides(lhs)) {
          return true;
        }
      }
      return false;
    }

  private:
    const VecSet<SparseVecType> &rhs_set_;
};

template <typename Var, typename Value, DIVIDER_TEMPLATE_TYPE Divider>
class SparseVecSimplifier {
  public:
    typedef SparseVec<Var, Value, Divider> SparseVecType;
    typedef SparseVec<Var, Value, DummyDivider> SparseVecTypeNoDiv;

    SparseVecSimplifier(const VecSet<SparseVecType> &s)
        : rhs_set_(s) {}

    static bool IsActive() { return true; }

    /*
     * Check, if lhs is a non-negative integer linear combination of the vectors
     * in rhs_set_
     */
    template <template <typename, typename> class AnyDivider>
    bool IsCovered(const SparseVec<Var, Value, AnyDivider> &lhs_any) {
      auto lhs = lhs_any.template Cast<DummyDivider>();
      /* Check the cheap and naive simplifier. */
      NaiveSimplifier<Var, Value, Divider> naive{rhs_set_};
      if (naive.IsCovered(lhs)) {
        return true;
      }
      return IsCovered_(lhs);
    }

  private:
    /* Dynamic programming/memoization */
    bool IsCovered_(const SparseVecTypeNoDiv &lhs) {
      auto iter = computed_.find(lhs);
      if (iter != computed_.end()) {
        return iter->second;
      }
      /* Should we go from the back? */
      for (auto &rhs : rhs_set_) {
        if (rhs.IsZero()) {
          assert(false);
          continue;
        }
        /* We don't want to do any simplification/division of the vector. */
        auto new_lhs = lhs - rhs.template Cast<DummyDivider>();
        if (new_lhs.IsValid() && (new_lhs.IsZero() || IsCovered_(new_lhs))) {
          computed_.emplace(lhs, true);
          return true;
        }
      }
      computed_.emplace(lhs, false);
      return false;
    }

    const VecSet<SparseVecType> &rhs_set_;
    std::unordered_map<SparseVecTypeNoDiv, bool> computed_;
};
