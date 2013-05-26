#pragma once

#include <cassert>
#include <set>


#include "../datastructs/sparse_vec.h"
#include "../datastructs/unique_set.h"
#include "../string_util.h"

#include "semilinear_util.h"

#define LIN_SIMPL_TEMPLATE_TYPE \
  template < \
    typename, typename, \
    template <typename, typename> class, \
    template <typename, typename, \
              template <typename, typename> class> class \
  > class


class VarId;

/* Default LinearSet includes the simplification with SparseVecSimplifier and
 * uses VarId and Counter for variable identifier and counter. */
template <typename Var = VarId,
          typename Value = Counter,
          DIVIDER_TEMPLATE_TYPE VecDivider = DummyDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl = SparseVecSimplifier2>
class LinearSet;

/* SimpleLinearSet does no simplification. */
typedef LinearSet<VarId, Counter, DummyDivider, DummyVecSimplifier2> SimpleLinearSet;

/* DivLinearSet additionally divides the SparseVec by its gcd.  NOTE: this
 * over-approximates and does not give a precise answer anymore.  */
typedef LinearSet<VarId, Counter, GcdDivider> DivLinearSet;


template <typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
class LinearSet {
  public:
    typedef SparseVec<Var, Value, DummyDivider> OffsetType;
    typedef SparseVec<Var, Value, VecDivider> GeneratorType;
    typedef VecSimpl<Var, Value, VecDivider> VecSimplType;

    LinearSet() : offset_(), generators_(builder_.New({})) {}

    LinearSet(const LinearSet &lset) = default;
    LinearSet(LinearSet &&lset) = default;

    LinearSet(const OffsetType &o)
        : offset_(o), generators_(builder_.New({})) {}

    LinearSet(OffsetType &&o)
        : offset_(std::move(o)), generators_(builder_.New({})) {}

    LinearSet(const OffsetType &o, const VecSet<GeneratorType> &vs)
        : offset_(o), generators_(builder_.New(vs)) {}

    LinearSet(OffsetType &&o, VecSet<GeneratorType> &&vs)
        : offset_(std::move(o)), generators_(builder_.New(std::move(vs))) {}

    LinearSet(const OffsetType &o, VecSet<GeneratorType> &&vs)
        : offset_(o), generators_(builder_.New(std::move(vs))) {}

    template <DIVIDER_TEMPLATE_TYPE OldVecDivider,
              VEC_SIMPL_TEMPLATE_TYPE OldVecSimpl>
    LinearSet(const LinearSet<Var, Value, OldVecDivider, OldVecSimpl> &s) {
      VecSet<GeneratorType> generators;
      for (const auto &g : s.GetGenerators()) {
        /* Since s.GetGenerators() is already a VecSet, then the elements are
         * already ordered and it's safe to use emplace_back. */
        generators.emplace_back(GeneratorType{g});
      }
      offset_ = s.GetOffset();
      generators_ = builder_.New(std::move(generators));
    }

    LinearSet& operator=(const LinearSet &s) = default;
    LinearSet& operator=(LinearSet &&s) = default;

    ~LinearSet() = default;

    bool operator==(const LinearSet &rhs) const {
      return offset_ == rhs.offset_ && generators_ == rhs.generators_;
    }

    bool operator<(const LinearSet &rhs) const {
      if (offset_ < rhs.offset_) {
        return true;
      } else if (offset_ == rhs.offset_) {
        return generators_ < rhs.generators_;
      }
      return false;
    }


    LinearSet operator+(const LinearSet &rhs) const {
      if (IsZero()) {
        return rhs;
      } else if (rhs.IsZero()) {
        return *this;
      }
      auto result_offset = offset_ + rhs.offset_;
      VecSet<GeneratorType> result_generators;

      // VecSetUnion(GetGenerators(), rhs.GetGenerators());
      std::set_union(GetGenerators().begin(), GetGenerators().end(),
                     rhs.GetGenerators().begin(), rhs.GetGenerators().end(),
                     std::back_inserter(result_generators));
                     // inserter(result_generators, result_generators.begin()));

      SimplifySet<VecSimplType>(result_generators);

      return LinearSet{std::move(result_offset),
                       builder_.New(std::move(result_generators))};
    }

    std::size_t Hash() const {
      std::size_t hash = 0;
      HashCombine(hash, offset_);
      HashCombine(hash, generators_);
      return hash;
    }

    friend std::ostream& operator<<(std::ostream &out, const LinearSet lset) {
      out << "<" << lset.offset_
          << " : "
          << ToStringSorted(lset.GetGenerators()) << ">";
      return out;
    }

    const OffsetType& GetOffset() const { return offset_; }
    const VecSet<GeneratorType>& GetGenerators() const {
      assert(Ok());
      return generators_->GetVecSet();
    }

    bool IsZero() const { return offset_.IsZero() && generators_->empty(); }

  private:
    LinearSet(OffsetType &&o, UniqueSetPtr<GeneratorType> s)
        : offset_(std::move(o)), generators_(s) {}

    bool Ok() const {
      return generators_ != nullptr;
    }

    OffsetType offset_;
    UniqueSetPtr<GeneratorType> generators_;

    /* TODO: Try to get rid of static... */
    static UniqueSetBuilder<GeneratorType> builder_;
};

namespace std {

template<typename S, typename V>
struct hash< LinearSet<S, V> > {
  inline std::size_t operator()(const LinearSet<S, V> &set) const {
    return set.Hash();
  }

};

}  /* namespace std */

template <typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
UniqueSetBuilder< SparseVec<Var, Value, VecDivider> >
  LinearSet<Var, Value, VecDivider, VecSimpl>::builder_;

template <typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
class LinearSetSimplifier {
  public:
    typedef LinearSet<Var, Value, VecDivider, VecSimpl> LinearSetType;

    LinearSetSimplifier(const std::set<LinearSetType> &s) : rhs_lsets_(s) {}

    static bool IsActive() { return true; }

    bool IsCovered(const LinearSetType &lset) {
      for (auto &rhs_lset : rhs_lsets_) {
        VecSimpl<Var, Value, VecDivider> vec_simpl{rhs_lset.GetGenerators()};
        auto new_offset = lset.GetOffset() - rhs_lset.GetOffset();
        bool covered = new_offset.IsValid() &&
          vec_simpl.IsCovered(new_offset);
        if (!covered) {
          continue;
        }
        for (auto &lhs_gen : lset.GetGenerators()) {
          assert(covered);
          covered = vec_simpl.IsCovered(lhs_gen);
          if (!covered) {
            break;
          }
        }
        if (covered) {
          return true;
        }
      }
      return false;
    }
  private:
    const std::set<LinearSetType> &rhs_lsets_;
};

template <typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
class DummyLinSimplifier {
  public:
    typedef LinearSet<Var, Value, VecDivider, VecSimpl> LinearSetType;

    DummyLinSimplifier(const std::set<LinearSetType> &s) {}

    static bool IsActive() { return false; }

    bool IsCovered(const LinearSetType &lset) { return false; }
};


/*
 * This is the older idea of simplification where we only tested subset
 * inclusion of generators.
 */
template <typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
class LinearSubsetSimplifier {
  public:
    typedef LinearSet<Var, Value, VecDivider, VecSimpl> LinearSetType;

    LinearSubsetSimplifier(const std::set<LinearSetType> &s) : rhs_lsets_(s) {}

    static bool IsActive() { return true; }

    bool IsCovered(const LinearSetType &lset) {
      for (auto &rhs_lset : rhs_lsets_) {
        auto new_offset = lset.GetOffset() - rhs_lset.GetOffset();
        VecSimpl<Var, Value, VecDivider> vec_simpl{rhs_lset.GetGenerators()};
        if (new_offset.IsValid() &&
            vec_simpl.IsCovered(new_offset) &&
            std::includes(rhs_lset.GetGenerators().begin(),
                          rhs_lset.GetGenerators().end(),
                          lset.GetGenerators().begin(),
                          lset.GetGenerators().end())) {
          return true;
        }
      }
      return false;
    }
  private:
    const std::set<LinearSetType> &rhs_lsets_;
};
