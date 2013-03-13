#pragma once

#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "key_wrapper.h"
#include "semilinear_util.h"
#include "sparse_vec.h"
#include "unique_set.h"
#include "var.h"

/* Default LinearSet includes the simplification with SparseVecSimplifier and
 * uses VarPtr and Counter for variable identifier and counter. */
template <typename Var = VarPtr,
          typename Value = Counter,
          typename VecDivider = DummyDivider,
          typename VecSimpl = SparseVecSimplifier<Var, Value, VecDivider> >
class LinearSet;

/* SimpleLinearSet does no simplification. */
typedef LinearSet<VarPtr, Counter, DummyDivider, DummySimplifier> SimpleLinearSet;

/* DivLinearSet additionally divides the SparseVec by its gcd.  NOTE: this
 * over-approximates and does not give a precise answer anymore.  */
typedef LinearSet<VarPtr, Counter, GcdDivider<VarPtr, Counter> > DivLinearSet;


template <typename Var,
          typename Value,
          typename VecDivider,
          typename VecSimpl>
class LinearSet {
  public:
    typedef SparseVec<Var, Value, DummyDivider> OffsetType;
    typedef SparseVec<Var, Value, VecDivider> GeneratorType;

    LinearSet() : offset_(), generators_(builder_.New({})) {}

    LinearSet(const LinearSet &lset) = default;
    LinearSet(LinearSet &&lset) = default;

    LinearSet(const OffsetType &o)
        : offset_(o), generators_(builder_.New({})) {}

    LinearSet(OffsetType &&o)
        : offset_(std::move(o)), generators_(builder_.New({})) {}

    LinearSet(const OffsetType &o, const std::set<GeneratorType> &vs)
        : offset_(o), generators_(builder_.New(vs)) {}

    LinearSet(OffsetType &&o, std::set<GeneratorType> &&vs)
        : offset_(std::move(o)), generators_(builder_.New(std::move(vs))) {}

    LinearSet(const OffsetType &o, std::set<GeneratorType> &&vs)
        : offset_(o), generators_(builder_.New(std::move(vs))) {}

    template <typename OldVecDivider, typename OldVecSimpl>
    LinearSet(const LinearSet<Var, Value, OldVecDivider, OldVecSimpl> &s) {
      std::set<GeneratorType> generators;
      for (const auto &g : s.GetGenerators()) {
        generators.insert(GeneratorType{g});
      }
      LinearSet(s.GetOffset(), std::move(generators));
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
      auto result_offset = offset_ + rhs.offset_;
      std::set<GeneratorType> result_generators;

      std::set_union(GetGenerators().begin(), GetGenerators().end(),
                     rhs.GetGenerators().begin(), rhs.GetGenerators().end(),
                     inserter(result_generators, result_generators.begin()));

      SimplifySet<VecSimpl>(result_generators);

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
      out << "<" << lset.offset_ << " : " << *lset.generators_ << ">";
      return out;
    }

    const OffsetType& GetOffset() const { return offset_; }
    const std::set<GeneratorType>& GetGenerators() const {
      return generators_->GetSet();
    }

  private:
    LinearSet(OffsetType &&o, UniqueSetPtr<GeneratorType> s)
        : offset_(std::move(o)), generators_(s) {}

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

template <typename Var, typename Value, typename VecDivider, typename VecSimpl>
UniqueSetBuilder< SparseVec<Var, Value, VecDivider> >
  LinearSet<Var, Value, VecDivider, VecSimpl>::builder_;

template <typename Var,
          typename Value,
          typename VecDivider,
          typename VecSimpl>
class LinearSetSimplifier {
  public:
    typedef LinearSet<Var, Value, VecDivider, VecSimpl> LinearSetType;

    LinearSetSimplifier(const std::set<LinearSetType> &s) : rhs_lsets_(s) {}

    static bool IsActive() { return true; }

    bool IsCovered(const LinearSetType &lset) {
      for (auto &rhs_lset : rhs_lsets_) {
        VecSimpl vec_simpl{rhs_lset.GetGenerators()};
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

/*
 * This is the older idea of simplification where we only tested subset
 * inclusion of generators.
 */
template <typename Var, typename Value, typename VecDivider, typename VecSimpl>
class LinearSubsetSimplifier {
  public:
    typedef LinearSet<Var, Value, VecDivider, VecSimpl> LinearSetType;

    LinearSubsetSimplifier(const std::set<LinearSetType> &s) : rhs_lsets_(s) {}

    static bool IsActive() { return true; }

    bool IsCovered(const LinearSetType &lset) {
      for (auto &rhs_lset : rhs_lsets_) {
        auto new_offset = lset.GetOffset() - rhs_lset.GetOffset();
        VecSimpl vec_simpl{rhs_lset.GetGenerators()};
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
