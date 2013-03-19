#pragma once

#include <set>
#include <sstream>

#include "semilinear_util.h"
#include "sparse_vec.h"

#include "debug_output.h"
#include "equations.h"
#include "string_util.h"

class VarId;

template <typename Var = VarId,
          typename Value = Counter,
          DIVIDER_TEMPLATE_TYPE VecDivider = DummyDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl = SparseVecSimplifier>
class PseudoLinearSet;

typedef PseudoLinearSet<VarId, Counter, GcdDivider> DivPseudoLinearSet;

static const std::uint_fast16_t simpl_freq = 1;

/*
 * An abstraction of semilinear sets that is quite similar to a linear set, but
 * with multiple offsets.  This allows to over-approximate the semilinear set.
 */
template <typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
class PseudoLinearSet : public Semiring< PseudoLinearSet<Var, Value, VecDivider,
                                                         VecSimpl> > {
  public:
    typedef VecSimpl<Var, Value, VecDivider> VecSimplType;

    typedef SparseVec<Var, Value, DummyDivider> OffsetType;
    typedef SparseVec<Var, Value, VecDivider> GeneratorType;

    PseudoLinearSet() = default;

    PseudoLinearSet(const Var &v, const Counter &c)
        : offsets_({ OffsetType{v, c} }) {}

    PseudoLinearSet(std::initializer_list<OffsetType> l) : offsets_(l) {}

    PseudoLinearSet(const PseudoLinearSet &s) = default;
    PseudoLinearSet(PseudoLinearSet &&s) = default;

    /* Allow construction from a LinearSet with a different divider and
     * simplifier. */
    template <DIVIDER_TEMPLATE_TYPE OldVecDivider,
              VEC_SIMPL_TEMPLATE_TYPE OldVecSimpl>
    PseudoLinearSet(const LinearSet<Var, Value, OldVecDivider, OldVecSimpl> &lset) {
      offsets_.insert(OffsetType{lset.GetOffset()});
      for (const auto &g : lset.GetGenerators()) {
        generators_.insert(GeneratorType{g});
      }
      Simplify();
    }

    PseudoLinearSet& operator=(const PseudoLinearSet &s) = default;
    PseudoLinearSet& operator=(PseudoLinearSet &&s) = default;

    ~PseudoLinearSet() = default;

    static PseudoLinearSet null() { return PseudoLinearSet{}; }
    static PseudoLinearSet one() { return PseudoLinearSet{{}}; }

    bool operator==(const PseudoLinearSet &rhs) const {
      return offsets_ == rhs.offsets_ && generators_ == rhs.generators_;
    }

    PseudoLinearSet& operator+=(const PseudoLinearSet &rhs) {
      std::set<OffsetType> result_offsets;
      std::set<GeneratorType> result_generators;

      std::set_union(offsets_.begin(), offsets_.end(),
                     rhs.offsets_.begin(), rhs.offsets_.end(),
                     std::inserter(result_offsets, result_offsets.begin()));

      std::set_union(generators_.begin(), generators_.end(),
                     rhs.generators_.begin(), rhs.generators_.end(),
                     std::inserter(result_generators, result_generators.begin()));

      offsets_ = std::move(result_offsets);
      generators_ = std::move(result_generators);

      without_simpl = std::max(without_simpl, rhs.without_simpl) + 1;
      Simplify();

      return *this;
    }

    PseudoLinearSet& operator*=(const PseudoLinearSet &rhs) {
      std::set<OffsetType> result_offsets;
      std::set<GeneratorType> result_generators;

      for (auto &vec_rhs : rhs.offsets_) {
        for (auto &vec_lhs : offsets_) {
          result_offsets.insert(vec_lhs + vec_rhs);
        }
      }

      std::set_union(generators_.begin(), generators_.end(),
                     rhs.generators_.begin(), rhs.generators_.end(),
                     std::inserter(result_generators, result_generators.begin()));

      offsets_ = std::move(result_offsets);
      generators_ = std::move(result_generators);

      without_simpl = std::max(without_simpl, rhs.without_simpl) + 1;
      Simplify();

      return *this;
    }

    PseudoLinearSet star() const {

      std::set<OffsetType> result_offsets = { OffsetType{} };
      std::set<GeneratorType> result_generators = generators_;

      for (auto &offset : offsets_) {
        /* Convert from OffsetType to GeneratorType. */
        result_generators.insert(GeneratorType{offset});
      }

      return PseudoLinearSet{ std::move(result_offsets),
                              std::move(result_generators) };
    }


    std::string string() const {
      std::stringstream ss;
      ss << "{";
      ss << ToStringSorted(offsets_, " ");
      ss << " : ";
      ss << ToStringSorted(generators_, " ");
      ss << " }";
      return ss.str();
    }

    friend std::ostream& operator<<(std::ostream &out,
                                    const PseudoLinearSet &set) {
      out << set.string();
      return out;
    }

    static const bool is_idempotent = true;
    static const bool is_commutative = true;

  private:
    PseudoLinearSet(std::set<OffsetType> &&os, std::set<GeneratorType> &&gs)
        : offsets_(os), generators_(gs) {}

    void Simplify() {
      if (!VecSimplType::IsActive() || without_simpl < simpl_freq) {
        return;
      }
      without_simpl = 0;

#ifdef DEBUG_OUTPUT
      DMSG("Offsets: {");
      for (auto &x : offsets_) {
        DMSG(x);
      }
      DMSG("}");
      DMSG("Generators: {");
      for (auto &x : generators_) {
        DMSG(x);
      }
      DMSG("}");
#endif

      /* Simplify generators. */
      SimplifySet<VecSimplType>(generators_);

      /* Simplify offsets.  This uses the same trick as SimplifySet -- erase in
       * std::set does not invalidate any iterators (except for the removed). */
      VecSimpl<Var, Value, VecDivider> simplifier{generators_};
      for (auto offset1_iter = offsets_.begin(); offset1_iter != offsets_.end(); ) {
        auto offset1 = *offset1_iter;
        /* Erase automatically advances the iterator to the next element. */
        offset1_iter = offsets_.erase(offset1_iter);

        bool necessary = true;
        for (auto &offset2 : offsets_) {
          auto new_offset = offset1 - offset2;
          if (new_offset.IsValid() &&
              simplifier.IsCovered(new_offset)) {
            necessary = false;
            break;
          }
        }

        if (necessary) {
          offsets_.insert(std::move(offset1));
        } else {
          DMSG("Removed offset: " << offset1);
        }
      }
    }

    std::set<OffsetType> offsets_;
    std::set<GeneratorType> generators_;

    std::uint_fast16_t without_simpl = 0;

    template <DIVIDER_TEMPLATE_TYPE OldVecDivider,
              VEC_SIMPL_TEMPLATE_TYPE OldVecSimpl>
    PseudoLinearSet SemilinearToPseudoLinear(
        const SemilinearSet<Var, Value, OldVecDivider, OldVecSimpl> &semilinear);
};


template <DIVIDER_TEMPLATE_TYPE NewVecDivider,
          VEC_SIMPL_TEMPLATE_TYPE NewVecSimpl,
          typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
PseudoLinearSet<Var, Value, NewVecDivider, NewVecSimpl> SemilinearToPseudoLinear(
    const SemilinearSet<Var, Value, VecDivider, VecSimpl> &semilinear) {
  PseudoLinearSet<Var, Value, NewVecDivider, NewVecSimpl> pseudo_lset;
  for (const auto &lset : semilinear) {
    pseudo_lset += PseudoLinearSet<Var, Value, NewVecDivider, NewVecSimpl>{lset};
  };
  return pseudo_lset;
}

template <DIVIDER_TEMPLATE_TYPE NewVecDivider,
          VEC_SIMPL_TEMPLATE_TYPE NewVecSimpl,
          typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
Equations< PseudoLinearSet<Var, Value, NewVecDivider, NewVecSimpl> >
SemilinearToPseudoLinearEquations(
    const Equations< SemilinearSet<Var, Value, VecDivider, VecSimpl> > &equations) {
  return MapEquations(equations,
    [](const SemilinearSet<Var, Value, VecDivider, VecSimpl> &s) {
      return SemilinearToPseudoLinear<NewVecDivider, NewVecSimpl>(s);
    });
}
