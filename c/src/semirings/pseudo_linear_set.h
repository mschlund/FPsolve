#pragma once

#include <set>
#include <sstream>

#include "../datastructs/sparse_vec.h"
#include "../datastructs/equations.h"

#include "../utils/debug_output.h"
#include "../utils/string_util.h"
#include "../utils/profiling-macros.h"

#include "semilinear_util.h"
#include "linear_set.h"
#include "semilinear_set.h"

class VarId;

template <typename Var = VarId,
          typename Value = Counter,
          DIVIDER_TEMPLATE_TYPE VecDivider = DummyDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl = SparseVecSimplifier>
class PseudoLinearSet;

typedef PseudoLinearSet<VarId, Counter, GcdDivider> DivPseudoLinearSet;

namespace {

static const std::uint_fast16_t simpl_freq = 0;


}  /* Anonymous namespace */


template <typename A>
using SetType = VecSet<A>;


/*
 * An abstraction of semilinear sets that is quite similar to a linear set, but
 * with multiple offsets.  This allows to over-approximate the semilinear set.
 */
template <typename Var,
          typename Value,
          DIVIDER_TEMPLATE_TYPE VecDivider,
          VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
class PseudoLinearSet : public StarableSemiring< PseudoLinearSet<Var, Value, VecDivider,
                                                         VecSimpl>,
                                                         Commutativity::Commutative, Idempotence::Idempotent > {
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
    PseudoLinearSet(const Var &v) : PseudoLinearSet(v, 1) {}


    /* Allow construction from a LinearSet with a different divider and
     * simplifier. */
    template <DIVIDER_TEMPLATE_TYPE OldVecDivider,
              VEC_SIMPL_TEMPLATE_TYPE OldVecSimpl>
    PseudoLinearSet(const LinearSet<Var, Value, OldVecDivider, OldVecSimpl> &lset) {
      offsets_.emplace_back(OffsetType{lset.GetOffset()});
      std::vector<GeneratorType> gens;
      for (const auto &g : lset.GetGenerators()) {
        gens.emplace_back(GeneratorType{g});
      }
      generators_ = VecSet<GeneratorType>{std::move(gens)};
      Simplify();
    }


    PseudoLinearSet(const SemilinearSet<Var, Value, VecDivider, VecSimpl>& slset) {
      std::vector<GeneratorType> gens;
      for(const auto &ls : slset) {
        offsets_.emplace_back(OffsetType{ls.GetOffset()});
        for (const auto &g : ls.GetGenerators()) {
          gens.emplace_back(GeneratorType{g});
        }
      }
      generators_ = VecSet<GeneratorType>{std::move(gens)};
      Simplify();
    }

    // parse a description of the form e.g. "<a:3,b:2>" or of the form "a"
    PseudoLinearSet(const std::string &str_val) : PseudoLinearSet(SemilinearSet<Var, Value, VecDivider, VecSimpl>(str_val)) { };


    PseudoLinearSet& operator=(const PseudoLinearSet &s) = default;
    PseudoLinearSet& operator=(PseudoLinearSet &&s) = default;

    ~PseudoLinearSet() = default;

    static PseudoLinearSet null() { return PseudoLinearSet{}; }
    bool IsZero() const { return offsets_.empty() && generators_.empty(); }
    static PseudoLinearSet one() { return PseudoLinearSet{{}}; }
    bool IsOne() const {
      return generators_.empty() &&
             offsets_.size() == 1 &&
             offsets_.begin()->IsZero();
    }



#ifdef USE_GENEPI
    bool operator==(const PseudoLinearSet &rhs) const {
      if (this->IsZero() && rhs.IsZero() || this->IsOne() && rhs.IsOne())
        return true;

      if (this->IsZero() && rhs.IsOne() || this->IsOne() && rhs.IsZero())
        return false;

      auto V1 = this->getVariables();
      auto V2 = rhs.getVariables();
      if (V1 != V2)
        return false;

      int k = V1.size();

      //std::cout << "numvars: " << k << std::endl;

      SemilinSetNdd::solver_init(k);
      //std::cout << *this << " ?== " << rhs << std::endl;
      //std::cout << "testing for eq of pseudolin in slsetNDD" << std::endl;
      bool eq = (SemilinSetNdd(this->offsets_, this->generators_) == SemilinSetNdd(rhs.offsets_, rhs.generators_));
      //std::cout << "eq solved: " << *this << " ?==" << rhs << ": " << eq << std::endl;
      SemilinSetNdd::solver_dealloc();
      return eq;
    }
#else

    bool operator==(const PseudoLinearSet &rhs) const {
      return offsets_ == rhs.offsets_ && generators_ == rhs.generators_;
    }
#endif

    PseudoLinearSet& operator+=(const PseudoLinearSet &rhs) {
      OPADD;
      if (rhs.IsZero()) {
        return *this;
      } else if (IsZero()) {
        *this = rhs;
        return *this;
      }

      SetType<OffsetType> result_offsets;
      SetType<GeneratorType> result_generators;

      std::set_union(offsets_.begin(), offsets_.end(),
                     rhs.offsets_.begin(), rhs.offsets_.end(),
                     std::back_inserter(result_offsets));

      std::set_union(generators_.begin(), generators_.end(),
                     rhs.generators_.begin(), rhs.generators_.end(),
                     std::back_inserter(result_generators));

      offsets_ = std::move(result_offsets);
      generators_ = std::move(result_generators);

      without_simpl = std::max(without_simpl, rhs.without_simpl) + 1;
      Simplify();
      return *this;
    }

    PseudoLinearSet& operator*=(const PseudoLinearSet &rhs) {
      OPMULT;
      if (IsZero() || rhs.IsZero()) {
        *this = null();
        return *this;
      }

      if (rhs.IsOne()) {
        return *this;
      } else if (IsOne()) {
        *this = rhs;
        return *this;
      }

      std::vector<OffsetType> result_offsets_vec;
      SetType<GeneratorType> result_generators;

      for (auto &vec_rhs : rhs.offsets_) {
        for (auto &vec_lhs : offsets_) {
          result_offsets_vec.emplace_back(vec_lhs + vec_rhs);
        }
      }

      std::set_union(generators_.begin(), generators_.end(),
                     rhs.generators_.begin(), rhs.generators_.end(),
                     std::back_inserter(result_generators));

      offsets_ = VecSet<OffsetType>{std::move(result_offsets_vec)};
      generators_ = std::move(result_generators);

      without_simpl = std::max(without_simpl, rhs.without_simpl) + 1;
      Simplify();
      return *this;
    }

    PseudoLinearSet star() const {
      DMSG("-> star");
      SetType<OffsetType> result_offsets = { OffsetType{} };

      auto result_generators = VecSetUnionWith(generators_, offsets_,
        [](const OffsetType &o) { return !o.IsZero(); });

      DMSG("<- star");
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
    PseudoLinearSet(SetType<OffsetType> &&os, SetType<GeneratorType> &&gs)
        : offsets_(os), generators_(gs) {
      assert(Ok());
    }

    bool Ok() const {
      for (auto &g : generators_) {
        if (g.IsZero()) {
          DMSG("Zero vector in generators.");
          return false;
        }
      }
      return true;
    }

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
        offset1_iter = offsets_.MarkErased(offset1_iter);

        bool necessary = true;
        for (auto &offset2 : offsets_) {
          auto new_offset = offset1 - offset2;
          if (new_offset.IsValid() && simplifier.IsCovered(new_offset)) {
            necessary = false;
            break;
          }
        }

        if (necessary) {
          offsets_.AbortErase();
          DMSG("AbortErase: " << offset1);
        } else {
          offsets_.CommitErase();
          DMSG("CommitErase: " << offset1);
        }
      }

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
    }

    std::set<Var> getVariables() const {
      std::set<Var> res;
      for (const auto &offset : offsets_) {
        auto ovars = offset.getVariables();
        res.insert(ovars.begin(), ovars.end());
      }
      for (const auto &g : generators_) {
        auto gvars = g.getVariables();
        res.insert(gvars.begin(),gvars.end());
      }
      return res;
    }

    SetType<OffsetType> offsets_;
    SetType<GeneratorType> generators_;

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
