#pragma once

#include <set>
#include <sstream>

#include "semilinear_util.h"
#include "sparse_vec.h"

#include "debug_output.h"


/* TODO:
 * - The simplification could be a bit better -- we currently only check if
 *   something (offset or generator) can be is included in the (remaining)
 *   generators.  What we could improve is to use the offsets in this
 *   calculation.  The challenge is that the current code for simplifiers
 *   assumes every element of the set can be used arbitrary many times, whereas
 *   we can only use offset once.
 */


/*
 * An abstraction of semilinear sets that is quite similar to a linear set, but
 * with multiple offsets.  This allows to over-approximate the semilinear set.
 */
template <typename Simpl, typename Var>
class PseudoLinearSet : public Semiring< PseudoLinearSet<Simpl, Var> > {
  typedef SparseVec<Var> SparseVec_;
  public:
    PseudoLinearSet() = default;

    PseudoLinearSet(const Var &v, const Counter &c) {
      offsets_.insert(SparseVec_(v, c));
    }

    PseudoLinearSet(std::initializer_list<SparseVec_> list) {
      for (auto &vec : list) {
        offsets_.insert(vec);
      }
    }

    PseudoLinearSet(const PseudoLinearSet &s) = default;
    PseudoLinearSet(PseudoLinearSet &&s) = default;

    template <typename Simpl2>
    PseudoLinearSet(const LinearSet<Simpl2, Var> &lset)
        : generators_(lset.GetGenerators()) {
      offsets_.insert(lset.GetOffset());

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
      std::set<SparseVec_> result_offsets;
      std::set<SparseVec_> result_generators;

      std::set_union(offsets_.begin(), offsets_.end(),
                     rhs.offsets_.begin(), rhs.offsets_.end(),
                     std::inserter(result_offsets, result_offsets.begin()));

      std::set_union(generators_.begin(), generators_.end(),
                     rhs.generators_.begin(), rhs.generators_.end(),
                     std::inserter(result_generators, result_generators.begin()));

      offsets_ = std::move(result_offsets);
      generators_ = std::move(result_generators);

      Simplify();

      return *this;
    }

    PseudoLinearSet& operator*=(const PseudoLinearSet &rhs) {
      std::set<SparseVec_> result_offsets;
      std::set<SparseVec_> result_generators;

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

      Simplify();

      return *this;
    }

    PseudoLinearSet star() const {

      std::set<SparseVec_> result_offsets = { SparseVec_{} };
      std::set<SparseVec_> result_generators;

      std::set_union(offsets_.begin(), offsets_.end(),
                     generators_.begin(), generators_.end(),
                     std::inserter(result_generators, result_generators.begin()));

      return PseudoLinearSet{ std::move(result_offsets),
                              std::move(result_generators) };
    }


    std::string string() const {
      std::stringstream ss;
      ss << "{";
      for (auto &vec : offsets_) {
        ss << " " << vec;
      }
      ss << " :";
      for (auto &vec : generators_) {
        ss << " " << vec;
      }
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
    PseudoLinearSet(std::set<SparseVec_> &&os, std::set<SparseVec_> &&gs)
        : offsets_(os), generators_(gs) {}

    void Simplify() {
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
      SimplifySet(simplifier_, generators_);

      // TODO: We could try a bit smarter approach by using IsCovered with an
      // additional set of elements that can be used only once (i.e., offsets)

      /* Simplify offsets.  This uses the same trick as SimplifySet -- erase in
       * std::set does not invalidate any iterators (except for the removed). */
      for (auto offset1_iter = offsets_.begin(); offset1_iter != offsets_.end(); ) {
        auto offset1 = std::move(*offset1_iter);
        /* Erase automatically advances the iterator to the next element. */
        offset1_iter = offsets_.erase(offset1_iter);

        bool necessary = true;
        for (auto &offset2 : offsets_) {
          auto new_offset = offset1 - offset2;
          if (new_offset.IsValid() &&
              simplifier_.IsCovered(new_offset, generators_)) {
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

    std::set<SparseVec_> offsets_;
    std::set<SparseVec_> generators_;

    // FIXME: Should that be static, pointer or just object???
    static Simpl simplifier_;

    template <typename Simpl1, typename Simpl2, typename Simpl3, typename V>
    friend PseudoLinearSet<Simpl1, V> SemilinearToPseudoLinear(
        const SemilinearSet<Simpl2, Simpl3, V> &semilinear);
};

template <typename Simpl, typename Var>
Simpl PseudoLinearSet<Simpl, Var>::simplifier_;


template <typename Simpl1, typename Simpl2, typename Simpl3, typename Var>
PseudoLinearSet<Simpl1, Var> SemilinearToPseudoLinear(
    const SemilinearSet<Simpl2, Simpl3, Var> &semilinear) {
  PseudoLinearSet<Simpl1, Var> pseudo_lset;
  semilinear.Iterate([&](const LinearSet<Simpl3, Var> &lset) -> void {
    pseudo_lset += PseudoLinearSet<Simpl1, Var>{lset};
  });
  return pseudo_lset;
}

typedef PseudoLinearSet<SparseVecSimplifier<VarPtr>, VarPtr> DefaultPseudoLinearSet;
