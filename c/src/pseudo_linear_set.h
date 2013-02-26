#pragma once

#include "sparse_vec.h"
#include "semilinear_util.h"

/*
 * An abstraction of semilinear sets that is quite similar to a linear set, but
 * with multiple offsets.  This allows to over-approximate the semilinear set.
 */
template <typename Simpl, typename Var>
class PseudoLinearSet {
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

    PseudoLinearSet& operator=(const PseudoLinearSet &s) = default;
    PseudoLinearSet& operator=(PseudoLinearSet &&s) = default;

    ~PseudoLinearSet() = default;

    static PseudoLinearSet null() { return PseudoLinearSet{}; }
    static PseudoLinearSet one() { return PseudoLinearSet{{}}; }

    bool operator==(const PseudoLinearSet &rhs) {
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

      SimplifySet(simplifier_, result_offsets);
      SimplifySet(simplifier_, result_generators);

      offsets_ = std::move(result_offsets);
      generators_ = std::move(result_generators);

      return *this;
    }

    PseudoLinearSet& operator*=(const PseudoLinearSet &rhs) {
      // std::cout << "-> SemilinearSet::operator*" << std::endl;
      // std::cout << *this << std::endl;
      // std::cout << rhs << std::endl;
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


      SimplifySet(simplifier_, result_offsets);
      SimplifySet(simplifier_, result_generators);

      offsets_ = std::move(result_offsets);
      generators_ = std::move(result_generators);

      // std::cout << *this << std::endl;
      // std::cout << "<- SemilinearSet::operator*" << std::endl;
      return *this;
    }

    PseudoLinearSet operator+(const PseudoLinearSet &rhs) const {
      PseudoLinearSet tmp = *this;
      return tmp += rhs;
    }

    PseudoLinearSet operator*(const PseudoLinearSet &rhs) const {
      PseudoLinearSet tmp = *this;
      return tmp *= rhs;
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
      ss << " : ";
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

  private:
    PseudoLinearSet(std::set<SparseVec_> &&os, std::set<SparseVec_> &&gs)
        : offsets_(os), generators_(gs) {}


    std::set<SparseVec_> offsets_;
    std::set<SparseVec_> generators_;

    // FIXME: Should that be static, pointer or just object???
    static Simpl simplifier_;
};

template <typename Simpl, typename Var>
Simpl PseudoLinearSet<Simpl, Var>::simplifier_;


