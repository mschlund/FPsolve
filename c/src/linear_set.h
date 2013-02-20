#include <cassert>
#include <set>

/*
 * Note that this should be std::set since we use the fact that inserting
 * doesn't invalidate any iterators and erase doesn't invalidate any iterators
 * except for the erased one.
 */
template <typename V>
using SparseVecSet = std::set< SparseVec<V> >;

template <typename Simplifier, typename V>
class LinearSet {
  public:
    LinearSet() = default;
    LinearSet(const SparseVec<V> &o, const SparseVecSet<V> &vs)
        : offset_(o), generators_(vs) {}
    LinearSet(SparseVec<V> &&o, SparseVecSet<V> &&vs)
        : offset_(std::move(o)), generators_(std::move(vs)) {}

    LinearSet(const SparseVec<V> &v) : offset_(v) {}
    LinearSet(SparseVec<V> &&v) : offset_(std::move(v)) {}

    // FIXME: do we need this?
    // LinearSet& operator=(const LinearSet &s) = default;
    // LinearSet& operator=(LinearSet &&s) = default;

    ~LinearSet() = default;

    LinearSet operator+(const LinearSet &rhs) const {
      SparseVec<V> result_offset{offset_ + rhs.offset_};

      SparseVecSet<V> result_generators;
      std::insert_iterator< SparseVecSet<V> >
        inserter(result_generators, result_generators.begin());

      std::set_union(generators_.begin(), generators_.end(),
                     rhs.generators_.begin(), rhs.generators_.end(), inserter);


      // FIXME: disable the whole following loop if the simplifier is dummy...

      /* This is a bit tricky.  We use here the fact that erase will return the
       * iterator to the next element (i.e., it's erase that's advancing iter).
       * Morevore std::set has the property that erase(iter) only invalidates
       * iter and insert doesn't invalidate any iterators. */
      if (simplifier_.IsActive()) {
        for (auto iter = result_generators.begin();
             iter != result_generators.end(); ) {
          auto vec = *iter;
          iter = result_generators.erase(iter);
          if (!simplifier_.IsCovered(vec, result_generators)) {
            result_generators.insert(std::move(vec));
          }
        }
      }

      return LinearSet{std::move(result_offset), std::move(result_generators)};
    }

    friend std::ostream& operator<<(std::ostream &out, const LinearSet lset) {
      out << lset.offset_;
      out << " : ";
      out << "{ ";
      for (const auto &v : lset.generators_) {
        out << v << " ";
      }
      out << " }";
      return out;
    }

  private:
    SparseVec<V> offset_;
    SparseVecSet<V> generators_;
    static Simplifier simplifier_;
};

template <typename Simplifier, typename V>
Simplifier LinearSet<Simplifier, V>::simplifier_;

class DummySimplifier {
  public:
    bool IsActive() const { return false; }
    template <typename V>
    bool IsCovered(const SparseVec<V> &lhs, const SparseVecSet<V> &rhs_set) {
      return false;
    }
};

