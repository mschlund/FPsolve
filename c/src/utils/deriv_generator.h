#ifndef DERIV_GENERATOR_H
#define DERIV_GENERATOR_H

/* This defines the generator that is able to create all possible combinations
 * of integers such they are smaller than max and their sum is between min_sum
 * and max_sum. */
class Generator {
  public:
    Generator(std::unordered_map<VarId, Degree> &max, Degree min_sum, Degree max_sum)
        : val_map_(max), max_(max), current_sum_(0),
          min_sum_(min_sum), max_sum_(max_sum) {
      for (auto &var_val : val_map_) {
        var_val.second = 0;
      }
      assert(0 < max_.size());
      assert(min_sum <= max_sum);
      assert(current_sum_ == CurrentSum());
      // FIXME: add assertion that sum of max >= min_sum
    }

    bool NextCombination() {
      bool added = false;
      bool valid = false;

      do {
        if (current_sum_ < min_sum_) {
          added = JumpMin();
        } else if (max_sum_ < current_sum_) {
          added = JumpMax();
        } else {
          added = AddOne();
        }
        assert(current_sum_ == CurrentSum() && BelowEqualMax());
        valid = min_sum_ <= current_sum_ && current_sum_ <= max_sum_;
        if (added && valid) {
          return true;
        }
      } while (added && !valid);

      return false;
    }

    const std::unordered_map<VarId, Degree>& GetMap() const {
      assert(current_sum_ == CurrentSum());
      assert(BelowEqualMax());
      assert(min_sum_ <= current_sum_ && current_sum_ <= max_sum_);

      return val_map_;
    }

  private:
    Degree Max(const VarId v) const {
      auto lookup = max_.find(v);
      if (lookup == max_.end()) {
        assert(false);
        return 0;
      }
      return lookup->second;
    }

    Degree CurrentSum() const {
      Degree result = 0;
      for (auto &var_val : val_map_) {
        result += var_val.second;
      }
      return result;
    }

    bool BelowEqualMax() const {
      for (auto &var_val : val_map_) {
        if (var_val.second > Max(var_val.first)) {
          return false;
        }
      }
      return true;
    }

    /* Add 1 to the current vector, wrap-around if some value is > max_.
     * Returns false if we cannot add 1 (i.e., the last element would
     * overflow). */
    bool AddOne() {
      return AddOne(val_map_.begin());
    }

    bool AddOne(typename std::unordered_map<VarId, Degree>::iterator start) {
      for (auto iter = start; iter != val_map_.end(); ++iter) {
        auto &integer = iter->second;
        const auto integer_max = Max(iter->first);

        if (integer_max == 0) {
          continue;
        }

        ++integer;
        ++current_sum_;
        if (integer <= integer_max) {
          return true;
        }
        current_sum_ -= integer;
        integer = 0;
      }
      return false;
    }

    /* Create the smallest vector that satisfies the min_sum_ requirement.  This
     * means adding (without overflowing) min_sum_ - current_sum_. */
    bool JumpMin() {
      assert(min_sum_ > current_sum_);
      auto remaining = min_sum_ - current_sum_;
      for (auto &elem : val_map_) {
        if (Max(elem.first) == 0) {
          continue;
        }
        auto &integer = elem.second;
        const auto integer_max = Max(elem.first);
        auto to_add = integer + remaining > integer_max ?
                      integer_max - integer : remaining;
        integer += to_add;
        current_sum_ += to_add;
        remaining -= to_add;
        if (remaining == 0) {
          return true;
        }
      }
      /* Added as much as we could, but still not enough... */
      return false;
    }

    /* Create the smallest vector that satisfies the max_sum_ requirement.  This
     * means that we add (with overflow) enough to get below max. */
    bool JumpMax() {
      assert(max_sum_ < current_sum_);
      for (auto iter = val_map_.begin(); iter != val_map_.end(); ++iter) {
        if (Max(iter->first) == 0) {
          continue;
        }

        auto &integer = iter->second;

        /* If integer == 0 then this is harmless. */
        current_sum_ -= integer;
        integer = 0;

        /* We're wrapping around integer and should add 1 to the next position.
         * Check if that is enough or whether we should try to wrap-around the
         * next position too.  This can happen when integer == 1. */
        if (current_sum_ - integer + 1 <= max_sum_) {
          ++iter;
          return AddOne(iter);
        }
      }
      return false;
    }

    std::unordered_map<VarId, Degree> val_map_;
    const std::unordered_map<VarId, Degree> &max_;
    Degree current_sum_;
    Degree min_sum_;
    Degree max_sum_;
};

#endif  /* DERIV_GENERATOR_H */
