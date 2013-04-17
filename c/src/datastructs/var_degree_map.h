#pragma once

#include <iosfwd>
#include <limits>
#include <map>

#include "var.h"

typedef std::uint_fast16_t Degree;

class VarDegreeMap {
  typedef std::map<VarId, Degree> Map_;

  public:
    typedef Map_::iterator iterator;
    typedef Map_::const_iterator const_iterator;

    VarDegreeMap() = default;
    VarDegreeMap(const VarDegreeMap &) = default;
    VarDegreeMap(VarDegreeMap &&) = default;
    VarDegreeMap& operator=(const VarDegreeMap &) = default;
    VarDegreeMap& operator=(VarDegreeMap &&) = default;

    iterator begin() { return map_.begin(); }
    iterator end() { return map_.end(); }
    const_iterator begin() const { return map_.begin(); }
    const_iterator end() const { return map_.end(); }

    iterator find(const VarId v) { return map_.find(v); }
    const_iterator find(const VarId v) const { return map_.find(v); }

    void clear() { map_.clear(); };
    bool empty() const { return map_.empty(); };

    bool operator<(const VarDegreeMap &rhs) const { return map_ < rhs.map_; }
    bool operator==(const VarDegreeMap &rhs) const { return map_ == rhs.map_; }

    /* Returns 0 if var is not in the map. */
    Degree GetDegreeOf(const VarId var) const;

    void Insert(const VarId var, Degree deg = 1);

    /* The VarId must be in the map. */
    void Erase(const VarId var, Degree deg = 1);
    void EraseAll(const VarId var) {
      Erase(var, std::numeric_limits<Degree>::max());
    }

    void Merge(const VarDegreeMap &to_merge);


  private:

    bool SanityCheck() const {
      for (auto &var_degree : map_) {
        if (var_degree.second == 0) {
          return false;
        }
      }
      return true;
    }


    Map_ map_;
};

std::ostream& operator<<(std::ostream &out, const VarDegreeMap &map);
