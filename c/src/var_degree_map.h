#pragma once

#include <iosfwd>
#include <limits>
#include <map>

#include "var.h"

/*
 * TODO:
 * - The std::map<VarPtr, Degree> should be a separate class, i.e., a thin
 *   wrapper around the std::map<VarPtr, Degree>...
 */

typedef std::uint_fast16_t Degree;

/*
Degree GetDegreeOf(const std::map<VarPtr, Degree> &map,
    const VarPtr var);

void EraseAll(std::map<VarPtr, Degree> &map, const VarPtr var);

void Insert(std::map<VarPtr, Degree> &map, const VarPtr var, Degree deg = 1);

void Erase(std::map<VarPtr, Degree> &map, const VarPtr var, Degree deg = 1);

void Merge(std::map<VarPtr, Degree> &to_modify,
    const std::map<VarPtr, Degree> &to_merge);

std::ostream& operator<<(std::ostream &out, const VarDegreeMap &map);
*/

class VarDegreeMap {
  typedef std::map<VarPtr, Degree> Map_;

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

    iterator find(const VarPtr &v) { return map_.find(v); }
    const_iterator find(const VarPtr &v) const { return map_.find(v); }

    void clear() { map_.clear(); };
    bool empty() const { return map_.empty(); };

    bool operator<(const VarDegreeMap &rhs) const { return map_ < rhs.map_; }
    bool operator==(const VarDegreeMap &rhs) const { return map_ == rhs.map_; }

    /* The VarPtr must be in the map. */
    Degree GetDegreeOf(const VarPtr var) const;

    void Insert(const VarPtr var, Degree deg = 1);

    /* The VarPtr must be in the map. */
    void Erase(const VarPtr var, Degree deg = 1);
    void EraseAll(const VarPtr var) {
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
