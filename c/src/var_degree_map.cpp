#include <cassert>
#include <iostream>

#include "var_degree_map.h"

Degree VarDegreeMap::GetDegreeOf(const VarPtr var) const {
  auto var_iter = map_.find(var);
  if (var_iter != map_.end()) {
    return var_iter->second;
  } else {
    return 0;
  }
}

void VarDegreeMap::Insert(const VarPtr var, std::uint16_t deg) {
  if (deg == 0) {
    return;
  }
  // GCC 4.7 is missing emplace
  // auto iter_inserted = map_.emplace(var, deg);
  auto iter_inserted = map_.insert(std::make_pair(var, deg));
  if (!iter_inserted.second) {
    iter_inserted.first->second += deg;
  }
  assert(SanityCheck());
}

void VarDegreeMap::Erase(const VarPtr var, std::uint16_t deg) {
  auto var_iter = map_.find(var);
  assert(var_iter != map_.end());
  if (deg >= var_iter->second) {
    map_.erase(var_iter);
  } else {
    var_iter->second -= deg;
  }
  assert(SanityCheck());
}

void VarDegreeMap::Merge(const VarDegreeMap &to_merge) {
  for (auto &var_degree : to_merge) {
    auto iter_inserted = map_.insert(var_degree);
    if (!iter_inserted.second) {
      iter_inserted.first->second = std::max(iter_inserted.first->second,
                                             var_degree.second);
    }
  }
  assert(SanityCheck());
}

std::ostream& operator<<(std::ostream &out, const VarDegreeMap &map) {
  for (const auto &pair : map) {
    out << "[" << pair.first << " |-> " << pair.second << "]";
  }
  return out;
}
