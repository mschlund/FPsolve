#include "polynomial.h"

std::uint16_t GetDegreeOf(const std::map<VarPtr, std::uint16_t> &map,
    const VarPtr var) {
  auto var_iter = map.find(var);
  if (var_iter != map.end()) {
    return var_iter->second;
  } else {
    return 0;
  }
}

void EraseAll(std::map<VarPtr, std::uint16_t> &map, const VarPtr var) {
  map.erase(var);
}


void Insert(std::map<VarPtr, std::uint16_t> &map, const VarPtr var,
    std::uint16_t deg) {
  auto iter_inserted = map.emplace(var, deg);
  if (!iter_inserted.second) {
    iter_inserted.first->second += deg;
  }
}

void Erase(std::map<VarPtr, std::uint16_t> &map, const VarPtr var,
    std::uint16_t deg) {
  auto var_iter = map.find(var);
  assert(var_iter != map.end());
  if (deg >= var_iter->second) {
    map.erase(var_iter);
  } else {
    var_iter->second -= deg;
  }
}
