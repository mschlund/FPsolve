#pragma once

#include <sstream>
#include <iostream>

#include <algorithm>
#include <string>

#include <unordered_map>
#include <map>
#include <set>




template <typename A>
std::string ToString(const A &a) {
  std::stringstream ss;
  ss << a;
  return ss.str();
}

template <typename A, typename B>
std::ostream& operator<<(std::ostream &out, const std::pair<A, B> &pair) {
  //return out << "(" << pair.first << "," << pair.second << ")";
  return out  << pair.first << ":" << pair.second;
}


template <typename A, typename B>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<A,B>& values) {
  for (auto value = values.begin(); value != values.end(); ++value) {
    os << value->first << "→" << value->second << ";";
  }
  return os;
}

template <typename A, typename B>
std::ostream& operator<<(std::ostream& os, const std::map<A,B>& values) {
  for (auto value = values.begin(); value != values.end(); ++value) {
    os << value->first << "→" << value->second << ";";
  }
  return os;
}

template <typename A>
std::ostream& operator<<(std::ostream& os, const std::set<A>& values) {
  os << "{" << ToStringSorted(values, ",") << "}";
  return os;
}

/* Print contents of a container in a sorted order using the given separator.
 * The sorting is done at the level of strings. */
template <typename Container>
std::string ToStringSorted(const Container &container,
                           const std::string &sep = "") {
  std::vector<std::string> strings;
  for (auto &x : container) {
    strings.emplace_back(ToString(x));
  }
  std::sort(strings.begin(), strings.end());
  std::stringstream ss;
  for (auto iter = strings.begin(); iter != strings.end(); ++iter) {
    if (iter != strings.begin()) {
      ss << sep;
    }
    ss << *iter;
  }
  return ss.str();
}
