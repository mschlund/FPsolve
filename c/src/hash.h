#pragma once

#include <functional>
#include <vector>

/* Taken from Boost. */
template <typename A>
inline void HashCombine(std::size_t &seed, const A &a) {
  std::hash<A> hash_value;
  seed ^= hash_value(a) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


namespace std {

template<typename A, typename B>
struct hash< std::pair<A, B> > {
  inline std::size_t operator()(const std::pair<A, B> &pair) const {
    std::size_t h = 0;
    HashCombine(h, pair.first);
    HashCombine(h, pair.second);
    return h;
  }

};

template<typename A>
struct hash< std::vector<A> > {
  inline std::size_t operator()(const std::vector<A> &vec) const {
    std::size_t h = 0;
    for (auto &x : vec) {
      HashCombine(h, x);
    }
    return h;
  }
};

template<typename A>
struct hash< std::set<A> > {
  inline std::size_t operator()(const std::set<A> &set) const {
    std::size_t h = 0;
    for (auto &x : set) {
      HashCombine(h, x);
    }
    return h;
  }
};


}  /* namespace std */
