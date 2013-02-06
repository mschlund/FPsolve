#pragma once

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

}  /* namespace std */
