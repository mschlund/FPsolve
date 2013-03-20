#pragma once

#include <functional>

/*
 * Wrapper around a pointer that can be used inside a map/unordered_map as a
 * key without incurring the copy of the object being pointed to.  So the copy
 * constructor and assignment operator only copy and assign the internal
 * pointer.  However, the ordering and equality uses the actual objects.
 * Similarly std::hash will hash the actual objects.
 */
template <typename A>
class KeyWrapper {
  public:
    explicit KeyWrapper(const A *a) : ptr_(a) {}
    explicit KeyWrapper(const A &a) : ptr_(&a) {}

    KeyWrapper(const KeyWrapper &k) = default;
    KeyWrapper(KeyWrapper &&k) = default;

    KeyWrapper& operator=(const KeyWrapper &k) = default;
    KeyWrapper& operator=(KeyWrapper &&k) = default;

    bool operator<(const KeyWrapper rhs) const { return *ptr_ < *rhs.ptr_; }
    bool operator==(const KeyWrapper &rhs) const {
      if (ptr_ == rhs.ptr_) {
        return true;
      }
      return *ptr_ == *rhs.ptr_;
    }

    A& Get() { return *ptr_; }
    const A& Get() const { return *ptr_; }

  private:
    const A *ptr_;
};

namespace std {

template<typename A>
struct hash< KeyWrapper<A> > {
  inline std::size_t operator()(const KeyWrapper<A> &k) const {
    std::hash<A> h;
    return h(k.Get());
  }
};

}  /* namespace std */
