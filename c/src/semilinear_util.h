#pragma once

template <typename Simplifier, typename Elem>
void SimplifySet(Simplifier &simplifier, std::set<Elem> &set) {

  /* This is a bit tricky.  We use here the fact that erase will return the
   * iterator to the next element (i.e., it's erase that's advancing iter).
   * Morevore std::set has the property that erase(iter) only invalidates
   * iter and insert doesn't invalidate any iterators. */
  if (simplifier.IsActive() && 1 < set.size()) {
    for (auto iter = set.begin(); iter != set.end(); ) {
      auto tmp_elem = std::move(*iter);
      /* Erase automatically advances the iterator to the next element. */
      iter = set.erase(iter);
      /* Add it back only if it's not "covered" by the set. */
      if (!simplifier.IsCovered(tmp_elem, set)) {
        // FIXME: GCC 4.7 does not have emplace
        set.insert(std::move(tmp_elem));
      }
    }
  }
}
