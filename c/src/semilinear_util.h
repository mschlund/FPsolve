#pragma once

#include "debug_output.h"

/*
 * Try removing elements from the first set that are included in the second one.
 * If you want to simplify just one set, it is ok to pass the same std::set for
 * both arguments.
 */
template <typename Simplifier, typename Elem>
void SimplifySet(Simplifier &simplifier, std::set<Elem> &to_simpl,
    const std::set<Elem> &available) {

  /* This is a bit tricky.  We use here the fact that erase will return the
   * iterator to the next element (i.e., it's erase that's advancing iter).
   * Morevore std::set has the property that erase(iter) only invalidates
   * iter and insert doesn't invalidate any iterators. */
  if (simplifier.IsActive() && 1 < to_simpl.size()) {

#ifdef DEBUG_OUTPUT
    DMSG("Lhs: {");
    for (auto &x : to_simpl) {
      DMSG(x);
    }
    DMSG("}");
    if (to_simpl != available) {
      DMSG("Rhs: {");
      for (auto &x : available) {
        DMSG(x);
      }
      DMSG("}");
    }
#endif

    for (auto iter = to_simpl.begin(); iter != to_simpl.end(); ) {
      Elem tmp_elem = std::move(*iter);
      /* Erase automatically advances the iterator to the next element. */
     iter = to_simpl.erase(iter);
      // Add it back only if it's not "covered" by the set.
      if (!simplifier.IsCovered(tmp_elem, available)) {
        DMSG("Cannot remove " << tmp_elem);
        // FIXME: GCC 4.7 does not have emplace
        to_simpl.insert(std::move(tmp_elem));
      } else {
        DMSG("Removing " << tmp_elem);
      }
    }
  }

}

template <typename Simplifier, typename Elem>
void SimplifySet(Simplifier &simplifier, std::set<Elem> &to_simpl) {
  SimplifySet(simplifier, to_simpl, to_simpl);
}
