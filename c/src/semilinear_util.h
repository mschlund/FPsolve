#pragma once

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

    std::cout << "-- SimplifySet" << std::endl;
    for (auto &x : to_simpl) {
      std::cout << "Lhs: " << x << std::endl;
    }
    if (to_simpl != available) {
      for (auto &x : available) {
        std::cout << "Rhs: " << x << std::endl;
      }
    }

    for (auto iter = to_simpl.begin(); iter != to_simpl.end(); ) {
      Elem tmp_elem = std::move(*iter);
      /* Erase automatically advances the iterator to the next element. */
      iter = to_simpl.erase(iter);
      /* Add it back only if it's not "covered" by the set. */
      if (!simplifier.IsCovered(tmp_elem, available)) {
        // FIXME: GCC 4.7 does not have emplace
        to_simpl.insert(std::move(tmp_elem));
      } else {
        std::cout << "SimplifySet: Removed " << tmp_elem << std::endl;
      }
    }
  }

}

template <typename Simplifier, typename Elem>
void SimplifySet(Simplifier &simplifier, std::set<Elem> &to_simpl) {
  SimplifySet(simplifier, to_simpl, to_simpl);
}
