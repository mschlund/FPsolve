#pragma once

#include "../debug_output.h"
#include "../datastructs/vec_set.h"

template <typename Simpl, typename Elem>
void SimplifySet(VecSet<Elem> &to_simpl) {


  if (!Simpl::IsActive() || to_simpl.size() <= 1) {
    return;
  }

#ifdef DEBUG_OUTPUT
  DMSG("SimplifySet: {");
  for (auto &x : to_simpl) {
    DMSG(x);
  }
  DMSG("}");
#endif

  /* This is a bit tricky.  We use here the fact that erase will return the
   * iterator to the next element (i.e., it's erase that's advancing iter).
   * Morevore std::set has the property that erase(iter) only invalidates iter
   * and insert doesn't invalidate any iterators. */

  for (auto iter = to_simpl.begin(); iter != to_simpl.end(); ) {
    Elem tmp_elem = *iter;
    /* Erase automatically advances the iterator to the next element. */
    iter = to_simpl.MarkErased(iter);

    /* Now we check if it's necessary to add it back. */

    Simpl simplifier{to_simpl};

    if (!simplifier.IsCovered(tmp_elem)) {
      DMSG("Cannot remove " << tmp_elem);
      // FIXME: GCC 4.7 does not have emplace
      to_simpl.AbortErase();
    } else {
      DMSG("Removing " << tmp_elem);
      to_simpl.CommitErase();
    }
  }

}
