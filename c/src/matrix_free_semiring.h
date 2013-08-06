#pragma once

#include "semirings/free-semiring.h"
#include "datastructs/matrix.h"


template <typename SR>
Matrix<SR> FreeSemiringMatrixEval(const Matrix<FreeSemiring> &matrix,
    const ValuationMap<SR> &valuation) {

  const std::vector<FreeSemiring> &elements = matrix.getElements();
  std::vector<SR> result;

  /* We have a single Evaluator that is used for all evaluations of the elements
   * of the original matrix.  This way if different elements refer to the same
   * FreeSemiring subexpression, we memoize the result and reuse it. */
  Evaluator<SR> evaluator{valuation};

  for(auto &elem : elements) {
    result.emplace_back(elem.Eval(evaluator));
  }

  return Matrix<SR>(matrix.getRows(), std::move(result));
}

/* FIXME: Temporary wrapper for compatibility with the old implementation. */
template <typename SR>
SR FreeSemiring_eval(FreeSemiring elem,
    ValuationMap<SR> *valuation) {
  return elem.Eval(*valuation);
}

/* FIXME: Temporary wrapper for compatibility with the old implementation. */
template <typename SR>
Matrix<SR> FreeSemiring_eval(Matrix<FreeSemiring> matrix,
    ValuationMap<SR> *valuation) {
  return FreeSemiringMatrixEval(matrix, *valuation);
}
