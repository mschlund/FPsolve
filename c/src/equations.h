#pragma once


template <typename A>
using Equations = std::vector< std::pair< VarPtr, Polynomial<A> > >;

template <typename A, typename F>
auto MapEquations(const Equations<A> &equations, F fun)
    -> Equations<typename std::result_of<F(A)>::type> {

  Equations<typename std::result_of<F(A)>::type> new_equations;

  for (auto &var_poly : equations) {
    new_equations.emplace_back(var_poly.first,
      var_poly.second.Map([&fun](const SemilinSetExp &slset) {
        return fun(slset);
      })
    );
  }

  return new_equations;
}
