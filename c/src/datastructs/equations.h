#pragma once

#include <vector>

#include "../polynomials/polynomial.h"
#include "../polynomials/non_commutative_polynomial.h"
#include "../semirings/free-semiring.h"

template <typename A>
using Equations = std::vector< std::pair< VarId, Polynomial<A> > >;

template <typename A>
using NCEquations = std::vector< std::pair< VarId, NonCommutativePolynomial<A> > >;


template <typename A, typename F>
auto MapEquations(const Equations<A> &equations, F fun)
    -> Equations<typename std::result_of<F(A)>::type> {

  Equations<typename std::result_of<F(A)>::type> new_equations;

  for (auto &var_poly : equations) {
    new_equations.emplace_back(var_poly.first,
      var_poly.second.Map([&fun](const A &coeff) {
        return fun(coeff);
      })
    );
  }

  return new_equations;
}

template <typename A, typename F>
auto MakeCommEquationsAndMap(const NCEquations<A> &equations, F fun)
    -> Equations<typename std::result_of<F(A)>::type> {

  Equations<typename std::result_of<F(A)>::type> new_equations;

  for (auto &var_poly : equations) {
    //std::cout << "COMM-POLY: " << var_poly.second.make_commutative().string() << std::endl;
    new_equations.emplace_back(var_poly.first,
      var_poly.second.make_commutative().Map([&fun](const A &coeff) {
        return fun(coeff);
      })
    );
  }

  return new_equations;
}

