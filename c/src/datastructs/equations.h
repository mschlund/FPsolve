#pragma once

#include <vector>

#include "../polynomials/commutative_polynomial.h"
#include "../polynomials/non_commutative_polynomial.h"
#include "../semirings/free-semiring.h"

template<template <typename> class Poly, typename SR>
using GenericEquations = std::vector< std::pair< VarId, Poly<SR> > > ;

template <typename A>
using Equations = GenericEquations<CommutativePolynomial, A>;

template <typename A>
using NCEquations = GenericEquations<NonCommutativePolynomial, A>;

template <typename A>
using NCEquationsBase = GenericEquations<NonCommutativePolynomialBase, A>;



template <template <typename> class Poly, typename SR, typename F>
auto MapEquations(const GenericEquations<Poly, SR> &equations, F fun)
    -> GenericEquations<Poly, typename std::result_of<F(SR)>::type> {

  GenericEquations<Poly, typename std::result_of<F(SR)>::type> new_equations;

  for (auto &var_poly : equations) {
    new_equations.emplace_back(var_poly.first,
      var_poly.second.Map([&fun](const SR &coeff) {
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


template <typename A>
auto MakeCommEquations(const NCEquations<A> &equations)
    -> Equations<A> {

  Equations<A> new_equations;

  for (auto &var_poly : equations) {
    //std::cout << "COMM-POLY: " << var_poly.second.make_commutative().string() << std::endl;
    new_equations.emplace_back(var_poly.first,
      var_poly.second.make_commutative());
  }

  return new_equations;
}
