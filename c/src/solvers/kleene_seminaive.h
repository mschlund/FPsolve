/*
 * kleene_seminaive.h
 *
 *  Created on: 31.03.2014
 *      Author: schlund
 */

#ifndef KLEENE_SEMINAIVE_H_
#define KLEENE_SEMINAIVE_H_

// Polynomial is parametrized by a semiring
#define POLY_TYPE template <typename> class


template <typename SR, POLY_TYPE Poly>
class Kleene {
  typedef std::vector< std::pair< VarId, Poly<SR> > > GenericEquations;
  public:
  ValuationMap<SR> solve_fixpoint(const GenericEquations& equations, int max_iter) {
    std::vector<Poly<SR>> F;
    std::vector<VarId> poly_vars;
    for (const auto &eq : equations) {
      poly_vars.push_back(eq.first);
      F.push_back(eq.second);
    }
    Matrix<SR> result = solve_fixpoint(F, poly_vars, max_iter);

    // repack everything and return it
    ValuationMap<SR> solution;
    auto result_vec = result.getElements();
    assert(result_vec.size() == poly_vars.size());
    for (std::size_t i = 0; i < result_vec.size(); ++i) {
      solution.insert(std::make_pair(poly_vars[i], result_vec[i]));
    }
    return solution;
  }

  /*
   * Init:
   * previous_values = 0
   * values = F(0)
   * update = F(0)
   *
   * Each iteration consists of :
   *
   * update = polynomials.HeightUnfoldingAt(values, update, previous_values)
   * previous_values = values
   * values = values + update
  */
  Matrix<SR> solve_fixpoint(
    const std::vector< Poly<SR> > &polynomials,
    const std::vector<VarId> &variables, std::size_t max_iter) {

    assert(polynomials.size() == variables.size());

    Matrix< Poly<SR> > F_mat{polynomials.size(), polynomials};

    Matrix<SR> values{polynomials.size(), 1};

    ValuationMap<SR> values;
    for (auto &variable : variables) {
      values.insert(std::make_pair(variable, SR::null()));
    }

    for (unsigned int i=0; i < max_iter; ++i) {

    }

    return values;
  }
};



#endif /* KLEENE_SEMINAIVE_H_ */
