/*
 * newton_generic.h
 *
 *  Created on: 11.04.2013
 *      Author: schlund
 */

#ifndef NEWTON_GENERIC_H_
#define NEWTON_GENERIC_H_

#include <algorithm>
#include <cstdint>

#include "datastructs/matrix.h"
#include "matrix_free_semiring.h"
#include "datastructs/var_degree_map.h"


#include "polynomials/polynomial.h"

#include "semirings/semiring.h"
#include "semirings/semilinear_set.h"
#include "semirings/float-semiring.h"
#include "semirings/free-semiring.h"


// Lin_Eq_Solver is parametrized by a semiring
#define LIN_EQ_SOLVER_TYPE template <typename> class

// DeltaGenerator is parametrized by a semiring
#define DELTA_GEN_TYPE template <typename> class



template <typename SR,
          LIN_EQ_SOLVER_TYPE LinEqSolverTemplate,
          DELTA_GEN_TYPE DeltaGeneratorTemplate>
class GenericNewton {
  typedef LinEqSolverTemplate<SR> LinEqSolver;
  typedef DeltaGeneratorTemplate<SR> DeltaGenerator;

  public:
  ValuationMap<SR> solve_fixpoint(const Equations<SR>& equations, int max_iter) {
    std::vector<Polynomial<SR>> F;
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
   * previous_newton_values = 0
   * newton_update = 0
   * delta = F(0)
   * newton_values = 0
   *
   * LinSolver = LinSolver(F)
   * DeltaGenerator = DeltaGenerator(F)
   *
   *Each iteration consists of :
   * newton_update = LinSolver.solve_linearization_at(newton_update,delta)
   * newton_values = (newton_values +) newton_update ( only if not itempotent))
   * if(!idempotent):
   *  delta = DeltaGenerator.delta_at(newton_update, previous_newton_values) (generates the "rhs" of linear system) //compute delta for the next iteration
   *
  */
  Matrix<SR> solve_fixpoint(
      const std::vector< Polynomial<SR> > &polynomials,
      const std::vector<VarId> &variables, std::size_t max_iter) {

    assert(polynomials.size() == variables.size());

    Matrix< Polynomial<SR> > F_mat{polynomials.size(), polynomials};

    Matrix<SR> newton_values{polynomials.size(), 1};
    Matrix<SR> previous_newton_values{polynomials.size(), 1};

    ValuationMap<SR> values;
    for (auto &variable : variables) {
      values.insert(std::make_pair(variable, SR::null()));
    }
    Matrix<SR> delta= Polynomial<SR>::eval(F_mat, values); //TODO: use delta_at here as well??

    Matrix<SR> newton_update{polynomials.size(), 1};

    LinEqSolver LinSolver = LinEqSolver(polynomials, variables);
    DeltaGenerator DeltaGen = DeltaGenerator(polynomials, variables);

    for (int i=0; i < max_iter; ++i) {
      newton_update = LinSolver.solve_lin_at(newton_values, delta, variables);

      /* No need to recompute delta if this is the last iteration... */
      if (!SR::IsIdempotent() && i + 1 < max_iter) {
        delta = DeltaGen.delta_at(newton_update, newton_values);
      }

      if (!SR::IsIdempotent()) {
        newton_values = newton_values + newton_update;
      } else {
        newton_values = newton_update;
      }

    }

    return newton_values;
  }
};

/*
 * TODO: Different LinSolvers:
 * 1) Comm-case: Compute star symbolically (F-W, recursive,...), store result, eval in each call to solve_linearization_at
 * 2) Approximate star by Kleene-iteration (dynamic number of iterations!)
 * 3) Comm-case: Newton-iteration on the Jacobian (a la Pivoteau et al)
 * 4) Faster algorithms that exploit idempotence?
 */

template <typename SR>
class CommutativeSymbolicLinSolver {
  public:
  // TODO: use initialization lists instead of temporaries+copying? see
  // http://www.parashift.com/c++-faq-lite/init-lists.html
  CommutativeSymbolicLinSolver(
      const std::vector< Polynomial<SR> >& F,
      const std::vector<VarId>& variables) {

    Matrix< Polynomial<SR> > jacobian = Polynomial<SR>::jacobian(F, variables);

    std::unordered_map<SR, VarId, SR> valuation_tmp;

    Matrix<FreeSemiring> jacobian_free = Polynomial<SR>::make_free(jacobian, &valuation_tmp);
    jacobian_star_ = new Matrix<FreeSemiring>(jacobian_free.star());

    // For benchmarking only ->
    std::cout << "Size of Jacobian: "
              << jacobian_star_->getRows()
              << " x "
              << jacobian_star_->getColumns()
              << std::endl;
    FreeSemiring::one().PrintStats();
    // <-

    //std::cout << "J*: " << *jacobian_star_ << std::endl;

    for (auto &pair : valuation_tmp) {
      valuation_.insert(std::make_pair(pair.second, pair.first));
      //std::cout << "valuation:  " << pair.second << ", " << pair.first << std::endl;
    }
  }

  virtual ~CommutativeSymbolicLinSolver(){
    delete jacobian_star_;
    jacobian_star_ = 0;
  }

  Matrix<SR> solve_lin_at(const Matrix<SR>& values, const Matrix<SR>& rhs,
                          const std::vector<VarId>& variables) {
    UpdateValuation(variables, values, valuation_);
    return FreeSemiringMatrixEval(*jacobian_star_, valuation_) * rhs;
  }

private:
  std::unordered_map<VarId, SR> valuation_;
  Matrix<FreeSemiring>* jacobian_star_;

  void UpdateValuation(const std::vector<VarId> &variables,
                       const Matrix<SR> &newton_values,
                       std::unordered_map<VarId, SR>& valuation) {
    //assert(!valuation.empty());
    assert(variables.size() == newton_values.getRows());
    assert(newton_values.getColumns() == 1);

    for (std::size_t i = 0; i < variables.size(); ++i) {
      valuation[variables[i]] = newton_values.At(i, 0);
    }
  }
};


template <typename SR>
class CommutativeConcreteLinSolver {
  public:
  CommutativeConcreteLinSolver(
      const std::vector< Polynomial<SR> >& F,
      const std::vector<VarId>& variables)
    : jacobian_(Polynomial<SR>::jacobian(F, variables)) {}

  Matrix<SR> solve_lin_at(const Matrix<SR>& values, const Matrix<SR>& rhs,
                          const std::vector<VarId>& variables) {
    assert(values.getColumns() == 1);

    assert(variables.size() == values.getRows());
    if (valuation_.size() > 0) {
      assert(valuation_.size() == variables.size());
    }

    for (std::size_t i = 0; i < variables.size(); ++i) {
      valuation_[variables[i]] = values.At(i, 0);
    }

    std::vector<SR> result_vec;
    for (auto &poly : jacobian_.getElements()) {
      result_vec.emplace_back(poly.eval(valuation_));
    }
    return Matrix<SR>{jacobian_.getRows(), std::move(result_vec)}.star()
           * rhs;
  }

  private:
    Matrix< Polynomial<SR> > jacobian_;
    /* We don't want to allocate the map every time, especially since the keys
     * do not change... */
    std::unordered_map<VarId, SR> valuation_;
};


template <typename SR>
class CommutativeDeltaGeneratorOld {
public:
  CommutativeDeltaGeneratorOld(
      const std::vector< Polynomial<SR> > &ps, const std::vector<VarId> &pvs)
      : polynomials(ps), poly_vars(pvs) {
    index_map_.reserve(poly_vars.size());
    for (std::size_t i = 0; i < poly_vars.size(); ++i) {
      index_map_.emplace(poly_vars[i], i);
      // GCC 4.7 doesn't have emplace for std::map
      zero_valuation_.insert(std::make_pair(poly_vars[i], SR::null()));
    }
  }

  Matrix<SR> delta_at(const Matrix<SR> &newton_update,
                      const Matrix<SR> &previous_newton_values) {

    assert(previous_newton_values.getColumns() == 1 && newton_update.getColumns() == 1);

    auto num_variables = poly_vars.size();
    assert(num_variables == previous_newton_values.getRows() &&
           num_variables == newton_update.getRows());

    std::vector<SR> delta_vector;

    std::unordered_map<VarId, Degree> current_max_degree;

    ValuationMap<SR> newton_update_map;

    for (std::size_t i = 0; i < num_variables; ++i) {
      current_valuation_[poly_vars[i]] = previous_newton_values.At(i, 0);
      newton_update_map.insert(
        std::make_pair(poly_vars[i], newton_update.At(i, 0)));
    }

    for (const auto &polynomial : polynomials) {
      SR delta = SR::null();
      Degree polynomial_max_degree = polynomial.get_degree();

      current_max_degree.clear();
      current_max_degree.insert(polynomial.GetVarDegreeMap().begin(),
                                polynomial.GetVarDegreeMap().end());

      if (polynomial_max_degree <= 1) {
        delta = SR::null();
      } else {
        /* We want to calculate all possible derivatives of at least second
         * order, but lower or equal to the degree of polynomial. */
        Generator generator{current_max_degree, 2, polynomial_max_degree};

        while (generator.NextCombination()) {

          SR prod = SR::one();

          for (const auto &var_val : generator.GetMap()) {
            assert(poly_vars[index_map_[var_val.first]] == var_val.first);
            for (int j = 0; j < var_val.second; ++j) {
              prod *= newton_update.At(index_map_[var_val.first], 0);
            }
          }

          SR polynomial_value =
            polynomial.DerivativeBinomAt(generator.GetMap(), current_valuation_);

          delta += polynomial_value * prod;
        }
      }
      delta_vector.emplace_back(std::move(delta));
    }
    return Matrix<SR>(delta_vector.size(), std::move(delta_vector));
  }

private:
  const std::vector< Polynomial<SR> > &polynomials;
  const std::vector<VarId> &poly_vars;
  std::unordered_map<VarId, std::size_t> index_map_;
  /* We cache the std::map so that we can avoid reallocating it every time.  And
   * we also make sure that we always overwrite everything before using it... */
  ValuationMap<SR> current_valuation_;
  ValuationMap<SR> zero_valuation_;
};




template <typename SR>
class CommutativeDeltaGenerator {
public:
  CommutativeDeltaGenerator(
      const std::vector< Polynomial<SR> > &ps, const std::vector<VarId> &pvs)
      : polynomials(ps), poly_vars(pvs) {
    index_map_.reserve(poly_vars.size());
    for (std::size_t i = 0; i < poly_vars.size(); ++i) {
      index_map_.emplace(poly_vars[i], i);
      // GCC 4.7 doesn't have emplace for std::map
      zero_valuation_.insert(std::make_pair(poly_vars[i], SR::null()));
    }
  }

  Matrix<SR> delta_at(const Matrix<SR> &newton_update,
                      const Matrix<SR> &previous_newton_values) {

    assert(previous_newton_values.getColumns() == 1 && newton_update.getColumns() == 1);

    auto num_variables = poly_vars.size();
    assert(num_variables == previous_newton_values.getRows() &&
           num_variables == newton_update.getRows());

    std::vector<SR> delta_vector;

    std::unordered_map<VarId, Degree> current_max_degree;

    std::unordered_map<VarId, SR> newton_update_map;

    for (std::size_t i = 0; i < num_variables; ++i) {
      current_valuation_[poly_vars[i]] = previous_newton_values.At(i, 0);
      newton_update_map.insert(
        std::make_pair(poly_vars[i], newton_update.At(i, 0)));
    }

    for (const auto &polynomial : polynomials) {
      SR delta = SR::null();
      Degree polynomial_max_degree = polynomial.get_degree();

      current_max_degree.clear();
      current_max_degree.insert(polynomial.GetVarDegreeMap().begin(),
                                polynomial.GetVarDegreeMap().end());

      if (polynomial_max_degree <= 1) {
        delta = SR::null();
      } else {
        delta =
          polynomial.AllNewtonDerivatives(current_valuation_, newton_update_map);
      }

      delta_vector.emplace_back(std::move(delta));
    }
    return Matrix<SR>(delta_vector.size(), std::move(delta_vector));
  }

private:
  const std::vector< Polynomial<SR> > &polynomials;
  const std::vector<VarId> &poly_vars;
  std::unordered_map<VarId, std::size_t> index_map_;
  /* We cache the std::map so that we can avoid reallocating it every time.  And
   * we also make sure that we always overwrite everything before using it... */
  std::unordered_map<VarId, SR> current_valuation_;
  std::unordered_map<VarId, SR> zero_valuation_;
};

/*
 * TODO: other signature for "solving" linsys ?? (with number of iterations?)
 * Be sure to call it only on semirings for which the iteration terminates...
 * TODO: implement round-robin iteration etc. (all the tricks from program-analysis.. see Tarjan's paper from '76)
 */
template <typename SR>
class SimpleKleeneLinSolver {

public:
  /*
   *
   */
  SimpleKleeneLinSolver(const std::vector< NonCommutativePolynomial<SR> >& F,  const std::vector<VarId>& variables) {
          polynomials_ = F;
  }

  /*
   * FIXME: very very ugly code.. lots of copying
   * compute f = "F.differential_at(values) + rhs" (=linear polynomial)
   * and then iterate f^n(0) until convergence or MAX_ITER is reached
   */
  Matrix<SR> solve_lin_at(const Matrix<SR>& values, const Matrix<SR>& rhs, const std::vector<VarId>& variables) {
          std::vector<NonCommutativePolynomial<SR> > F_lin;

          ValuationMap<SR> val_map = make_valmap(variables, values);


          // linearize F at point "values", obtain linear polynomials F_lin
          for(NonCommutativePolynomial<SR>& p : polynomials_) {
                  F_lin.push_back(p.differential_at(val_map));
          }

          int iteration=0;
          bool converged = false;

          ValuationMap<SR> iter_values = make_valmap(variables,rhs);
          ValuationMap<SR> iter_values_new;

          // with start s = rhs, compute s = F_lin(s) + rhs until convergence or iteration bound has been hit

          while(!converged && iteration < MAX_ITER) {
                  for(std::size_t j = 0; j < variables.size(); ++j) {
                          iter_values_new[variables[j]] = F_lin[j].eval(iter_values) + rhs.At(j,0);
                  }
                  if(iter_values == iter_values_new)
                          converged = true;
                  else
                          iter_values = iter_values_new;
                  ++iteration;
          }
          return Matrix<SR>(variables.size(),std::vector<SR>(iter_values.begin(), iter_values.end()));
  }

private:
  const int MAX_ITER = 10; // FIXME: for testing only---TODO: change interface of LinSolver ...
  std::vector< NonCommutativePolynomial<SR> >& polynomials_;

  ValuationMap<SR> make_valmap(const std::vector<VarId>& variables, const Matrix<SR>& rhs) {
    ValuationMap<SR> val_map;
    for (std::size_t i = 0; i < variables.size(); ++i) {
      // FIXME: GCC 4.7 is missing emplace
      val_map.insert(std::make_pair(variables[i], rhs.At(i, 0)));
    }
    return val_map;
  }
};


template <typename SR>
class NonCommutativeDeltaGenerator {
public:
        NonCommutativeDeltaGenerator(const std::vector< Polynomial<SR> >& polynomials, const std::vector<VarId> &poly_vars) {
    //FIXME: no need for copying... just keep a pointer to the polynomial system and the vars
    this->polynomials_ = polynomials;
    this->poly_vars = poly_vars;
  }

  virtual ~NonCommutativeDeltaGenerator(){ }

  Matrix<SR> delta_at(const Matrix<SR>& newton_update,
                                          const Matrix<SR>& previous_newton_values) {
        assert(previous_newton_values.getColumns() == 1 && newton_update.getColumns() == 1);

        std::vector<SR> delta;

          for(NonCommutativePolynomial<SR>& p : polynomials_) {
                  delta.push_back(p.calculate_delta(make_valmap(poly_vars,previous_newton_values),
                                                        make_valmap(poly_vars,newton_update)));

          }
          return Matrix<SR>(delta);
  }

private:
  std::vector< Polynomial<SR> > polynomials_;
  std::vector<VarId> poly_vars;

  ValuationMap<SR> make_valmap(const std::vector<VarId>& variables, const Matrix<SR>& rhs) {
    ValuationMap<SR> val_map;
    for (std::size_t i = 0; i < variables.size(); ++i) {
      // FIXME: GCC 4.7 is missing emplace
      val_map.insert(std::make_pair(variables[i], rhs.At(i, 0)));
    }
    return val_map;
  }
};



// compatability with old implementation
template <typename SR>
using Newton =
  GenericNewton<SR, CommutativeSymbolicLinSolver, CommutativeDeltaGenerator>;

template <typename SR>
using NewtonCL =
  GenericNewton<SR, CommutativeConcreteLinSolver, CommutativeDeltaGenerator>;

// default NonCommutative Newton implementation using naive Kleene iteration
template <typename SR>
using NonCommutativeNewton =
  GenericNewton<SR, SimpleKleeneLinSolver, NonCommutativeDeltaGenerator>;

#endif /* NEWTON_GENERIC_H_ */
