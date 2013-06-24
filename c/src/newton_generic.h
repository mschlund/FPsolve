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


/* This defines the generator that is able to create all possible combinations
 * of integers such they are smaller than max and their sum is between min_sum
 * and max_sum. */
class Generator {
  public:
    Generator(const std::vector<Degree> &max, Degree min_sum, Degree max_sum)
        : vector_(max.size()), max_(max), current_sum_(0),
          min_sum_(min_sum), max_sum_(max_sum) {
      assert(0 < vector_.size());
      assert(min_sum <= max_sum);
      assert(current_sum_ == CurrentSum());
      assert(min_sum <= std::accumulate(max_.begin(), max_.end(),
                                        static_cast<Degree>(0)));
    }

    bool NextCombination() {
      bool added = false;
      bool valid = false;

      do {
        if (current_sum_ < min_sum_) {
          added = JumpMin();
        } else if (max_sum_ < current_sum_) {
          added = JumpMax();
        } else {
          added = AddOne();
        }
        assert(current_sum_ == CurrentSum() && BelowEqualMax());
        valid = min_sum_ <= current_sum_ && current_sum_ <= max_sum_;
        if (added && valid) {
          return true;
        }
      } while (added && !valid);

      return false;
    }

    const std::vector<std::uint16_t>& GetVectorRef() const {
      assert(current_sum_ == CurrentSum());
      assert(BelowEqualMax());
      assert(min_sum_ <= current_sum_ && current_sum_ <= max_sum_);

      return vector_;
    }

  private:
    std::uint16_t CurrentSum() const {
      return std::accumulate(vector_.begin(), vector_.end(), 0);
    }

    bool BelowEqualMax() const {
      for (std::size_t i = 0; i < vector_.size(); ++i) {
        if (vector_[i] > max_[i]) {
          return false;
        }
      }
      return true;
    }


    /* Add 1 to the current vector, wrap-around if some value is > max_.
     * Returns false if we cannot add 1 (i.e., the last element would
     * overflow). */
    bool AddOne(std::size_t start_index = 0) {
      for (auto i = start_index; i < vector_.size(); ++i) {
        if (max_[i] == 0) {
          continue;
        }

        auto &integer = vector_[i];
        ++integer;
        ++current_sum_;
        if (integer <= max_[i]) {
          return true;
        }
        current_sum_ -= integer;
        integer = 0;
      }
      return false;
    }

    /* Create the smallest vector that satisfies the min_sum_ requirement.  This
     * means adding (without overflowing) min_sum_ - current_sum_. */
    bool JumpMin() {
      assert(min_sum_ > current_sum_);
      auto remaining = min_sum_ - current_sum_;
      for (std::size_t i = 0; i < vector_.size(); ++i) { // auto &integer : vector_) {
        if (max_[i] == 0) {
          continue;
        }
        auto &integer = vector_[i];
        auto integer_max = max_[i];
        auto to_add = integer + remaining > integer_max ?
                      integer_max - integer : remaining;
        integer += to_add;
        current_sum_ += to_add;
        remaining -= to_add;
        if (remaining == 0) {
          return true;
        }
      }
      /* Added as much as we could, but still not enough... */
      return false;
    }

    /* Create the smallest vector that satisfies the max_sum_ requirement.  This
     * means that we add (with overflow) enough to get below max. */
    bool JumpMax() {
      assert(max_sum_ < current_sum_);
      for (std::size_t i = 0; i < vector_.size(); ++i) {
        if (max_[i] == 0) {
          continue;
        }

        auto &integer = vector_[i];

        /* If integer == 0 then this is harmless. */
        current_sum_ -= integer;
        integer = 0;

        /* We're wrapping around integer and should add 1 to the next position.
         * Check if that is enough or whether we should try to wrap-around the
         * next position too.  This can happen when integer == 1. */
        if (current_sum_ - integer + 1 <= max_sum_) {
          return AddOne(i + 1);
        }
      }
      return false;
    }

    std::vector<std::uint16_t> vector_;
    const std::vector<Degree> &max_;
    Degree current_sum_;
    Degree min_sum_;
    Degree max_sum_;
};


template <typename SR,
          LIN_EQ_SOLVER_TYPE LinEqSolverTemplate,
          DELTA_GEN_TYPE DeltaGeneratorTemplate>
class GenericNewton {
  typedef LinEqSolverTemplate<SR> LinEqSolver;
  typedef DeltaGeneratorTemplate<SR> DeltaGenerator;

  public:
  std::map<VarId,SR> solve_fixpoint(const Equations<SR>& equations, int max_iter) {
    std::vector<Polynomial<SR>> F;
    std::vector<VarId> poly_vars;
    for (const auto &eq : equations) {
      poly_vars.push_back(eq.first);
      F.push_back(eq.second);
    }
    Matrix<SR> result = solve_fixpoint(F, poly_vars, max_iter);

    // repack everything and return it
    std::map<VarId,SR> solution;
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

    std::map<VarId,SR> values;
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
    return Matrix<SR>{jacobian_.getRows(), std::move(result_vec)}.FloydWarshall()
           * rhs;
  }

  private:
    Matrix< Polynomial<SR> > jacobian_;
    /* We don't want to allocate the map every time, especially since the keys
     * do not change... */
    std::map<VarId, SR> valuation_;
};





template <typename SR>
class CommutativeDeltaGenerator {
public:
  CommutativeDeltaGenerator(const std::vector< Polynomial<SR> >& polynomials,
                            const std::vector<VarId> &poly_vars) {
    //FIXME: no need for copying... just keep a pointer to the polynomial system and the vars
    this->polynomials = polynomials;
    this->poly_vars = poly_vars;
  }

  virtual ~CommutativeDeltaGenerator(){ }

  Matrix<SR> delta_at(const Matrix<SR> &newton_update,
                      const Matrix<SR> &previous_newton_values) {

    assert(previous_newton_values.getColumns() == 1 && newton_update.getColumns() == 1);

    auto num_variables = poly_vars.size();
    assert(num_variables == previous_newton_values.getRows() &&
           num_variables == newton_update.getRows());

    std::vector<SR> delta_vector;

    std::vector<Degree> current_max_degree(num_variables);

    for (const auto &polynomial : polynomials) {
      SR delta = SR::null();
      Degree polynomial_max_degree = polynomial.get_degree();

      for (std::size_t i = 0; i < num_variables; ++i) {
              current_max_degree[i] = polynomial.GetMaxDegreeOf(poly_vars[i]);
      }


      if (polynomial_max_degree <= 1) {
        for (auto &variable : poly_vars) {
          tmp_valuation_[variable] = SR::null();
        }
        delta = polynomial.eval(tmp_valuation_);

      } else {
        /* We want to calculate all possible derivatives of at least second
         * order, but lower or equal to the degree of polynomial. */
        Generator generator{current_max_degree, 2, polynomial_max_degree};

        while (generator.NextCombination()) {
          VarDegreeMap deriv_variables;

          SR prod = SR::one();
          for (std::size_t i = 0; i < num_variables; ++i) {
                  assert(deriv_variables.GetDegreeOf(poly_vars[i]) == 0);
                  /* If we don't differentiate over the current variable, just
                   * continue with the next one. The map deriv_variables will
                   * implicitly return 0 for the current variable. */
                  if (generator.GetVectorRef()[i] == 0) {
                          continue;
                  }
                  deriv_variables.Insert(poly_vars[i], generator.GetVectorRef()[i]);
                  for (int j = 0; j < generator.GetVectorRef()[i]; ++j) {
                          prod *= newton_update.At(i, 0);
                  }
          }

          for (std::size_t i = 0; i < num_variables; ++i) {
            // FIXME: GCC 4.7 is missing emplace
            tmp_valuation_[poly_vars[i]] = previous_newton_values.At(i, 0);
          }
          SR polynomial_value =
            polynomial.DerivativeBinomAt(deriv_variables, tmp_valuation_);

          delta += polynomial_value * prod;
        }
      }

      delta_vector.emplace_back(std::move(delta));
    }
    return Matrix<SR>(delta_vector.size(), std::move(delta_vector));
  }

private:
  std::vector< Polynomial<SR> > polynomials;
  std::vector<VarId> poly_vars;
  /* We cache the std::map so that we can avoid reallocating it every time.  And
   * we also make sure that we always overwrite everything before using it... */
  std::map<VarId, SR> tmp_valuation_;
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

          std::map<VarId, SR> val_map = make_valmap(variables, values);


          // linearize F at point "values", obtain linear polynomials F_lin
          for(NonCommutativePolynomial<SR>& p : polynomials_) {
                  F_lin.push_back(p.differential_at(val_map));
          }

          int iteration=0;
          bool converged = false;

          std::map<VarId, SR> iter_values = make_valmap(variables,rhs);
          std::map<VarId, SR> iter_values_new;

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

  std::map<VarId, SR> make_valmap(const std::vector<VarId>& variables, const Matrix<SR>& rhs) {
          std::map<VarId, SR> val_map;
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

  std::map<VarId, SR> make_valmap(const std::vector<VarId>& variables, const Matrix<SR>& rhs) {
          std::map<VarId, SR> val_map;
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
