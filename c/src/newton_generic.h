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


template <typename SR, LIN_EQ_SOLVER_TYPE LinEqSolverTemplate, DELTA_GEN_TYPE DeltaGeneratorTemplate>
class GenericNewton {
  typedef LinEqSolverTemplate<SR> LinEqSolver;
  typedef DeltaGeneratorTemplate<SR> DeltaGenerator;
public:

  std::map<VarId,SR> solve_fixpoint(
      const std::vector<std::pair<VarId, Polynomial<SR>>>& equations,
      int max_iter) {
    std::vector<Polynomial<SR>> F;
    std::vector<VarId> poly_vars;
    for (auto equation_it = equations.begin(); equation_it != equations.end(); ++equation_it) {
      poly_vars.push_back(equation_it->first);
      F.push_back(equation_it->second);
    }
    Matrix<SR> result = this->solve_fixpoint(F, poly_vars, max_iter);

    // repack everything and return it
    std::map<VarId,SR> solution;
    auto result_vec = result.getElements();
    auto var_it = poly_vars.begin();
    for (auto result_it = result_vec.begin(); result_it != result_vec.end();
         ++result_it) {
      solution.insert(solution.begin(),
                      std::pair<VarId,SR>(*var_it, *result_it));
      ++var_it;
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
  Matrix<SR> solve_fixpoint(const std::vector< Polynomial<SR> >& polynomials,
                            const std::vector<VarId>& variables,
                            int max_iter) {

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
    DeltaGenerator DeltaGen = DeltaGenerator(F_mat, variables);

    for (int i=0; i<max_iter; ++i) {
      newton_update = LinSolver.solve_lin_at(newton_values,delta,variables);

      previous_newton_values = newton_values;

      if (!SR::IsIdempotent())
        newton_values = newton_values + newton_update;
      else
        newton_values = newton_update;

      if (!SR::IsIdempotent())
        delta = DeltaGen.delta_at(newton_update, previous_newton_values);
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
  //TODO: use initialization lists instead of temporaries+copying? see http://www.parashift.com/c++-faq-lite/init-lists.html
  CommutativeSymbolicLinSolver(const std::vector< Polynomial<SR> >& F,  const std::vector<VarId>& variables) {
    Matrix< Polynomial<SR> > jacobian = Polynomial<SR>::jacobian(F, variables);

    std::unordered_map<SR, VarId, SR> valuation_tmp;

    Matrix<FreeSemiring> jacobian_free = Polynomial<SR>::make_free(jacobian, &valuation_tmp);
    jacobian_star_ = new Matrix<FreeSemiring>(jacobian_free.star());

    std::cout << "J*: " << *jacobian_star_ << std::endl;

    for (auto &pair : valuation_tmp) {
      valuation_.insert(std::make_pair(pair.second, pair.first));
      std::cout << "valuation:  " << pair.second << ", " << pair.first << std::endl;
    }
  }

  virtual ~CommutativeSymbolicLinSolver(){
    delete jacobian_star_;
    jacobian_star_ = 0;
  }

  Matrix<SR> solve_lin_at(const Matrix<SR>& values, const Matrix<SR>& rhs, const std::vector<VarId>& variables) {
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
class CommutativeDeltaGenerator {
public:
  CommutativeDeltaGenerator(const Matrix< Polynomial<SR> >& F_mat, const std::vector<VarId> &poly_vars) {
    //FIXME: no need for copying... just keep a pointer to the polynomial system and the vars
    this->F_mat = new Matrix< Polynomial<SR> >(F_mat);
    this->poly_vars = poly_vars;
  }

  virtual ~CommutativeDeltaGenerator(){
    delete F_mat;
    F_mat =0;
  }

  /*
   * FIXME: get rid of vectors, always work on matrices -> saves conversion-time
   */

  Matrix<SR> delta_at(const Matrix<SR>& newton_update, const Matrix<SR>& previous_newton_values) {
    auto num_variables = newton_update.getRows();

    assert(num_variables == previous_newton_values.getRows() &&
           num_variables == poly_vars.size());


    std::vector<SR> delta;
    std::vector<Degree> current_max_degree(num_variables);


    for (std::size_t i = 0; i < num_variables; ++i) {
       SR delta_i = SR::null();
       Polynomial<SR> f = F_mat->At(i,0);
       Degree poly_max_degree = f.get_degree();

       for (std::size_t j = 0; j < num_variables; ++j) {
         current_max_degree[j] = f.GetMaxDegreeOf(poly_vars[j]);
       }

       /* We want to calculate all possible derivatives of at least second
        * order, but lower or equal to the degree of polynomial. */
       Generator generator{current_max_degree, 2, poly_max_degree};

       while (generator.NextCombination()) {
         std::map<VarId, Degree> dx;
         SR prod = SR::one();

         for (std::size_t index = 0; index < num_variables; ++index) {
           dx[poly_vars[index]] = generator.GetVectorRef()[index];
           for (int j = 0; j < generator.GetVectorRef()[index]; ++j) {
             prod *= newton_update.At(index,0);
           }
         }

         // eval f.derivative(dx) at v
         std::map<VarId, SR> value_map;
         for (std::size_t index = 0; index < num_variables; ++index) {
           // FIXME: GCC 4.7 is missing emplace
           // value_map.emplace(poly_vars[index], previous_newton_values[index]);
           value_map.insert(std::make_pair(poly_vars[index], previous_newton_values.At(index,0)));
         }

         SR f_eval = f.derivative_binom_at(dx,value_map);

         delta_i = delta_i + f_eval * prod;
       }

       delta.emplace_back(std::move(delta_i));
     }

     std::cout << "Delta: " << Matrix<SR >(delta.size(), delta)<<std::endl;

    return Matrix<SR>(delta.size(), delta);
  }

private:
  Matrix< Polynomial<SR> >* F_mat;
  std::vector<VarId> poly_vars;
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

  }

  /*
   * compute f = "F.differential_at(values) + rhs" (=linear polynomial)
   * and then iterate f^n(0) until convergence or MAX_ITER is reached
   */
  Matrix<SR> solve_lin_at(const Matrix<SR>& values, const Matrix<SR>& rhs, const std::vector<VarId>& variables) {

  }

private:

};




// compatability with old implementation
template <typename SR> class Newton : public GenericNewton<SR, CommutativeSymbolicLinSolver, CommutativeDeltaGenerator >{ };



#endif /* NEWTON_GENERIC_H_ */


