#ifndef NEWTON_H
#define NEWTON_H

#include <algorithm>
#include <cstdint>

#include "datastructs/matrix.h"
#include "datastructs/var_degree_map.h"
#include "polynomials/polynomial.h"

#include "matrix_free_semiring.h"



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

template <typename SR>
class Newton {
  private:
    Matrix<Polynomial<SR> > compute_symbolic_delta(
        const std::vector<VarId> &v,
        const std::vector<VarId> &v_upd,
        const std::vector<Polynomial<SR> > &F,
        const std::vector<VarId> &poly_vars) {

      auto num_variables = v.size();
      assert(num_variables == v_upd.size() &&
             num_variables == poly_vars.size());

      std::vector<Polynomial<SR> > delta;

      std::vector<Degree> current_max_degree(num_variables);

      for (std::size_t i = 0; i < num_variables; ++i) {
        Polynomial<SR> delta_i = Polynomial<SR>::null();
        Polynomial<SR> f = F.at(i);
        Degree poly_max_degree = f.get_degree();

        for (std::size_t j = 0; j < num_variables; ++j) {
          current_max_degree[j] = f.GetMaxDegreeOf(poly_vars[j]);
        }

        /* We want to calculate all possible derivatives of at least second
         * order, but lower or equal to the degree of polynomial. */
        Generator generator{current_max_degree, 2, poly_max_degree};

        while (generator.NextCombination()) {
          //std::vector<VarId> dx;
          VarDegreeMap dx;
          Polynomial<SR> prod{SR::one()};

          for (std::size_t index = 0; index < num_variables; ++index) {
            dx.Insert(poly_vars[index], generator.GetVectorRef()[index]);
            for (int j = 0; j < generator.GetVectorRef()[index]; ++j) {
              prod *= v_upd[index];
            }
          }

          // eval f.derivative(dx) at v
          std::map<VarId, VarId> values;
          for (std::size_t index = 0; index < v.size(); ++index) {
            // FIXME: GCC 4.7 is missing emplace
            // values.emplace(poly_vars[index], v[index]);
            values.insert(std::make_pair(poly_vars[index], v[index]));
          }
          Polynomial<SR> f_eval = f.derivative_binom(dx).subst(values);

          delta_i = delta_i + f_eval * prod;
        }

        delta.emplace_back(std::move(delta_i));
      }
      std::cout << "Delta: " << Matrix<Polynomial<SR> >(delta.size(), delta)<<std::endl;
      return Matrix<Polynomial<SR> >(delta.size(), delta);
    }

    // return a vector of fresh anonymous variables
    std::vector<VarId> get_symbolic_vector(int size) {
      std::vector<VarId> ret;
      for (int i=0; i<size; i++)
      {
        ret.push_back(Var::GetVarId());
      }
      return ret;
    }

  public:
    void UpdateValuation(const std::vector<VarId> &variables,
                         const Matrix<SR> &newton_values,
                         std::unordered_map<VarId, SR> *valuation) {
      assert(valuation != nullptr);
      assert(variables.size() == newton_values.getRows());
      assert(newton_values.getColumns() == 1);
      for (std::size_t i = 0; i < variables.size(); ++i) {
        (*valuation)[variables[i]] = newton_values.At(i, 0);
      }
    }

    /* Calculate the next newton iterand and update the valuation mapping. */
    Matrix<SR> step(const std::vector<VarId> &variables,
                    const Matrix<SR> &newton_values,
                    const Matrix<FreeSemiring> &jacobian_star,
                    const Matrix<SR> &delta,
                    std::unordered_map<VarId, SR> *valuation) {

      assert(valuation != nullptr);
      assert(variables.size() == newton_values.getRows());
      assert(newton_values.getColumns() == 1);

      UpdateValuation(variables, newton_values, valuation);

      return FreeSemiringMatrixEval(jacobian_star, *valuation) * delta;
    }

    // this is just a wrapper function at the moment
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
      // Matrix<SR> result = SolveNoSd(F, poly_vars, max_iter);

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

    // iterate until convergence
    // TODO: seems to be 2 iterations off compared to sage-impl..
    Matrix<SR> solve_fixpoint(const std::vector< Polynomial<SR> >& polynomials,
                              const std::vector<VarId>& variables,
                              int max_iter) {

      Matrix< Polynomial<SR> > F_mat{polynomials.size(), polynomials};
      Matrix< Polynomial<SR> > jacobian =
        Polynomial<SR>::jacobian(polynomials, variables);

      std::unordered_map<SR, VarId, SR> valuation_tmp;
      Matrix<FreeSemiring> jacobian_free =
        Polynomial<SR>::make_free(jacobian, &valuation_tmp);
      Matrix<FreeSemiring> jacobian_star = jacobian_free.star();

      //std::cout << "Jacobian: " << std::endl;
      //std::cout << jacobian << std::endl;

      //std::cout << "Jacobian (free): " << std::endl;
      //std::cout << jacobian_free << std::endl;

      std::unordered_map<VarId, SR> valuation;
      for (auto &pair : valuation_tmp) {
        valuation.insert(std::make_pair(pair.second, pair.first));
      }

      std::vector<VarId> u = this->get_symbolic_vector(variables.size());
      std::vector<VarId> u_upd =
        this->get_symbolic_vector(variables.size());

      // newton_values^0 = 0
      Matrix<SR> newton_values{polynomials.size(), 1};

      // d^0 = polynomials(0)
      std::map<VarId,SR> values;
      for (auto &variable : variables) {
        values.insert(std::make_pair(variable, SR::null()));
      }

      Matrix<SR> delta_new = Polynomial<SR>::eval(F_mat, values);

      Matrix<SR> newton_update = step(variables, newton_values, jacobian_star,
                                      delta_new, &valuation);

      // FIXME: ugly since we do not have a simple standard Matrix-Constructor
      // without any arguments
      Matrix<Polynomial<SR> > delta{1,1};

      if (!SR::IsIdempotent())
        delta = compute_symbolic_delta(u, u_upd, polynomials, variables);

      //start with i=1 as we have already done one iteration explicitly
      for (int i=1; i<max_iter; ++i) {
        if (!SR::IsIdempotent()) {
          values.clear();
          for (unsigned int i = 0; i<u.size(); i++) {
            values.insert(values.begin(),
                std::pair<VarId,SR>(u_upd.at(i), newton_update.getElements().at(i)));
            values.insert(values.begin(),
                std::pair<VarId,SR>(u.at(i), newton_values.getElements().at(i)));
          }
          delta_new = Polynomial<SR>::eval(delta,values);
        }
        std::cout << "Delta_new: " << delta_new << std::endl;


        if (SR::IsIdempotent())
          // for idempotent SRs we do not have do perform the addition (terms
          // are already accounted for in u_upd)!
          newton_values = newton_update;
        else
          newton_values = newton_values + newton_update;

        std::cout << "v: " << newton_values << std::endl;
        newton_update = step(variables, newton_values, jacobian_star,
                             delta_new, &valuation);
        std::cout << "v_upd: " << newton_update << std::endl;
      }

      if (SR::IsIdempotent())
        newton_values = newton_update;
      else
        newton_values = newton_values + newton_update;

      return newton_values;
    }


    Matrix<SR> ComputeDelta(
        const std::vector<VarId> &variables,
        const std::vector< Polynomial<SR> > &polynomials,
        const Matrix<SR> &newton_values,
        const Matrix<SR> &newton_update) {

      assert(newton_values.getColumns() == 1 && newton_update.getColumns() == 1);

      auto num_variables = variables.size();
      assert(num_variables == newton_values.getRows() &&
             num_variables == newton_update.getRows());

      std::vector<SR> delta_vector;

      std::vector<Degree> current_max_degree(num_variables);

      for (const auto &polynomial : polynomials) {
        SR delta = SR::null();
        Degree polynomial_max_degree = polynomial.get_degree();

        for (std::size_t i = 0; i < num_variables; ++i) {
          current_max_degree[i] = polynomial.GetMaxDegreeOf(variables[i]);
        }

        /* We want to calculate all possible derivatives of at least second
         * order, but lower or equal to the degree of polynomial. */
        Generator generator{current_max_degree, 2, polynomial_max_degree};

        while (generator.NextCombination()) {
          VarDegreeMap deriv_variables;

          SR prod = SR::one();
          for (std::size_t i = 0; i < num_variables; ++i) {
            assert(deriv_variables.GetDegreeOf(variables[i]) == 0);
            /* If we don't differentiate over the current variable, just
             * continue with the next one.  The map deriv_variables will
             * implicitly return 0 for the current variable. */
            if (generator.GetVectorRef()[i] == 0) {
              continue;
            }
            deriv_variables.Insert(variables[i], generator.GetVectorRef()[i]);
            for (int j = 0; j < generator.GetVectorRef()[i]; ++j) {
              prod *= newton_update.At(i, 0);
            }
          }

          std::map<VarId, SR> values;
          for (std::size_t i = 0; i < num_variables; ++i) {
            // FIXME: GCC 4.7 is missing emplace
            values.insert(std::make_pair(variables[i], newton_values.At(i, 0)));
          }
          SR polynomial_value =
            polynomial.DerivativeBinomAt(deriv_variables, values);

          delta += polynomial_value * prod;
        }

        delta_vector.emplace_back(std::move(delta));
      }
      return Matrix<SR>(delta_vector.size(), std::move(delta_vector));
    }



    /* Solve without using symbolic delta. */
    Matrix<SR> SolveNoSd(const std::vector< Polynomial<SR> > &polynomials,
                         const std::vector<VarId> &variables,
                         std::uint_fast64_t max_iter) {

      Matrix< Polynomial<SR> > polynomials_matrix{polynomials.size(), polynomials};
      Matrix< Polynomial<SR> > jacobian =
        Polynomial<SR>::jacobian(polynomials, variables);

      std::unordered_map<SR, VarId, SR> valuation_tmp;
      Matrix<FreeSemiring> jacobian_free =
        Polynomial<SR>::make_free(jacobian, &valuation_tmp);
      Matrix<FreeSemiring> jacobian_star = jacobian_free.star();

      std::unordered_map<VarId, SR> valuation;
      for (auto &pair : valuation_tmp) {
        valuation.insert(std::make_pair(pair.second, pair.first));
      }

      /* This is values^0 (i.e. \_ -> 0). */
      Matrix<SR> newton_values{polynomials.size(), 1};

      /* Temporary mapping just to get the first delta. */
      std::map<VarId, SR> tmp_values;
      for (auto &variable : variables) {
        // FIXME: GCC 4.7 doesn't have emplace
        tmp_values.insert(std::make_pair(variable, SR::null()));
      }

      /* This is delta^1. */
      Matrix<SR> delta = Polynomial<SR>::eval(polynomials_matrix, tmp_values);
      tmp_values.clear();

      /* This is update^1. */
      Matrix<SR> newton_update = step(variables, newton_values, jacobian_star,
                                      delta, &valuation);

      for (std::uint_fast16_t i = 1; i < max_iter; ++i) {

        /* Here we have newton_values^(i - 1), newton_update^i and delta^i. */

        /* This allows us to compute the delta^(i + 1) using:
         * - newton_update^i
         * - newton_values^{i - 1} */

        if (!SR::IsIdempotent()) {
          delta = ComputeDelta(variables, polynomials, newton_values, newton_update);
        }

        /* Then we compute the newton_values^i by (optionally) adding:
         * - newton_values^{i - 1}
         * - newton_update^i */

        if (SR::IsIdempotent()) {
          /* For idempotent SRs we do not have do perform the addition. */
          newton_values = std::move(newton_update);
        } else {
          newton_values = newton_values + newton_update;
        }

        /* Finally we compute the newton_update^{i + 1} using:
         * - newton_values^i
         * - delta^{i + 1} */

        newton_update = step(variables, newton_values, jacobian_star,
                             delta, &valuation);

        /* Here we have values^i, update^{i + 1} and delta^{i + 1}. */

      }

      if (SR::IsIdempotent()) {
        /* For idempotent SRs we do not have do perform the addition. */
        newton_values = std::move(newton_update);
      } else {
        newton_values = newton_values + newton_update;
      }

      return newton_values;
    }

};

#endif
