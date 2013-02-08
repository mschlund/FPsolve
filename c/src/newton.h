#ifndef NEWTON_H
#define NEWTON_H

#include <cstdint>
#include <algorithm>

#include "matrix.h"
#include "polynomial.h"

#ifndef OLD_FREESEMIRING
#include "free-semiring.h"
#else
#include "free-semiring-old.h"
#endif  /* OLD_FREESEMIRING */


/* This defines the generator that is able to create all possible combinations
 * of integers such they are smaller than max and their sum is between min_sum
 * and max_sum. */
class Generator {
  public:
    Generator(std::vector<std::uint32_t>::size_type size, std::uint32_t max,
              std::uint32_t min_sum, std::uint32_t max_sum)
        : vector_(size), max_(max), min_sum_(min_sum), max_sum_(max_sum) {
      assert(0 < vector_.size());
      assert(0 < max);
      assert(min_sum <= max_sum);
      assert(min_sum <= vector_.size() * max);
    }

    bool NextCombination() {
      std::uint32_t sum = 0;
      bool added = false;
      bool valid = false;

      do {
        added = AddOne();
        sum = CurrentSum();
        valid = min_sum_ <= sum && sum <= max_sum_;
        if (added && valid) {
          return true;
        }
      } while (added && !valid);

      return false;
    }

    const std::vector<std::uint32_t>& GetVectorRef() const { return vector_; }

  private:
    std::uint32_t CurrentSum() const {
      assert(std::all_of(vector_.begin(), vector_.end(),
                         [this](std::uint32_t i) { return i <= max_; }));
      return std::accumulate(vector_.begin(), vector_.end(), 0);
    }

    /* Add 1 to the current vector, wrap-around if some value is > max.  Returns
     * false if we cannot add 1 (i.e., the last element would overflow). */
    bool AddOne() {
      for (auto &integer : vector_) {
        ++integer;
        if (integer <= max_) {
          return true;
        }
        integer = 0;
      }
      return false;
    }


    std::vector<std::uint32_t> vector_;
    std::uint32_t max_;
    std::uint32_t min_sum_;
    std::uint32_t max_sum_;
};

template <typename SR>
class Newton {
  private:
    Matrix<Polynomial<SR> > compute_symbolic_delta(
        const std::vector<VarPtr>& v,
        const std::vector<VarPtr>& v_upd,
        const std::vector<Polynomial<SR> >& F,
        const std::vector<VarPtr>& poly_vars) {

      auto num_variables = v.size();
      assert(num_variables == v_upd.size() &&
             num_variables == poly_vars.size());

      std::vector<Polynomial<SR> > delta;

      for (int i = 0; i < num_variables; ++i) {
        Polynomial<SR> delta_i = Polynomial<SR>::null();
        Polynomial<SR> f = F.at(i);
        int degree = f.get_degree();

        /* We want to calculate all possible derivatives of at least second
         * order, but lower or equal to the degree of polynomial. */
        Generator generator{num_variables, degree, 2, degree};

        while (generator.NextCombination()) {
          std::vector<VarPtr> dx;
          Polynomial<SR> prod{SR::one()};

          for (auto index = 0; index < num_variables; ++index) {
            for (int j = 0; j < generator.GetVectorRef()[index]; ++j) {
              dx.push_back(poly_vars[index]);
              prod *= v_upd[index];
            }
          }

          // eval f.derivative(dx) at v
          std::map<VarPtr,VarPtr> values;
          int j = 0;
          for (std::vector<VarPtr>::const_iterator poly_var = poly_vars.begin();
               poly_var != poly_vars.end(); ++poly_var) {
            values.insert(values.begin(),
                          std::pair<VarPtr,VarPtr>((*poly_var),v.at(j++)));
          }
          Polynomial<SR> f_eval = f.derivative(dx).subst(values);

          delta_i = delta_i + f_eval * prod;
        }

        delta.push_back(delta_i);
      }

      return Matrix<Polynomial<SR> >(1, delta.size(), delta);
    }

    std::vector<VarPtr> get_symbolic_vector(int size, std::string prefix) {
      // define new symbolic vector [u1,u2,...,un] TODO: this is ugly...
      std::vector<VarPtr> ret;
      for (int i=0; i<size; i++)
      {
        std::stringstream ss;
        ss << prefix << "_" << i;
        ret.push_back(Var::getVar(ss.str()));
      }
      return ret;
    }

  public:
    // calculate the next newton iterand
    Matrix<SR> step(const std::vector<VarPtr>& poly_vars,
        const Matrix<FreeSemiring>& J_s,
        std::unordered_map<VarPtr, SR>* valuation,
        const Matrix<SR>& v, const Matrix<SR>& delta) {
      assert(poly_vars.size() == (unsigned int)v.getRows());
      int i=0;
      for (std::vector<VarPtr>::const_iterator poly_var = poly_vars.begin();
          poly_var != poly_vars.end(); ++poly_var) {
        SR sr_elem = v.getElements().at(i++);
        valuation->erase(*poly_var); // clean the old variables from the map
        valuation->insert(valuation->begin(),
            std::pair<VarPtr,SR>(*poly_var,sr_elem));
      }
      Matrix<SR> J_s_new = FreeSemiring_eval<SR>(J_s, valuation);

      //std::cout << "Jacobian (evaluated): " << std::endl;
      //std::cout << J_s_new << std::endl;

      Matrix<SR> result = J_s_new * delta;
      return result;
    }

    // this is just a wrapper function at the moment
    std::map<VarPtr,SR> solve_fixpoint(
        const std::vector<std::pair<VarPtr, Polynomial<SR>>>& equations,
        int max_iter) {
      std::vector<Polynomial<SR>> F;
      std::vector<VarPtr> poly_vars;
      for (auto equation_it = equations.begin(); equation_it != equations.end(); ++equation_it) {
        poly_vars.push_back(equation_it->first);
        F.push_back(equation_it->second);
      }
      Matrix<SR> result = this->solve_fixpoint(F, poly_vars, max_iter);

      // repack everything and return it
      std::map<VarPtr,SR> solution;
      auto result_vec = result.getElements();
      auto var_it = poly_vars.begin();
      for (auto result_it = result_vec.begin(); result_it != result_vec.end();
           ++result_it) {
        solution.insert(solution.begin(),
                        std::pair<VarPtr,SR>(*var_it, *result_it));
        ++var_it;
      }
      return solution;
    }

    // iterate until convergence
    // TODO: seems to be 2 iterations off compared to sage-impl..
    Matrix<SR> solve_fixpoint(const std::vector<Polynomial<SR> >& F,
                              const std::vector<VarPtr>& poly_vars, int max_iter) {
      Matrix<Polynomial<SR> > F_mat = Matrix<Polynomial<SR> >(1,F.size(),F);
      Matrix<Polynomial<SR> > J = Polynomial<SR>::jacobian(F, poly_vars);
      auto valuation_tmp = new std::unordered_map<SR, VarPtr, SR>();
      Matrix<FreeSemiring> J_free = Polynomial<SR>::make_free(J, valuation_tmp);

      auto valuation = new std::unordered_map<VarPtr, SR>();
      // insert null and one valuations into the map
      // valuation->insert(valuation->begin(), std::pair<FreeSemiring,SR>(FreeSemiring::null(), SR::null()));
      // valuation->insert(valuation->begin(), std::pair<FreeSemiring,SR>(FreeSemiring::one(), SR::one()));
      for (auto v_it = valuation_tmp->begin(); v_it != valuation_tmp->end();
           ++v_it) {
        valuation->insert(valuation->begin(),
                          std::pair<VarPtr, SR>(v_it->second, v_it->first));
      }

      // std::cout << "Jacobian (with vars): " << std::endl;
      // std::cout << J << std::endl;

      Matrix<FreeSemiring> J_s = J_free.star();

      // define new symbolic vectors [u1,u2,...,un] TODO: this is ugly...
      std::vector<VarPtr> u = this->get_symbolic_vector(poly_vars.size(), "u");
      std::vector<VarPtr> u_upd =
        this->get_symbolic_vector(poly_vars.size(), "u_upd");

      Matrix<SR> v = Matrix<SR>(1,(int)F.size()); // v^0 = 0

      // d^0 = F(0)
      std::map<VarPtr,SR> values;
      for (std::vector<VarPtr>::const_iterator poly_var = poly_vars.begin();
           poly_var != poly_vars.end(); ++poly_var) {
        values.insert(values.begin(), std::pair<VarPtr,SR>(*poly_var, SR::null()));
      }
      Matrix<SR> delta_new = Polynomial<SR>::eval(F_mat, values);

      Matrix<SR> v_upd = step(poly_vars, J_s, valuation, v, delta_new);

      // FIXME: ugly since we do not have a simple standard Matrix-Constructor
      // without any arguments
      Matrix<Polynomial<SR> > delta = Matrix<Polynomial<SR> >(1,1);

      if (!SR::is_idempotent)
        delta = compute_symbolic_delta(u, u_upd, F, poly_vars);

      //start with i=1 as we have already done one iteration explicitly
      for (int i=1; i<max_iter; ++i) {
        if (!SR::is_idempotent) {
          values.clear();
          for (unsigned int i = 0; i<u.size(); i++) {
            values.insert(values.begin(),
                std::pair<VarPtr,SR>(u_upd.at(i), v_upd.getElements().at(i)));
            values.insert(values.begin(),
                std::pair<VarPtr,SR>(u.at(i), v.getElements().at(i)));
          }
          delta_new = Polynomial<SR>::eval(delta,values);
        }

        if (SR::is_idempotent)
          // for idempotent SRs we do not have do perform the addition (terms
          // are already accounted for in u_upd)!
          v = v_upd;
        else
          v = v + v_upd;

        v_upd = step(poly_vars, J_s, valuation, v, delta_new);
      }

      if (SR::is_idempotent)
        v = v_upd;
      else
        v = v + v_upd;

      delete valuation_tmp;
      delete valuation;

      return v;
    }
};

#endif
