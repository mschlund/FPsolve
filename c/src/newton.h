#ifndef NEWTON_H
#define NEWTON_H

#include <cstdio>
#include "matrix.h"
#include "polynomial.h"

#ifndef OLD_FREESEMIRING
#include "free-semiring.h"
#else
#include "free-semiring-old.h"
#endif  /* OLD_FREESEMIRING */

template <typename SR>
class Newton {
  private:
    int sum(int* array, int n) {
      int s = 0;
      for (int i=0; i<n; i++) {
        s += array[i];
      }
      return s;
    };

    // generate vectors of vectors [[x0,x1,...xn,]] such that each possible permutation
    // of integers with 0 <= xi <max is in the result. also each permutation has a sum with
    // min_sum <= sum <= max_sum
    std::vector<std::vector<int> > genIdx(int max, int n, int min_sum,
                                          int max_sum) {
      std::vector<std::vector<int> > result;
      // initialize base array
      int* base = new int[n];

      while (base[n-1] <= max) {
        int k = 0;
        base[k]++;
        if (base[k]>max) {
          while (k<n && base[k]>max && base[n-1] <= max) {
            base[k] = 0;
            k++;
            base[k]++;
          }
        }
        int s = sum(base, n);
        if (min_sum <= s && s <= max_sum) {
          std::vector<int> tmp(base,base+n);
          result.push_back(tmp);
        }
      }

      delete[] base;
      return result;
    }

    Matrix<Polynomial<SR> > compute_symbolic_delta(
        const std::vector<VarPtr>& v,
        const std::vector<VarPtr>& v_upd,
        const std::vector<Polynomial<SR> >& F,
        const std::vector<VarPtr>& poly_vars) {

      std::vector<Polynomial<SR> > delta;
      int n = v.size();

      for (int i=0; i<n; ++i) {
        Polynomial<SR> delta_i = Polynomial<SR>::null();
        Polynomial<SR> f = F.at(i);
        int deg = f.get_degree();

        // create [[0,...,0],[0,...,deg],[deg,...,deg]]
        std::vector<std::vector<int> > p = genIdx(deg,n,2,deg);

        //iterate over (x,y,...,z) with x,y,...,z \in [0,deg+1]
        for (std::vector<std::vector<int> >::const_iterator it = p.begin();
             it != p.end(); ++it) {
          std::vector<VarPtr> dx;
          Polynomial<SR> prod = Polynomial<SR>(SR::one());

          std::vector<VarPtr>::const_iterator var = poly_vars.begin();
          std::vector<VarPtr>::const_iterator elem = v_upd.begin();
          for (std::vector<int>::const_iterator z = it->begin(); z != it->end(); ++z) {
            for (int j=0; j<(*z); ++j) {
              dx.push_back(*var); // generate a multiset of variables like "xxxyyzzzz";
              prod = prod * (*elem);	// prod = prod * elem^z ... (elem is a variable)
            }
            ++var; // next variable
            ++elem; // next element of vector v_upd
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
