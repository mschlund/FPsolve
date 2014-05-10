/*
 * kleene_seminaive.h
 *
 *  Created on: 31.03.2014
 *      Author: schlund
 */

#ifndef KLEENE_SEMINAIVE_H_
#define KLEENE_SEMINAIVE_H_


#include <vector>
#include "../datastructs/var.h"

// Polynomial is parametrized by a semiring
#define POLY_TYPE template <typename> class

template <typename SR, POLY_TYPE Poly>
class KleeneSeminaive {

public:
  /*
   * Init:
   * previous_values = 0
   * values = F(0)
   * update = F(0)
   *
   * precompute compute polynomials.HeightUnfolding()
   *
   * Each iteration consists of :
   *
   * update = evaluate the height unfolding at respective values
   * previous_values = values
   * values = values + update
  */
  ValuationMap<SR> solve_fixpoint(const GenericEquations<Poly, SR>& equations, int max_iter) {

    std::vector<Poly<SR>> F;
    std::vector<VarId> poly_vars;

    for (const auto &eq : equations) {
      poly_vars.push_back(eq.first);
      F.push_back(eq.second);
    }

    std::vector<Poly<SR> > unfolded_polys;

    // we use the original poly_vars for the "X^{=h+1}"-variables
    SubstitutionMap prev_val_map; // X^{<h}
    SubstitutionMap val_map;      // X^{<h+1}

    //std::cout << "Unfolding polys" << std::endl;
    for (Poly<SR> f : F) {
      //the new anonymous unfolding-variables are accumulated in the two maps!
      Poly<SR> t = f.HeightUnfolding(prev_val_map, val_map);
      unfolded_polys.push_back(t);
      //std::cout << t << std::endl;
    }

    //std::cout << "Unfolding done!" << std::endl;
    //for (int i=0; i<poly_vars.size(); i++)
    //  std::cout << poly_vars[i] << std::endl;

    //std::cout << prev_val_map<< std::endl;
    //std::cout << val_map<< std::endl;

    //TODO: use vectors for better cache efficiency?

    // valuations of X^{<h},  X^{<h+1}, and X (representing X^{=h+1})
    // keep all in one map to pass it easily to eval
    ValuationMap<SR> all_values;
    ValuationMap<SR> updates;


    ValuationMap<SR> zeros;
    for (unsigned int i=0; i<poly_vars.size(); ++i) {
      zeros.insert({poly_vars[i], SR::null()});
    }

    // Init all values, note that they all talk about different variable names (connected via the two maps above)

    // FIXME: if a variable DOES NOT appear on the rhs (e.g. start symbol S for a PCFG) this will give an exception!!!


    for (unsigned int i=0; i<poly_vars.size(); ++i) {
      SR f_0 = F[i].eval(zeros);
      all_values.insert({prev_val_map.at(poly_vars[i]), SR::null()});
      all_values.insert({val_map.at(poly_vars[i]), f_0});
      all_values.insert({poly_vars[i], f_0});
    }


    for (unsigned int n=0; n < max_iter; ++n) {

      for (unsigned int i=0; i<poly_vars.size(); ++i) {
      // compute new update, note that we cannot modify all_values, yet!
        updates[poly_vars[i]] = unfolded_polys[i].eval(all_values);
      }

      // now all effects have been computed -> update the valuation
      for (unsigned int i=0; i<poly_vars.size(); ++i) {
        all_values.at(poly_vars[i]) = updates.at(poly_vars[i]);
        all_values.at(prev_val_map.at(poly_vars[i])) = all_values.at(val_map.at(poly_vars[i])); //safe vals
        all_values.at(val_map.at(poly_vars[i])) += all_values.at(poly_vars[i]); // add update to vals
      }
    }

    ValuationMap<SR> result;
    for (unsigned int i=0; i<poly_vars.size(); ++i) {
      result.insert({poly_vars[i], all_values.at(val_map.at(poly_vars[i]))});
    }

    return result;
  }
};

// compatability with old implementation
template <typename SR>
using KleeneComm =
    KleeneSeminaive<SR, CommutativePolynomial>;




#endif /* KLEENE_SEMINAIVE_H_ */
