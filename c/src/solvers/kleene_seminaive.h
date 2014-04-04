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
class Kleene {
  typedef std::vector< std::pair< VarId, Poly<SR> > > GenericEquations;

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
  ValuationMap<SR> solve_fixpoint(const GenericEquations& equations, int max_iter) {

    std::vector<Poly<SR>> F;
    std::vector<VarId> poly_vars;

    for (const auto &eq : equations) {
      poly_vars.push_back(eq.first);
      F.push_back(eq.second);
    }

    std::vector<Poly<SR> > unfolded_polys;

    // we use the original poly_vars for the "X^{=h+1}"-variables
    std::unordered_map<VarId, VarId> prev_val_map; // X^{<h}
    std::unordered_map<VarId, VarId> val_map;      // X^{<h+1}


    for (Poly<SR> f : F) {
      //the new anonymous unfolding-variables are accumulated in the two maps!
      unfolded_polys.push_back(f.HeightUnfolding(prev_val_map, val_map));
    }

    //TODO: use vectors for better cache efficiency?
    ValuationMap<SR> prev_values; // valuations of X^{<h}
    ValuationMap<SR> values; // valuations of X^{<h+1}
    ValuationMap<SR> updates; // valuations of X^{=h+1}

    ValuationMap<SR> zeros;
    for (unsigned int i=0; i<poly_vars.size(); ++i) {
      zeros.insert({poly_vars[i], SR::null()});
    }

    // Init all values, note that they all talk about different variable names (connected via the two maps above)
    for (unsigned int i=0; i<poly_vars.size(); ++i) {
      SR f_0 = F[i].eval(zeros);
      prev_values.insert({prev_val_map.at(poly_vars[i]), SR::null()});
      values.insert({val_map.at(poly_vars[i]), f_0});
      updates.insert(poly_vars[i], f_0);
    }


    for (unsigned int i=0; i < max_iter; ++i) {

    }


  }
};



#endif /* KLEENE_SEMINAIVE_H_ */
