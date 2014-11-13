/*
 * semilinSetNdd.h
 *
 *  Created on: 22.01.2013
 *      Author: Michael Kerscher
 */

//#define USE_GENEPI

#ifdef USE_GENEPI
#ifndef SEMILINSETNDD_H_
#define SEMILINSETNDD_H_
#include <set>
#include <string>
#include <memory>
#include <genepi/genepi.h>
#include <genepi/genepi-loader.h>
#include "semiring.h"
#include "../datastructs/var.h"
#include "../datastructs/sparse_vec.h"
#include "../datastructs/vec_set.h"
#include "../utils/string_util.h"
#include "linear_set.h"
#include "semilinSetNdd_util.h"

class SemilinSetNdd : public StarableSemiring<SemilinSetNdd, Commutativity::Commutative, Idempotence::Idempotent>
{
public:
	static std::shared_ptr<SemilinSetNdd> elem_null;
	static std::shared_ptr<SemilinSetNdd> elem_one;

	static bool genepi_init();
	static bool genepi_dealloc();

	static bool solver_init(int number_of_variables);
  static bool solver_dealloc();


  static genepi_solver* solver;
  static genepi_plugin* plugin;
  static int k;
  static std::unordered_map<VarId, int> var_map;

	SemilinSetNdd();
	SemilinSetNdd(int zero);
  SemilinSetNdd(VarId var);
	SemilinSetNdd(VarId var, int cnt);
	SemilinSetNdd(const SemilinSetNdd& expr);
	virtual ~SemilinSetNdd();
  SemilinSetNdd operator =  (const SemilinSetNdd& term);
	SemilinSetNdd operator += (const SemilinSetNdd& term);
	SemilinSetNdd operator *= (const SemilinSetNdd& term);
	bool operator < (const SemilinSetNdd& term) const;
	bool operator == (const SemilinSetNdd& term) const;
	SemilinSetNdd star () const;
	static SemilinSetNdd null();
	static SemilinSetNdd one();
	std::string string() const;
	static const bool is_idempotent;
	static const bool is_commutative;

	//Generate an NDD from a semilinear set described explicitly as set of linear sets
		template <DIVIDER_TEMPLATE_TYPE Divider, VEC_SIMPL_TEMPLATE_TYPE VecSimpl>
	SemilinSetNdd(const VecSet<LinearSet<VarId, Counter, Divider, VecSimpl> >& linsets) {
		set = Genepi(solver, k, false);

		//std::cout << "building " << k <<"-NDD for: " << ToStringSorted(linsets, "\n ") << std::endl;

	  // for each linear set we create an NDD (resp. genepi-set) and compute the union over them
	  for(auto ls : linsets) {
	    SparseVec<VarId, Counter, Divider> offset = ls.GetOffset();
	    VecSet<SparseVec<VarId, Counter, Divider> > generators = ls.GetGenerators();

      std::vector<std::vector<int> > generator_set;
      for (const auto &g : generators) {
        generator_set.push_back(toIntVec(g));
      }

      std::vector<int> off = toIntVec(offset);

      this->set = this->set.union_op(Genepi(solver, Genepi::createLinSet(solver, off, generator_set)));
      this->offsets.push_back(off);
	  }
    //std::cout << "... done"  << std::endl;
	}


private:
  Genepi set;
  // offset[x][0] counts the number of occurences of the first variable (for order see var_map)
  // FIXME: offsets are not offsets anymore but the minimal base of offsets
  std::vector<std::vector<int>> offsets;
  SemilinSetNdd(Genepi set, std::vector<std::vector<int>> offsets);
  bool isGenerator(const std::vector<int>& offset, const std::vector<int>& candidate, int n) const;
  bool isGeneratorFor(const std::vector<int>& offset1, const std::vector<int>& candidate, const std::vector<int>& offset2) const;
  std::vector<std::vector<int>> getIndependentOffsets(const std::vector<std::vector<int>>& offsets) const;


  template <DIVIDER_TEMPLATE_TYPE Divider>
  std::vector<int> toIntVec(const SparseVec<VarId, Counter, Divider>& sv) {
      std::vector<std::pair<VarId, Counter>> var_val_vec = sv.getVector();
      std::vector<int> tmpvec = std::vector<int>(k,0);

      for(const std::pair<VarId, Counter> vv : var_val_vec) {
        auto v = var_map.find(vv.first);
        if(v == var_map.end())
        {
          var_map.insert(std::make_pair(vv.first, var_map.size()));
          v = var_map.find(vv.first);
        }

        int position = v->second;
        tmpvec.at(position) = vv.second;
      }

      return tmpvec;
  }

};

#endif /* SEMILINSETNDD_H_ */
#endif /* USE_GENEPI */
