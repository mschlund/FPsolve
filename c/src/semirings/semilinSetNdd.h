/*
 * semilinSetNdd.h
 *
 *  Created on: 22.01.2013
 *      Author: Michael Kerscher
 */

#ifndef SEMILINSETNDD_H_
#define SEMILINSETNDD_H_
#include <set>
#include <string>
#include <memory>
#include "../datastructs/var.h"
#include "semiring.h"
#include <genepi/genepi.h>
#include <genepi/genepi-loader.h>
#include "semilinSetNdd_util.h"

class SemilinSetNdd : public StarableSemiring<SemilinSetNdd, Commutativity::Commutative, Idempotence::Idempotent>
{
public:
private:
        Genepi set;
        // offset[x][0] counts the number of occurences of the first variable (for order see var_map)
        // FIXME: offsets are not offsets anymore but the minimal base of offsets
        std::vector<std::vector<int>> offsets;
        SemilinSetNdd(Genepi set, std::vector<std::vector<int>> offsets);
        bool isGenerator(const std::vector<int>& offset, const std::vector<int>& candidate, int n) const;
        bool isGeneratorFor(const std::vector<int>& offset1, const std::vector<int>& candidate, const std::vector<int>& offset2) const;
        std::vector<std::vector<int>> getIndependentOffsets(const std::vector<std::vector<int>>& offsets) const;

public:
	static std::shared_ptr<SemilinSetNdd> elem_null;
	static std::shared_ptr<SemilinSetNdd> elem_one;

        static bool genepi_init(std::string plugin, int number_of_variables);
        static bool genepi_dealloc();
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
};

#endif /* SEMILINSETNDD_H_ */
