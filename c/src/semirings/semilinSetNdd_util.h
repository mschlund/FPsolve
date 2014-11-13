#ifndef SEMILINSETNDD_UTIL_H_
#define SEMILINSETNDD_UTIL_H_

#include <genepi/genepi.h>
#include <vector>

#include "../utils/string_util.h"

class Genepi {
  // TODO: less flags and more enums
private:
  genepi_set *set = nullptr;
  genepi_solver* solver = nullptr; // TODO: maybe this does not have to be in here...
public:
  Genepi(); // should not be called
  Genepi(genepi_solver* solver, genepi_set* set);
  // constructor for empty or full set
  Genepi(genepi_solver* solver, int dimensions, bool complete);
  // constructor for creating vectors or generators
  Genepi(genepi_solver* solver, std::vector<int> x, bool generator);
  // constructor for linear equality sets
  Genepi(genepi_solver* solver, std::vector<int> alpha, int constant);
  Genepi(const Genepi& g);
  ~Genepi();

  Genepi& operator=(const Genepi& g);
  Genepi intersect(const Genepi& g) const;
  Genepi union_op(const Genepi& g) const;
  // TODO: ints needed? smaller type?
  Genepi project(const std::vector<int>& selection) const;
  Genepi invproject(const std::vector<int>& selection) const;
  bool isSolution(const std::vector<int>& solution) const;
  std::vector<std::vector<int>> getOffsets() const;
  bool operator<(const Genepi& g) const;
  bool operator==(const Genepi& g) const;
  int getSize() const;
  std::string output(std::string prefix, int i, std::string postfix) const;

  static genepi_set* createVector(genepi_solver* solver, const std::vector<int>& x);
  static genepi_set* createGenerator(genepi_solver* solver, const std::vector<int>& x);
  static genepi_set* createLinSet(genepi_solver* solver, const std::vector<int>& offset, const std::vector<std::vector<int>>& generators);
};

#endif
