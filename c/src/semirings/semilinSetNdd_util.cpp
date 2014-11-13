#include <iomanip>
#include <iostream>
#include <sstream>
#include "semilinSetNdd_util.h"

Genepi::Genepi() : solver(0), set(0) {
  // invalid state ...
}

Genepi::Genepi(genepi_solver* solver, genepi_set* set) : solver(solver), set(set) {
}

// constructor for empty or full set
Genepi::Genepi(genepi_solver* solver, int dimensions, bool complete) {
  this->solver = solver;
  if(complete)
    this->set = genepi_set_top_N(this->solver, dimensions);
  else
    this->set = genepi_set_bot(this->solver, dimensions);
}

// constructor for creating vectors or generators
Genepi::Genepi(genepi_solver* solver, std::vector<int> x, bool generator) {
  this->solver = solver;
  if(generator)
    this->set = createGenerator(this->solver, x);
  else
    this->set = createVector(this->solver, x);
}

// constructor for linear equality sets
Genepi::Genepi(genepi_solver* solver, std::vector<int> alpha, int constant) {
  this->solver = solver;
  this->set = genepi_set_linear_equality(this->solver, alpha.data(), alpha.size(), constant);
}

Genepi::Genepi(const Genepi& g) {
  this->solver = g.solver;
  this->set = g.set;
  genepi_set_add_reference(this->solver, this->set);
}

Genepi::~Genepi() {
  //if(this->solver != 0 && this->set != 0)
    genepi_set_del_reference(this->solver, this->set);
}

Genepi& Genepi::operator=(const Genepi& g) {
  this->solver = g.solver;
  // if(this->set != 0)
  //   genepi_set_del_reference(this->solver, this->set);
  this->set = g.set;
  genepi_set_add_reference(this->solver, this->set);
  return *this;
}

Genepi Genepi::intersect(const Genepi& g) const {
  return Genepi(this->solver, genepi_set_intersection(this->solver, this->set, g.set));
}

Genepi Genepi::union_op(const Genepi& g) const {
  return Genepi(this->solver, genepi_set_union(this->solver, this->set, g.set));
};

// TODO: ints needed? smaller type?
Genepi Genepi::project(const std::vector<int>& selection) const {
  return Genepi(this->solver, genepi_set_project(this->solver, this->set, selection.data(), selection.size()));
}

Genepi Genepi::invproject(const std::vector<int>& selection) const {
  return Genepi(this->solver, genepi_set_invproject(this->solver, this->set, selection.data(), selection.size()));
}

bool Genepi::isSolution(const std::vector<int>& solution) const {
  // TODO: figure out if 1 is correct for the xden parameter
  return genepi_set_is_solution(this->solver, this->set, solution.data(), solution.size(), 1);
}

bool Genepi::operator<(const Genepi& g) const {
  return genepi_set_compare(this->solver, this->set, GENEPI_LT, g.set);
}

bool Genepi::operator==(const Genepi& g) const {
  return genepi_set_equal(this->solver, this->set, g.set);
}

int Genepi::getSize() const {
  return genepi_set_get_data_structure_size(this->solver, this->set);
}

// FIXME: rewrite this. A lot of assumptions on path structure and available tools...
std::string Genepi::output(std::string prefix, int i, std::string postfix) const {
  std::stringstream ret;
  std::stringstream ss_filename;
  ss_filename << "automata/" << prefix << std::setfill('0') << std::setw(2) << std::hex << i << postfix << ".dot";
  // std::cout << ss_filename.str() << std::endl;
  FILE* file = fopen(ss_filename.str().c_str(), "w");
  genepi_set_display_data_structure(solver, set, file);
  fclose(file);
  std::stringstream ss_sed;
  ss_sed << "sed -i \"4d\" automata/" << prefix << std::setfill('0') << std::setw(2) << std::hex << i << postfix << ".dot";
  int ignore = system(ss_sed.str().c_str());
  std::stringstream ss;
  ss << "dot -Tsvg automata/" << prefix << std::setfill('0') << std::setw(2) << std::hex << i << postfix << ".dot > svgs/" << prefix << std::setfill('0') << std::setw(2) << std::hex << i << postfix << ".svg";
  ignore = system(ss.str().c_str());
  ret << "Automaton written to file : " << ss_filename.str();
  return ret.str();
}

genepi_set* Genepi::createVector(genepi_solver* solver, const std::vector<int>& x) {
  genepi_set *aux1;
  genepi_set *aux2;
  genepi_set *result = genepi_set_top_N (solver, x.size());
  std::vector<int> alpha(x.size());

  for (unsigned int i = 0; i < x.size(); i++)
  {
    alpha.at(i) = 1;
    aux1 = genepi_set_linear_equality (solver, alpha.data(), x.size(), x.at(i));
    aux2 = genepi_set_intersection (solver, aux1, result);
    genepi_set_del_reference (solver, result);
    genepi_set_del_reference (solver, aux1);
    result = aux2;
    alpha.at(i) = 0;
  }

  return result;
}
/*
 * For the input x = (x_1,...,x_n) we first compute
 * the set of all (y_1,...,y_n,µ) such that
 * y_i = µ * x_i (for all i) (equivalently (-1)*y_i + µ * x_i = 0)
 * then we project away the last component (µ)
 */
genepi_set* Genepi::createGenerator(genepi_solver* solver, const std::vector<int>& x) {
  genepi_set *aux1;
  genepi_set *aux2;
  std::vector<int> alpha(x.size()+1);
  genepi_set *result = genepi_set_top_N (solver, alpha.size());

  for (unsigned int i = 0; i < x.size(); i++)
  {
    alpha.at(i) = -1; // coefficient of y_i
    alpha.at(x.size()) = x.at(i);
    aux1 = genepi_set_linear_equality (solver, alpha.data(), alpha.size(), 0);
    aux2 = genepi_set_intersection (solver, aux1, result);
    genepi_set_del_reference (solver, result);
    genepi_set_del_reference (solver, aux1);
    result = aux2;
    alpha.at(i) = 0;
  }

  // we have one component (µ) too much. Project it away
  std::vector<int> selection(alpha.size(),0);
  selection.at(alpha.size()-1) = 1;
  aux1 = genepi_set_project(solver, result, selection.data(), selection.size());
  genepi_set_del_reference (solver, result);
  result = aux1;

  return result;
}

/*
 * For the input offset = (x_1,...,x_n) and generators = (g_1,...,g_k)
 * we first compute
 * the set of all (y_1,...,y_n,µ_1,...,µ_k) such that
 *  y_i = x_i + µ_1 * g_1,i + ... + µ_k * g_k,i
 * <=> (-1)*x_i = (-1)*y_i + µ_1 * g_1,i + ... + µ_k * g_k,i  (for all i)
 * then we project away the last k components (µ_1,...,µ_k)
 */
genepi_set* Genepi::createLinSet(genepi_solver* solver, const std::vector<int>& offset, const std::vector<std::vector<int>>& generators) {

  int n = offset.size();
  int k = generators.size();

  genepi_set *result = genepi_set_top_N(solver, n+k);
  genepi_set *aux1;
  genepi_set *aux2;

  for (unsigned int i = 0; i < n; i++)
  {
    std::vector<int> alpha(n + k, 0);
    alpha.at(i) = -1; // coefficient of y_i
    for (int j=0; j<k; j++) {
      alpha.at(n+j) = generators.at(j).at(i); // coefficient of g_j,i
    }
    aux1 = genepi_set_linear_equality (solver, alpha.data(), alpha.size(), -offset.at(i));
    aux2 = genepi_set_intersection (solver, aux1, result);
    genepi_set_del_reference (solver, result);
    genepi_set_del_reference (solver, aux1);
    result = aux2;
  }

  //project away all the µ's
  std::vector<int> selection(n+k,0); // an entry of "1" means we throw it AWAY in the projection ... yes it is counterintuitive
  for (int j=n; j<n+k; j++) {
    selection.at(j) = 1; // throw away µ_j
  }

  aux1 = genepi_set_project(solver, result, selection.data(), selection.size());
  genepi_set_del_reference (solver, result);
  result = aux1;

  return result;
}



