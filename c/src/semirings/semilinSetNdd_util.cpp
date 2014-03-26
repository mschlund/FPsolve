#include <iomanip>
#include <iostream>
#include <sstream>
#include "semilinSetNdd_util.h"
#include "semilinSetNdd_getoffsets.h"

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

bool Genepi::operator<(const Genepi& g) const {
  return genepi_set_compare(this->solver, this->set, GENEPI_LT, g.set);
}

bool Genepi::operator==(const Genepi& g) const {
  return genepi_set_equal(this->solver, this->set, g.set);
}

int Genepi::getSize() const {
  return genepi_set_get_data_structure_size(this->solver, this->set);
}
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

std::vector<std::vector<int>> Genepi::getOffsets() const {
  // TODO: make this more beautiful...
  std::string filename = "read_offsets_from_automaton.dot";
  FILE* file = fopen(filename.c_str(), "w");
  genepi_set_display_data_structure(solver, set, file);
  fclose(file);
  int dummy = system("bash convert.sh read_offsets_from_automaton.dot > read_offsets_from_automaton.parseable");
  int k = genepi_set_get_width(this->solver, this->set);
  Parsed information = parse_automaton(k, "read_offsets_from_automaton.parseable");

  std::vector<std::vector<int>> offsets;
  for(auto f : information.finals) {
    std::vector<int> visited;
    visited.push_back(information.init);
    auto result = DepthFirst(&(information.graph), visited, f, k);
    offsets.insert(offsets.end(), result.begin(), result.end());
  }

  unlink("read_offsets_from_automaton.parseable");
  unlink(filename.c_str());
  return offsets;
}

genepi_set* Genepi::createVector(genepi_solver* solver, std::vector<int> x) {
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
// TODO: describe how alpha is structured (last component ist µ)
genepi_set* Genepi::createGenerator(genepi_solver* solver, std::vector<int> x) {
  genepi_set *aux1;
  genepi_set *aux2;
  std::vector<int> alpha(x.size()+1);
  genepi_set *result = genepi_set_top_N (solver, alpha.size());

  for (unsigned int i = 0; i < x.size(); i++)
  {
    alpha.at(i) = -1;
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
genepi_set* Genepi::createGenerators(genepi_solver* solver, std::vector<std::vector<int>> sets) {
  genepi_set *result = genepi_set_top_N(solver, sets.at(0).size());
  genepi_set *aux1;
  genepi_set *aux2;
  for(auto set : sets) {
    aux1 = createGenerator(solver, set);
    aux2 = genepi_set_intersection(solver, result, aux1);
    genepi_set_del_reference(solver, aux1);
    genepi_set_del_reference(solver, result);
    result = aux2;
  }
  return result;
}

