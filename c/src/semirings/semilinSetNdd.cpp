#include <iostream>
#include "semilinSetNdd.h"

int static_i = 0;

// insert an offset if it is unique in offsets
std::vector<std::vector<int>>& insert_offset(std::vector<std::vector<int>>& offsets, const std::vector<int>& offset) {
  if(std::find(offsets.begin(), offsets.end(), offset) == offsets.end())
  {
    offsets.push_back(offset);
  }
  return offsets;
}

std::vector<std::vector<int>> multiply_offsets(const std::vector<std::vector<int>>& offsets1, const std::vector<std::vector<int>>& offsets2) {
  std::vector<std::vector<int>> offsets;
  for(auto o1 : offsets1) {
    for(auto o2 : offsets2) {
      assert(o1.size() == o2.size());
      std::vector<int> o(o1.size());
      for(int i=0; i<o1.size(); ++i) {
        o.at(i) = o1.at(i) + o2.at(i);
      }
      offsets = insert_offset(offsets, o);
    }
  }
  return offsets;
}

bool SemilinSetNdd::genepi_init(std::string pluginname, int number_of_variables)
{
  genepi_loader_init();
  genepi_loader_load_default_plugins();
  plugin = genepi_loader_get_plugin(pluginname.c_str());
  std::cout << "Plugin: " << genepi_plugin_get_name(plugin) << std::endl;
  // TODO: use of cache???
  solver = genepi_solver_create(plugin, 0, 0, 0);
  k = number_of_variables;
  return true; // TODO: return correct value!
}

bool SemilinSetNdd::genepi_dealloc() {
  genepi_solver_delete(solver);
  genepi_plugin_del_reference(plugin);
  genepi_loader_terminate();
  return true;
}

SemilinSetNdd::SemilinSetNdd() : set(Genepi(solver, k, false)) {
}

SemilinSetNdd::SemilinSetNdd(int zero) : set(Genepi(solver, k, false)) {
  assert(zero == 0);
}

SemilinSetNdd::SemilinSetNdd(VarId var) : SemilinSetNdd(var, 1) {
}

SemilinSetNdd::SemilinSetNdd(VarId var, int cnt) {
  auto v = var_map.find(var);
  if(v == var_map.end())
  {
    var_map.insert(std::make_pair(var, var_map.size()));
    v = var_map.find(var);
  }

  int position = v->second; // which position is this variable on?
  std::vector<int> alpha(k);
  alpha.at(position) = cnt;

  this->set = Genepi(solver, alpha, false);
  this->offsets.push_back(alpha);
}
SemilinSetNdd::SemilinSetNdd(Genepi set, std::vector<std::vector<int>> offsets) : set(set), offsets(offsets){
}

SemilinSetNdd::SemilinSetNdd(const SemilinSetNdd& expr) : set(expr.set), offsets(expr.offsets) {
}

SemilinSetNdd::~SemilinSetNdd() {
}

SemilinSetNdd SemilinSetNdd::operator=(const SemilinSetNdd& term) {
  this->set = term.set;
  this->offsets = term.offsets;
  return *this;
}

SemilinSetNdd SemilinSetNdd::operator+=(const SemilinSetNdd& term) {
  this->set = this->set.union_op(term.set);
  for(auto offset : term.offsets)
    insert_offset(this->offsets, offset);
  return *this;
}


SemilinSetNdd SemilinSetNdd::operator*=(const SemilinSetNdd& term) {
  std::vector<int> sel_a(3*k, 1);
  std::vector<int> sel_b(3*k, 1);
  std::vector<int> sel_res(3*k, 1);
  for(int i = 0; i < k; i++)
  {
    sel_a[3*i]     = 0;
    sel_b[3*i+1]   = 0;
    sel_res[3*i+2] = 0;
  }

  // inverse projection from original dimensions to an extended version, which is intersected with the natural numbers
  Genepi a_ext(this->set.invproject(sel_a).intersect(Genepi(solver, 3*k, true)));
  Genepi b_ext(term.set.invproject(sel_b).intersect(Genepi(solver, 3*k, true)));

  std::vector<int> generic_alpha = {1, 1, -1}; // 1*a[i]+1*b[i]-1*c[i]=0
  Genepi generic_sum(solver, generic_alpha, 0);

  std::vector<int> component_selection(3*k, 1);

  // initialize the result automaton with (a_i,b_i,N) for i in 0..k-1
  Genepi result(solver, 3*k, true); // natural numbers
  result = result.intersect(a_ext).intersect(b_ext);

  // one run for each variable
  for(int i = 0; i < k; i++)
  {
    // TODO: figure out, if precalculating this and apply inv_project has better complexity
    // than always create a new sum automaton for different components
    component_selection[3*i]   = 0;
    component_selection[3*i+1] = 0;
    component_selection[3*i+2] = 0;
    // inverse projection on generic sum automaton to create automaton for component i
    Genepi component_sum(generic_sum.invproject(component_selection));
    // use this sum automaton on the intermediate result
    result = result.intersect(component_sum);

    component_selection[3*i]   = 1; // reset component_selection
    component_selection[3*i+1] = 1;
    component_selection[3*i+2] = 1;
  }
  this->offsets = multiply_offsets(this->offsets, term.offsets);
  this->set = result.project(sel_res);
  return *this;

}

bool SemilinSetNdd::operator < (const SemilinSetNdd& term) const {
  return this->set < term.set;
}
bool SemilinSetNdd::operator == (const SemilinSetNdd& term) const {
  return this->set == term.set;
}
SemilinSetNdd SemilinSetNdd::star () const {
  SemilinSetNdd offset_star = one();
  for(auto offset : this->offsets) {
    offset_star *= SemilinSetNdd(Genepi(this->solver, offset, true),{std::vector<int>(k,0)}); // offsets are used and converted to generators
  }

  SemilinSetNdd result = one(); // result = 1

  SemilinSetNdd temp = one();
  for(int i = 1; i <= k; i++) { // 1..k
    temp *= *this; // temp = S^i
    result += temp;
  }
  result = one() + (result * offset_star); // S^i * offset_star
  return result;
}

SemilinSetNdd SemilinSetNdd::null() {
  if(!SemilinSetNdd::elem_null)
    SemilinSetNdd::elem_null = std::shared_ptr<SemilinSetNdd>(new SemilinSetNdd());
  return *SemilinSetNdd::elem_null;
}

SemilinSetNdd SemilinSetNdd::one() {
  if(!SemilinSetNdd::elem_one)
    SemilinSetNdd::elem_one = std::shared_ptr<SemilinSetNdd>(new SemilinSetNdd(Genepi(solver, std::vector<int>(k,0), false), {std::vector<int>(k,0)}));
  return *SemilinSetNdd::elem_one;
}

std::string serialize_offsets(const std::vector<std::vector<int>>& offsets) {
  std::stringstream result;
  result << "[";
  for(const auto& offset : offsets) {
    if(&offset != &offsets.at(0))
      result << ",";
    result << "(";
    for(const auto& o : offset) {
      if(&o != &offset.at(0))
        result << ",";
      result << o;
    }
    result << ")";
  }
  result << "]";
  return result.str();
}

std::string SemilinSetNdd::string() const {
  std::stringstream result;
  result << "explicit offsets:\t" << serialize_offsets(this->offsets) << std::endl;
  auto calculated_offsets = this->set.getOffsets();
  result << "calculated offsets:\t" << serialize_offsets(calculated_offsets) << std::endl;
  result << "automaton written:\t" << this->set.output("result", static_i++, "") << std::endl;
  return result.str();
}

const bool SemilinSetNdd::is_idempotent = true;
const bool SemilinSetNdd::is_commutative = true;
std::shared_ptr<SemilinSetNdd> SemilinSetNdd::elem_null;
std::shared_ptr<SemilinSetNdd> SemilinSetNdd::elem_one;
genepi_solver* SemilinSetNdd::solver;
genepi_plugin* SemilinSetNdd::plugin;
std::unordered_map<VarId, int> SemilinSetNdd::var_map;
int SemilinSetNdd::k = 0;
