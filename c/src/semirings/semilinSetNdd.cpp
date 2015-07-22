#include <iostream>
#include "semilinSetNdd.h"
#include "semilinear_util.h"

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

// return all linear independent offsets
// Important: we have to keep the [0,…,0] vector if it is there
std::vector<std::vector<int>> SemilinSetNdd::getIndependentOffsets(const std::vector<std::vector<int>>& offsets) const {
  std::vector<SparseVec<VarId, Counter, DummyDivider>> offset_;
  bool has_zero = false;

  for(auto offset:offsets) {
    std::vector<std::pair<VarId, Counter>> tmp;
    for(auto var_pos : var_map) {
      VarId varid = var_pos.first;
      Counter value = offset.at(var_pos.second);
      if(value != 0) {
        tmp.push_back({varid, value});
      }
    }
    if(tmp.size() != 0) {
      SparseVec<VarId, Counter, DummyDivider> tmp_(std::move(tmp));
      offset_.push_back(tmp_);
    }
    else {
      has_zero = true;
    }
  }

  if(offset_.size() == 0){
    if(has_zero) {
      return {std::vector<int>(k,0)}; // add the zero-vector back if it was there!
    }
    else {
      return {};
    }
  }

  VecSet<SparseVec<VarId, Counter, DummyDivider>> offset_vec_set{std::move(offset_)};

  SimplifySet<SparseVecSimplifier<VarId,Counter,DummyDivider>>(offset_vec_set);

  std::vector<std::vector<std::pair<VarId, Counter>>> temp_result;

  for(auto offset : offset_vec_set) {
    temp_result.push_back(offset.getVector());
  }

  std::vector<std::vector<int>> returned_result;

  for(auto offset : temp_result) {
    std::vector<int> o(k,0);
    for(auto map : offset) {
      o.at( var_map.find(map.first)->second ) = map.second;
    }
    returned_result.push_back(o);
  }

  // if [0,…,0] was in the set before, add it back -- otherwise incorrect (see test-semilinSetNdd.cpp)!
  if(has_zero) {
    std::vector<int> o(k,0);
    returned_result.push_back(o);
  }
  return returned_result;
}

// init/dealloc to main, here: just init/dealloc solver!
bool SemilinSetNdd::genepi_init() {
  genepi_loader_init();
  genepi_loader_load_default_plugins();
  plugin = genepi_loader_get_plugin("mona");
  //std::cout << "Plugin: " << genepi_plugin_get_name(plugin) << std::endl;
  return true;
}

bool SemilinSetNdd::genepi_dealloc() {
  genepi_plugin_del_reference(plugin);
  genepi_loader_terminate();
  return true;
}

bool SemilinSetNdd::solver_init(int number_of_variables)
{
  // TODO: use of cache???
  solver = genepi_solver_create(plugin, 0, 0, 0);
  k = number_of_variables;
  var_map = std::unordered_map<VarId, int>();
  //SemilinSetNdd::elem_null = std::shared_ptr<SemilinSetNdd>(new SemilinSetNdd());
  //SemilinSetNdd::elem_one = std::shared_ptr<SemilinSetNdd>(new SemilinSetNdd(Genepi(solver, std::vector<int>(k,0), false), {std::vector<int>(k,0)}));
  return true; // TODO: return correct value!
}

bool SemilinSetNdd::solver_dealloc() {
  genepi_solver_delete(solver);
  solver=nullptr;
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

SemilinSetNdd::SemilinSetNdd(Genepi set, std::vector<std::vector<int>> offsets) : set(set), offsets(offsets) {
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
  for(auto offset : term.offsets) {
    insert_offset(this->offsets, offset);
  }

  this->offsets = this->getIndependentOffsets(this->offsets);
  return *this;
}

SemilinSetNdd SemilinSetNdd::operator*=(const SemilinSetNdd& term) {
  std::vector<int> sel_a(3*k, 1);
  std::vector<int> sel_b(3*k, 1);
  for(int i = 0; i < k; i++)
  {
    sel_a[3*i]     = 0;
    sel_b[3*i+1]   = 0;
  }

  // inverse projection from original dimensions to an extended version, which is intersected with the natural numbers
  Genepi a_ext(this->set.invproject(sel_a).intersect(Genepi(solver, 3*k, true)));
  Genepi b_ext(term.set.invproject(sel_b).intersect(Genepi(solver, 3*k, true)));

  std::vector<int> generic_alpha = {1, 1, -1}; // 1*a[i]+1*b[i]-1*c[i]=0
  Genepi generic_sum(solver, generic_alpha, 0);

  std::vector<int> component_selection(3*k, 1); // TODO: maybe we only use k+2 or so elements
  std::vector<int> component_projection(3*k, 0); // we want to project away component 3i and 3i+1

  // initialize the result automaton with (a_i,b_i,N) for i in 0..k-1
  Genepi result(solver, 3*k, true); // natural numbers
  result = result.intersect(a_ext).intersect(b_ext);

  // one run for each variable
  // at the end of each loop run we delete the used a and b component
  // the new sum component will be at position i at the end of the loop run
  // (a,b,a+b) component relevant for current run is at position (i,i+1,i+2) at the begin of a run
  for(int i = 0; i < k; i++)
  {
    component_selection.resize(3*k-2*i);
    component_projection.resize(3*k-2*i);
    // TODO: figure out, if precalculating this and apply inv_project has better complexity
    // than always create a new sum automaton for different components
    component_selection[i]   = 0;
    component_selection[i+1] = 0;
    component_selection[i+2] = 0;
    // inverse projection on generic sum automaton to create automaton for component i
    Genepi component_sum(generic_sum.invproject(component_selection));
    // use this sum automaton on the intermediate result
    result = result.intersect(component_sum);

    component_selection[i]   = 1; // reset component_selection
    component_selection[i+1] = 1;
    component_selection[i+2] = 1;

    component_projection[i]   = 1; // project away the already used a and b component at position (i,i+1)
    component_projection[i+1] = 1;
    result = result.project(component_projection);
    component_projection[i]   = 0;
    component_projection[i+1] = 0;
  }
  this->set = result;
  this->offsets = multiply_offsets(this->offsets, term.offsets);
  this->offsets = this->getIndependentOffsets(this->offsets);
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
  // INFO: this is the point where we use the saved (minimal base) of offsets
  // if our simplifier kills null-vectors we loose information (0^*=1...)
  // and our loop will never reach the fixed point
  for(auto offset : this->offsets) {
    offset_star *= SemilinSetNdd(Genepi(this->solver, offset, true),{std::vector<int>(k,0)});
  }

  SemilinSetNdd result = one(); // result = 1

  SemilinSetNdd final_result = one();
  SemilinSetNdd final_result_old;

  SemilinSetNdd S_i = one();


//  do {
//    final_result_old = final_result;
//
//    S_i *= *this;  // = S^i
//    result += S_i; // result = sum_{i=0}^n S^i
//    final_result = one() + (result * offset_star); // 1 + (sum_{i=0}^n S^i) * offset_star
//  } while (final_result != final_result_old);

  //Alternative implementation based on repeated squaring (which is sound since we have idempotent addition!)

  //std::cout << "size (offset-star): " << offset_star.set.getSize() << std::endl;
  //std::cout << "offsets: " << this->offsets << std::endl;
  S_i = one() + *this;
  final_result = one() + *this*offset_star;
  do {
    final_result_old = final_result;
    //std::cout << "size (intermediate automaton -- star-loop-begin): " << final_result.set.getSize() << std::endl;
    //std::cout << "size (S_i -- star-loop-begin): " << S_i.set.getSize() << std::endl;
    S_i *= S_i;
    final_result = S_i * offset_star;
    //std::cout << "size (intermediate automaton -- star-loop-end): " << final_result.set.getSize() << std::endl;
  } while (final_result != final_result_old);

  return final_result;
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


std::string SemilinSetNdd::string() const {
  std::stringstream result;
  // auto calculated_offsets = this->offsets;
  // result << "calculated offsets:\t\t" << serialize_offsets(calculated_offsets) << std::endl;
  // auto cleaned_calculated_offsets = this->getIndependentOffsets(calculated_offsets);
  // result << "calculated offsets (clean):\t" << serialize_offsets(cleaned_calculated_offsets) << std::endl;
  result << std::endl;
  result << "size:\t" << this->set.getSize() << std::endl;
  // result << "automaton written:\t" << this->set.output("result", static_i++, "") << std::endl;
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
