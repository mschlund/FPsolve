/*
 * semilinSetExp.cpp
 *
 *  Created on: 20.09.2012
 *      Author: maxi
 */

#include <algorithm>
#include <cassert>
#include <iostream>

#ifdef LPSOLVE_OPT
#include <lpsolve/lp_lib.h>
#endif

#include "semilinSetExp.h"

// adding two Var-maps componentwise... could be put in a util-class ?
VecSparse operator+(const VecSparse &a, const VecSparse &b) {
  VecSparse result{a};
  for(auto &pair : b) {
    auto iter_bool = result.insert(pair);
    if (!iter_bool.second) {
      iter_bool.first->second += pair.second;
    }
  }
  return result;
}

// check, if a <= b for two sparse vectors a,b
bool operator<=(const VecSparse &a, const VecSparse &b) {
  if(a.size() > b.size())
    return false;

  // keys of a must be a subset of the keys of b
  // and the respective values must be <=
  for(auto &pair : a) {
    auto it = b.find(pair.first);
    if(it != b.end()) {
      if(pair.second > it->second)
        return false;
    }
    else
      return false;
  }
  return true;
}

std::ostream& operator<<(std::ostream &os, const VecSparse &v) {
  os << "<";
  for (auto &pair : v) {
    os << pair.first << ":" << pair.second << ", ";
  }
  os << ">";
  return os;
}

std::ostream& operator<<(std::ostream &os, const LinSet &ls) {
  os << ls.first;
  for (auto &veclist : ls.second) {
    os << "+" << veclist;
  }
  return os;
}

std::ostream& operator<<(std::ostream &os, const std::set<VecSparse> gens) {
  os << "{";
  for (auto &veclist : gens) {
    os << veclist << ",";
  }
  os << "}";
  return os;
}


// subtracting two Var-maps componentwise (a-b is only defined if b <=a) ... could be put in a util-class ?
VecSparse operator-(const VecSparse &a, const VecSparse &b) {
  assert(b <= a);

  VecSparse result{a};
  for(auto &pair : b) {
    auto it = result.find(pair.first);
    if(it != result.end()){
      it->second -= pair.second;
      if(0 == it->second)
        result.erase(it);
    }
  }
  return result;
}

/*
 * check if b = k*a for a natural number k.
 *  TODO: replace "count > 0" by find (reuse value in subsequent checks)
 */
bool divides(const VecSparse &a, const VecSparse &b) {
  if(a.size() != b.size())
    return false;

  unsigned int k=0;

  for(auto &pair : b) {
    if(a.count(pair.first) > 0) {
      // have we already obtained a candidate for a multiple?
      if(0 != k) {
        if(pair.second != k*a.at(pair.first) )
          return false;
      }
      // no candidate for k yet
      else if(pair.second % a.at(pair.first) != 0)
        return false;
      else
        // division leaves no remainder -> candidate k found
        k = pair.second / a.at(pair.first);
    }
    else {
      // Domains are different! => a does not divide b
      return false;
    }
  }
  return true;
}

// check, if v is a positive integer combination of the vectors in gens
bool is_spanned_by(const VecSparse& v, const std::set<VecSparse>& gens) {
#ifdef LPSOLVE_OPT
  // std::cout << "Check: is " << v << "spanned by " << gens << "?...";

  if(gens.empty()) {
    return false;
  }

  if(v.empty()) {
    return true;
  }

  if(gens.find(v) != gens.end()) {
    //std::cout << "trivially spanned!"<<std::endl;
    return true;
  }

  lprec *lp;
  unsigned int N = gens.size();
  // dim =  number of variables that "gens" is talking about (altogether)

  //collect the variables of all the sparse vectors (v and gens) and number them
  unsigned int num_vars = 1; //we start counting at 1, since the 0-th row will be the values of the objective function
  std::map<VarId, unsigned int> vars;

  for(auto &pair : v) {
    if(vars.find(pair.first) == vars.end()) {
      vars[pair.first] = num_vars;
      num_vars++;
    }
  }

  for(auto &g : gens) {
    for(auto &pair : g) {
      if(vars.find(pair.first) == vars.end()) {
        vars[pair.first] = num_vars;
        num_vars++;
      }
    }
  }

  unsigned int dim = num_vars-1;
  lp = make_lp(dim, N); //first row = objective function (will be zeros in our case)!

  if(lp==NULL) {
    std::cerr << "ERROR: could not create LP (make_lp) of size " << dim+1 << "," << N << std::endl;
    return false;
  }

  set_verbose(lp, IMPORTANT);

  int colno = 1;
  // set the columns of the LP to the generators
  for (auto &g : gens) {
    int num_entries = g.size();

    //std::cout << "num_entries: " << num_entries << std::endl;

    REAL *sparsecolumn = new REAL[num_entries]; //non-zero-vals
    int *rowno = new int[num_entries]; //numbers of the non-zero rows

    int idx = 0;
    //assemble the columns and write them to the LP
    for(auto &pair : g) {
        sparsecolumn[idx] = pair.second;
        rowno[idx] = vars[pair.first];
        idx++;
    }
    set_columnex(lp, colno, num_entries, sparsecolumn, rowno);
    colno++;

    delete[] sparsecolumn;
    delete[] rowno;
  }

  // set the rhs to v
  REAL *rhs = new REAL[dim+1];
  for(auto pair : v) {
    rhs[vars[pair.first]] = pair.second;
  }
  set_rh_vec(lp, rhs);
  delete[] rhs;

  // set all variables to be integers and >=0
  for(int i=1; i<=N; i++){
    set_int(lp, i, TRUE);
    // set_lowbo(lp, i, 0); // the default should be 0 anyways
  }

  for(int i=1; i<=dim; i++){
      set_constr_type(lp, i, EQ);
  }


  set_maxim(lp);

  //print_lp(lp);

  int ret = solve(lp);

  if(ret == 0) {
    // std::cout << "yes!" << std::endl;
    //feasible solution found
/*    REAL *row = new REAL[N];
    get_variables(lp, row);
    for(int j = 0; j < N; j++)
      printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);

    delete[] row;
*/
    return true;
  }
  else if(ret == 2) {
    // std::cout << "no!" << std::endl;
    // LP is infeasible
    return false;
  }
  else {
    std::cerr << "ERROR: lp_solve returned: " << ret << std::endl;
    return false;
  }
#else
  return false;
#endif

}


// simple inclusion check: just check if offset of ls1 is reachable in ls2 and if generators of ls1
// are included (set inclusion) in those of ls2
// TODO: even simpler: check if offsets are the SAME and if generators are subsets
bool is_included_simple(const LinSet &ls1, const LinSet &ls2) {

  //std::cout << "incl check: " << ls1 << " and " << ls2 <<std::endl;

  // definitely a necessary condition
  if(not(ls2.first <= ls1.first))
    return false;

  // here, some first approximation happens. This check is sound of course, but not complete.
  if(!std::includes(ls2.second.begin(), ls2.second.end(), ls1.second.begin(), ls1.second.end())) {
    // std::cout << ls1 << " is NOT (simply) included in " << ls2 << std::endl;
    return false;
  }

  /*
   *  now that we are sure that generators of ls1 are a subset of those of ls2
   *  let's check if the offset of ls1 can be reached from the offset of ls2 (this amounts to solving an ILP!)
  */
  // std::cout << "ILP-based check.. ";
  if (is_spanned_by(ls1.first - ls2.first, ls2.second)) {
    // std::cout << ls1 << " IS (simply) included in " << ls2 << std::endl;
    return true;
  }
  else {
    // std::cout << ls1 << " IS (simply) included in " << ls2 << std::endl;
    return false;
  }
}


// check for inclusion between linsets via solving an ILP
bool is_included(LinSet &ls1, LinSet &ls2) {
  assert(false);
  // TODO: implement... :)
  return false;
}

/*
 * compute the (unique, minimal) Hilbert-basis for the linear set and replace the generators by the basis-elements
 * (this only reduces the size of the generator-set!)
 *
*/
void min_hilbert(LinSet &ls) {
  // TODO: implement... :)
  assert(false);
}

/*
 * Perform simple optimization:
 * if h = k*g for some generators g!=h and natural number k, then remove h from generators
 */
void clean_generators_simple(LinSet &ls) {
  // check for all pairs of generators (g,h) if divides(h,g)
  for(auto g = ls.second.begin(); g!=ls.second.end(); ++g) {
    for(auto h = ls.second.begin(); h!=ls.second.end();) {
      if(g==h){
        ++h;
        continue;
      }
      if(divides(*g,*h)) {
        h = ls.second.erase(h);
      }
      else
        ++h;
    }
  }
}

void clean_generators(LinSet &ls) {
  clean_generators_simple(ls);

  VecSparse v;

  for(auto g = ls.second.begin(); g!=ls.second.end();) {
    v=*g;
    g = ls.second.erase(g);
    if(!is_spanned_by(v,ls.second)) {
      ls.second.insert(v);
    }
  }
}

// check for pairwise inclusions between the linsets,
// if inclusion found: remove the subset
void SemilinSetExp::clean_slset() {
  for(auto l = val.begin(); l!=val.end(); ++l) {
    for(auto m = val.begin(); m!=val.end();) {
      if(l==m){
        ++m;
        continue;
      }

      if(is_included_simple(*m,*l)) {
        m = val.erase(m);
      }
      else
        ++m;
    }
  }
}

LinSet operator*(const LinSet &ls1, const LinSet &ls2) {

  LinSet result;

  /* Add the offsets... */
  result.first = ls1.first + ls2.first;

  /* ... and union on the generators */
  std::insert_iterator< std::set<VecSparse> > it(result.second,
                                                 result.second.begin());
  set_union(ls1.second.begin(), ls1.second.end(),
            ls2.second.begin(), ls2.second.end(), it);

  clean_generators(result);
  return result;
}

SemilinSetExp::SemilinSetExp() : val() { }

SemilinSetExp::SemilinSetExp(VarId var) : val() {
  VecSparse offset = { std::make_pair(var, 1) };
  LinSet ls{};
  ls.first = offset;
  val.insert(std::move(ls));
}

SemilinSetExp::SemilinSetExp(VarId var, unsigned int cnt) : val() {
  if(0 != cnt) {
    VecSparse offset = { std::make_pair(var, cnt) };
    LinSet ls{};
    ls.first = offset;
    val.insert(std::move(ls));
  }
  else {
    std::cerr << "[INFO] SL-set: tried to generate slset having a variable with count zero.. ignoring it." << std::endl;
  }
}

SemilinSetExp::SemilinSetExp(const std::set<LinSet> &v) {
  val = v;
}

SemilinSetExp::~SemilinSetExp() {
  // do NOT delete static pointers!!!
}

SemilinSetExp SemilinSetExp::null() {
  if(!SemilinSetExp::elem_null) {
    SemilinSetExp::elem_null =
      std::make_shared<SemilinSetExp>(std::set<LinSet>());
  }
  return *SemilinSetExp::elem_null;
}

SemilinSetExp SemilinSetExp::one() {
  if (!SemilinSetExp::elem_one) {
    std::set<LinSet> elone = { LinSet{} };
    SemilinSetExp::elem_one = std::make_shared<SemilinSetExp>(elone);
  }
  return *SemilinSetExp::elem_one;
}

// TODO: check for obvious inclusions and remove them
SemilinSetExp SemilinSetExp::operator+=(const SemilinSetExp &sl) {
  std::set<LinSet> result;
  std::insert_iterator< std::set<LinSet> > it(result, result.begin());
  std::set_union(val.begin(), val.end(), sl.val.begin(), sl.val.end(), it);
  val = std::move(result);

  //TODO: for testing purposes only right now :)
  clean_slset();

  return *this;
}

SemilinSetExp SemilinSetExp::operator*=(const SemilinSetExp &sl) {
  std::set<LinSet> result;
  for(auto &lin_set_rhs : sl.val) {
    for(auto &lin_set_lhs : val) {
      result.insert(lin_set_rhs * lin_set_lhs);
    }
  }
  val = std::move(result);

  //TODO: just for testing-purposes right now :)
  clean_slset();

  return *this;
}

// TODO: semantic equivalence check or at least some more sophisticated check
bool SemilinSetExp::operator == (const SemilinSetExp &sl) const {
  return (val == sl.val);
}

std::set<LinSet> SemilinSetExp::star(const LinSet &ls) {

  /* If we do not have any generators, i.e.,
   *   ls = w  (for some word w)
   * just return
   *   w*
   * instead of 1 + ww*.  */
  if (ls.second.empty()) {
    LinSet r = ls;
    /* If w is not the one-element, move w to the generators. */
    if (ls.first != VecSparse()) {
      r.second.insert(ls.first);
      r.first = VecSparse();
    }
    std::set<LinSet> res = { r };
    return res;
  }

  SemilinSetExp tmp_one{one()};
  std::set<LinSet> v{tmp_one.val};
  std::set<LinSet> result{v};

  /* Star of a linear set is a semilinear set:
   * (w_0.w_1*.w_2*...w_n*)* = 1 + (w_0.w_0*.w_1*.w_2*...w_n*)
   */

  LinSet ls_tmp = ls;
  ls_tmp.second.insert(ls.first);
  result.insert(std::move(ls_tmp));

  return result;
}

SemilinSetExp SemilinSetExp::star() const {
  SemilinSetExp result = SemilinSetExp::one();
  for (auto &ls : val) {
    result *= SemilinSetExp(star(ls));
  }
  return result;
}

std::string SemilinSetExp::string() const {
  std::stringstream ss;
  ss << "[";
  for (auto &ls : val) {
    ss << ls << "\n";
  }
  ss << "]";
  return ss.str();
}

std::ostream& SemilinSetExp::operator<<(std::ostream &os) const {
  return os << string();
}


const bool SemilinSetExp::is_idempotent = true;
const bool SemilinSetExp::is_commutative = true;
std::shared_ptr<SemilinSetExp> SemilinSetExp::elem_null;
std::shared_ptr<SemilinSetExp> SemilinSetExp::elem_one;
