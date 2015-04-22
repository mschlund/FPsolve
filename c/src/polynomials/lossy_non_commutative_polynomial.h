/*
 * lossy_non_commutative_polynomial.h
 *
 *  Created on: 22.04.2015
 *      Author: schlund
 */

#ifndef LOSSY_NON_COMMUTATIVE_POLYNOMIAL_H_
#define LOSSY_NON_COMMUTATIVE_POLYNOMIAL_H_

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <list>
#include <map>
#include <queue>
#include <string>
#include <climits>


#include "../semirings/lossy-finite-automaton.h"
#include "non_commutative_polynomial.h"


class LossyNonCommutativePolynomial : public NonCommutativePolynomial<LossyFiniteAutomaton> {

public:

  LossyNonCommutativePolynomial(const NonCommutativePolynomial<LossyFiniteAutomaton> & p) {
    for(auto m : p.monomials_) {
      this->monomials_[LossyNonCommutativeMonomial(m.first)] = m.second;
    }
  }


  /*
   * Extracts the terminal letters used in this polynomial; only recognizes alphanumeric characters.
   *
   * Use this as little as possible; cache it if you need it.
   */
  std::set<unsigned char> get_terminals() const {
      std::set<unsigned char> terminals;

      for(auto const &monomial: monomials_) {
          auto monomial_terminals = monomial.first.get_terminals();
          terminals.insert(monomial_terminals.begin(), monomial_terminals.end());
      }

      return terminals;
  }

  /*
   * As with get_variables(), make sure you cache this if you use it.
   *
   * Only considers monomials with exactly two variables; used in generating the lossy quadratic normal form
   * of a system, see LossyFiniteAutomaton::downwardClosureDerivationTrees(..).
   */
  std::set<VarId> get_variables_quadratic_monomials() const {
      std::set<VarId> vars;

      for(auto const &monomial : monomials_) {
          if(monomial.first.get_degree() == 2) {
              auto tmp = monomial.first.get_variables();
              vars.insert(tmp.begin(), tmp.end());
          }
      }

      return vars;
  }

  /*
   * Finds all the constants appearing in nonterminal monomials of this polynomial and appends each new one to the vector.
   *
   * Only checks linear monomials if checkLinearTerms is true.
   */
  void findConstantsInNonterminalMonomials(std::set<SR> &constants, bool checkLinearTerms) const {

      // delegate to the monomials
      for(auto &monomial: monomials_) {
          monomial.first.findConstantsInNonterminalMonomials(constants, checkLinearTerms);
      }
  }


  SR sumOfConstantMonomials() const {
      SR sum = SR::null();

      for(auto &monomial: monomials_) {
          if(monomial.first.get_degree() == 0) {
              sum = sum + monomial.first.getLeadingSR();
          }
      }

      return sum;
  }

  /*
   * Replaces all constants in the nonterminal monomials of the polynomial with their respective VarId mapping.
   * Used in transformation to Chomsky Normal Form.
   *
   * Only replaces constants in linear monomials if replaceInLinearMonomials is true.
   */
  NonCommutativePolynomial<SR> replaceConstants(std::map<SR, VarId> &constantsToVariables,
          bool replaceInLinearMonomials) const {
      NonCommutativePolynomial<SR> temp = null();

      // delegate to the monomials
      for(auto &monomial: monomials_) {

          // we only replace constants in nonterminal monomials, all terminal monomials are added to the polynomial as is
          // only replace constants in linear monomials if replaceInLinearMonomials is true
          if((monomial.first.get_degree() != 0)
                  && (replaceInLinearMonomials || monomial.first.get_degree() != 1)) {
              temp += monomial.first.replaceConstants(constantsToVariables);
          } else {
              InsertMonomial(temp.monomials_, monomial.first, monomial.second);
          }
      }

      return temp;
  }

  /*
   * Removes all epsilon monomials from a polynmoial.
   * An epsilon monomial is a monomial consisting only of epsilon.
   */
  NonCommutativePolynomial<SR> removeEpsilonMonomials() const {
      NonCommutativePolynomial<SR> temp = null();

      for(auto monomial: monomials_) {

          // if the monomial has at least one variable or it doesn't have a variable but isn't
          // equal to epsilon, then it's not an epsilon monomial
          if(!monomial.first.isEpsilonMonomial()) {
              temp += monomial.first * monomial.second;
          }
      }

      return temp;
  }

  /*
   * Removes all epsilon productions from the system.
   */
  static void removeEpsilonProductions(std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &equations) {
      for(auto equation: equations) {
          equation.second = equation.second.removeEpsilonMonomials();
      }
  }



  /*
   * Cleans the polynomial system, i.e. removes variables that are unproductive or unreachable
   * from the set of nonterminals given in the worklist; the worklist needs to contain the
   * variables we want to start derivations from, i.e. the set of initial nonterminals.
   */
  static std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> cleanSystem
      (const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations,
              std::queue<VarId> &worklist) {

      std::queue<VarId> temp; // to store the variables that still need checking
      std::map<VarId, int> encounteredVariables; // to remember whether a variable was encountered and check for
                                                 // for being productive yet;
                                                 // 0 means "never encountered",
                                                 // 1 means "encountered but not checked",
                                                 // 2 means "checked at least once"
      std::map<VarId, bool> productiveVariables; // to map variables to "is this variable known to be productive?"
      std::map<VarId, LossyNonCommutativePolynomial> productions; // to map variables to their productions

      VarId var;
      std::set<VarId> vars;
      // build the data structures that will store the info about which variables still
      // need checking and which variables are known to be productive
      for(auto &equation: equations) {
          vars = equation.second.get_variables();

          // just to make sure we catch all variables in case someone inputs a grammar with a terminal
          // that only appears in a lefthand side of a production
          for(auto &prodVar: vars) {
              encounteredVariables.insert(std::make_pair(prodVar, 0));
              productiveVariables.insert(std::make_pair(prodVar, false));
          }

          encounteredVariables.insert(std::make_pair(equation.first, 0));
          productiveVariables.insert(std::make_pair(equation.first, false));
          productions.insert(std::make_pair(equation.first, equation.second));

          // prepare the info for the initial variables
          while(!worklist.empty()) {
              var = worklist.front();
              worklist.pop();
              encounteredVariables[var] = 1;
              temp.push(var);
          }
          worklist.swap(temp);
      }

      // keep checking the variables that are not yet known to be productive
      // until no further variable becomes known to be productive; we do this
      // using BFS so unreachable variables will be removed as well since they
      // can never be found to be productive this way
      bool update;
      do {
          update = false;

          // terminates once there is nothing more to do; specifically, this terminates at the latest once
          // all reachable variables are known to be productive because then they won't be put in the
          // worklist again; it terminates at the earliest once all reachable variables have been checked
          // without finding a productive one since the worklist is being filled with all reachable variables
          // before the end of the first iteration
          while(!worklist.empty()) {
              var = worklist.front();
              worklist.pop();

              // if we have not yet checked this variable, remember that we now have and put all
              // the unencountered variables in its productions into the worklist; this way, we only
              // put reachable variables into the worklist once
              if(encounteredVariables[var] == 1) {
                  encounteredVariables[var] = 2;

                  for(auto &monomial: productions[var].monomials_) {
                      for(auto &reachableVar: monomial.first.get_variables()) {

                          // only put variables into the worklist that have not been seen before
                          if(encounteredVariables[reachableVar] == 0) {
                              worklist.push(reachableVar);
                              encounteredVariables[reachableVar] = 1;
                          }
                      }
                  }
              }

              // if the variable was found to be productive, remember the update and change the
              // flag for the variable; we know a variable is productive once one of the productions
              // associated with it becomes known to be productive; since there is no need to recheck a
              // productive variable, don't put it in the worklist again
              if(productions[var].isProductive(productiveVariables)) {
                  productiveVariables[var] = true;
                  update = true;
              } else { // if the variable isn't known to be productive yet, put it in the worklist
                       // for the next iteration
                   temp.push(var);
              }
          }

          // don't do work we don't need
          if(update) {
              worklist.swap(temp);
          }
      } while(update);

      // build the clean system: only use variables that are both productive and reachable;
      // remove unproductive monomials
      std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> cleanEquations;
      for(auto &equation: equations) {
          if(productiveVariables[equation.first]) {

              // at this point we know that the polynomial we are associating with the
              // variable is not the null polynomial since at least one of its monomials was
              // found to be productive
              cleanEquations.push_back(std::make_pair(equation.first,
                      equation.second.removeUnproductiveMonomials(productiveVariables)));
          }
      }

      return cleanEquations;
  }

  /*
   * Equations can contain productions that are sentential forms, and variablefiedEquations will afterwards hold
   * productions that are either terminal or only contain nonterminals (i.e. the monomials in each polynomial
   * are either semiring elements or only contain variables).
   *
   * Linear terms (i.e. such with exactly one variable in them) will have their terminals replaced by variables
   * iff eliminateInLinearTerms is true. This way, we can produce a quadratic normal form without  blowing up
   * linear terms, if desired.
   */
  static std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> eliminateTerminalsInNonterminalProductions
              (const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations,
                     std::map<VarId, SR> &variablesToConstants, bool eliminateInLinearTerms) {

      std::set<SR> constants; // will hold all terminals appearing in nonterminal productions

      // find all constants in nonterminal productions in the system
      for(auto &equation: equations) {
          equation.second.findConstantsInNonterminalMonomials(constants, eliminateInLinearTerms);
      }

      // return value
      std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> variablefiedEquations;

      // introduce a new variable for each constant found, add the respective equation to the system
      VarId var;
      std::map<SR, VarId> constantsToVariables; // will hold a mapping from terminals to variables that produce them
      for(auto &constant: constants) {
          var = Var::GetVarId();
          variablefiedEquations.push_back(std::make_pair(var, LossyNonCommutativePolynomial(constant)));
          constantsToVariables.insert(std::make_pair(constant, var));
          variablesToConstants.insert(std::make_pair(var,constant));
      }

      // replace the constants in each nonterminal production with the new constant variables
      LossyNonCommutativePolynomial allVariablesPoly;
      for(auto &equation: equations) {
          allVariablesPoly = equation.second.replaceConstants(constantsToVariables, eliminateInLinearTerms);
          variablefiedEquations.push_back(std::make_pair(equation.first, allVariablesPoly));
      }

      return variablefiedEquations;
  }

  /*
   * Binarizes nonterminal productions the way it is done when constructing a CNF of a grammar.
   */
  static std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> binarizeProductions
              (const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations,
                     std::map<VarId, SR> &variablesToConstants) {

      // stores mappings between suffixes of monomials and the variables that are introduced
      // during binarization to produce those suffixes
      std::map<std::string, VarId> binarizationVariables;

      // stores the productions (i.e. equations) that are introduced while producing the
      // Chomsky Normal Form
      std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> binarizationVariablesEquations;

      // to store the result
      std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> binarizedEquations;

      // determine the necessary new variables to give all non-terminal productions that need binarizing
      // (i.e. monomials of degree > 2) the form "X = YZ"
      LossyNonCommutativePolynomial binarizedPoly;
      for(auto &equation: equations) {
          binarizedPoly =
                  equation.second.binarize(binarizationVariables, binarizationVariablesEquations, variablesToConstants);
          binarizedEquations.push_back(std::make_pair(equation.first, binarizedPoly));
      }

      // finally, add the productions to the system that were introduced during the last step
      for(auto &equation: binarizationVariablesEquations) {
          binarizedEquations.push_back(std::make_pair(equation.first, equation.second));
      }

      return binarizedEquations;
  }

  /*
   * Checks whether this polynomial is productive, depending on the set of variables
   * that are already known to be productive.
   */
  bool isProductive(const std::map<VarId, bool> &productiveVariables) const {

      // check if any monomial in this polynomial is productive; that will be enough
      // for the polynomial to be productive
      for(auto &monomial: monomials_) {
          if(monomial.first.isProductive(productiveVariables)) {
              return true;
          }
      }

      // if there is no productive monomial, then the polynomial isn't known to be productive so far
      return false;
  }

  /*
   * Decides whether the variables in a strongly connected component of a polynomial system can
   * duplicate themselves via this polynomial. In order to decide that, the map "varToComponent"
   * needs to hold the information in which component each variable lies, while "component" is
   * the component we are interested in.
   *
   * Used in the algorithm by Courcelle, see LossyFiniteAutomaton::downwardClosureCourcelle.
   */
  bool componentIsSquarable(std::map<VarId, int> &varToComponent, int component) const {
      bool squarable = false;

      for(auto &monomial: monomials_) {

          // skip monomials that don't have at least two variables
          if(monomial.first.get_degree() >= 2) {
              squarable |= monomial.first.componentIsSquarable(varToComponent, component);

              // don't do work we don't need
              if(squarable) {
                  break;
              }
          }
      }

      return squarable;
  }

  /*
   * Used in the algorithm by Courcelle, see LossyFiniteAutomaton::downwardClosureCourcelle.
   *
   * Assumes that all productions have been binarized.
   */
  void mapQuadraticLHStoRHS(std::map<int, std::map<int, std::set<int>>> &quadraticLHStoRHS,
          std::map<VarId, int> &varToComponent, int component) const {
      for(auto &monomial: monomials_) {
          monomial.first.mapQuadraticLHStoRHS(quadraticLHStoRHS, varToComponent, component);
      }
  }

  /*
   * Used in the algorithm by Courcelle, see LossyFiniteAutomaton::downwardClosureCourcelle.
   *
   * Assumes that all productions have been binarized.
   */
  void calculateLowerComponentVariables(
              std::map<int, std::set<int>> &lhsLowerComponentVariables,
              std::map<int, std::set<int>> &rhsLowerComponentVariables,
              std::vector<std::vector<std::pair<VarId, LossyNonCommutativePolynomial>>> &components,
              std::map<VarId, int> &varToComponent,
              int component) {

      for(auto &monomial: monomials_) {
          monomial.first.calculateLowerComponentVariables(lhsLowerComponentVariables, rhsLowerComponentVariables,
                  components, varToComponent, component);
      }
  }

  /*
   * Used in the algorithm by Courcelle, see LossyFiniteAutomaton::downwardClosureCourcelle.
   *
   * Assumes grammar is in quadratic normal form.
   *
   * WARNING: will break if used with anything but NonCommutativePolynomial<LossyFiniteAutomaton>
   */
  void calculateSameComponentLetters(
          std::map<int, std::set<unsigned char>> &lhsSameComponentLetters,
          std::map<int, std::set<unsigned char>> &rhsSameComponentLetters,
          std::vector<std::vector<std::pair<VarId, LossyNonCommutativePolynomial>>> &components,
          std::map<VarId, int> &varToComponent,
          int component) {

      for(auto &monomial: monomials_) {
          monomial.first.calculateSameComponentLetters(lhsSameComponentLetters, rhsSameComponentLetters,
                  components, varToComponent, component);
      }
  }

  /*
   * Used during the generation of the intersection between a CFG and a FiniteAutomaton.
   * Generates all productions for a triple in states x states x nonterminals, represented
   * by their indices in the vectors "states" and "oldGrammar",
   * where "oldGrammar[nonterminalIndex].first" is the nonterminal we care about.
   *
   * See FiniteAutomaton::intersectionWithCFG(..) for details.
   *
   * WARNING: if you use an SR where the string representations of elements are NOT strings over
   * the English alphabet (e.g. they may contain 1s to represent epsilons, or they may use a
   * different character set), then this function will possibly not do what you want it to do.
   * The problem we run into there is that we have to process all productions of the form
   * A -> X_1 X_2 X_3 ... X_m where the (X_i)s are all elements of either the terminal or the
   * nonterminal alphabet, so having any symbols in the grammar that don't belong to either
   * will give you trouble
   *
   * WARNING: will currently break if used with anything but NonCommutativePolynomial<LossyFiniteAutomaton>
   */
  LossyNonCommutativePolynomial intersectionPolynomial(std::vector<unsigned long> &states,
          std::map<unsigned long, std::map<unsigned char, std::forward_list<unsigned long>>> &transitionTable,
          unsigned long &startState,
          unsigned long &targetState,
          std::map<unsigned long, unsigned long> &statesToIndices,
          std::map<VarId, unsigned long> &oldVariablesToIndices,
          std::vector<std::vector<std::vector<VarId>>> &newVariables) const {

      LossyNonCommutativePolynomial result = LossyNonCommutativePolynomial::null();
      bool epsilonAdded = false;

      // delegate to the monomials
      for(auto &monomial: monomials_) {

          // if the nonterminal can produce epsilon, then the new grammar can produce epsilon only without
          // the FA changing state (since the FA does not have epsilon transitions)
          if(monomial.first.isEpsilonMonomial()) {
              if(startState == targetState && !epsilonAdded) {
                  result += LossyNonCommutativePolynomial::one();
                  epsilonAdded = true;
              }
          } else { // if this is a non-epsilon production, we calculate the new productions that derive from it
              result += monomial.first.intersectionPolynomial
                      (states, transitionTable, startState, targetState,
                              statesToIndices, oldVariablesToIndices, newVariables);
          }
      }

      return result;
  }

  /*
   * Finds a shortest word in the language described by the grammar starting from "startSymbol" and having
   * productions "productions". Will be null if the language is empty. Undefined behavior if startSymbol is
   * not productive in the grammar defined by productions - "productions" must define a clean grammar starting
   * from "startSymbol".
   */
  static SR shortestWord(std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &productions, VarId &startSymbol) {

      if(productions.size() == 0) {
          return SR::null();
      }

      // check that the grammar contains the start symbol
      bool startSymbolIsProductive = false;
      for(auto &equation: productions) {
          if(equation.first == startSymbol) {
              startSymbolIsProductive = true;
          }
      }

      if(!startSymbolIsProductive) {
          return SR::null();
      }

      std::map<VarId, unsigned long> lengthOfShortestWords;
      std::map<VarId, NonCommutativeMonomial<SR>> productionsForShortestWords;

      for(auto &production: productions) {
          lengthOfShortestWords.insert(std::make_pair(production.first, ULONG_MAX));
      }

      // while some shorter derivable word was found, check whether we can now find some other shorter
      // derivable word; terminates once we know a shortest derivable word for every nonterminal
      bool update;
      do {
          update = false;

          for(auto &equation: productions) {
              for(auto &monomial: equation.second.monomials_) {
                  update |= monomial.first.findLengthOfDerivableStrings(lengthOfShortestWords, productionsForShortestWords, equation.first);
              }
          }
      } while(update);

      assert(productionsForShortestWords.find(startSymbol) != productionsForShortestWords.end());
      return productionsForShortestWords[startSymbol].shortestDerivableTerminal(productionsForShortestWords);
  }

  /*
   * Finds all terminal symbols that can be the first letter of some word in the language of the grammar.
   */
  static std::set<char> getDerivableFirstLettes(const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &productions, const VarId &startSymbol) {
      std::map<VarId, std::set<char>> varToDerivableFirstLetters;
      std::set<VarId> varsWithEpsilon;
      if(productions.size() != 0) {

          // while we find new derivable first letters, continue updating
          bool update;
          do {
              update = false;

              for(auto &equation: productions) {
                  for(auto &monomial: equation.second.monomials_) {
                      update |= monomial.first.getDerivableFirstLettes(varToDerivableFirstLetters, varsWithEpsilon, equation.first);
                  }
              }
          } while(update);
      }

      return varToDerivableFirstLetters[startSymbol];
  }

  /*
   * Returns a cleaned version of this polynomial, i.e. a version that had all monomials with unproductive
   * variables eliminated.
   */
  LossyNonCommutativePolynomial removeUnproductiveMonomials(const std::map<VarId, bool> &productiveVariables) const {
      LossyNonCommutativePolynomial cleanPoly = null();

      // add all monomials to the clean polynomial that only contain productive variables
      for(auto &monomial: monomials_) {
          if(monomial.first.isProductive(productiveVariables)) {
              InsertMonomial(cleanPoly.monomials_, monomial.first, monomial.second);
          }
      }

      return cleanPoly;
  }

private:

  /*
   * Binarizes this polynomial, as one would when calculating a CNF of a grammar. Make sure all monomials
   * either have no terminals (i.e. semiring elements) or no nonterminals (i.e. variables) as factors.
   *
   * If a monomial consists of exactly one variable, it will not be changed (chain productions are not eliminated).
   *
   * Make sure there are no productions with both terminals and nonterminals as factors; those will be eliminated if
   * they have at least 3 variables as factors. Call eliminateTerminalsInNonterminalProductions(..) first to avoid this.
   *
   * The variables and productions that are introduced during that process will be stored in binarizationVariables
   * and binarizationVariablesEquations, respectively.
   */
  LossyNonCommutativePolynomial binarize(std::map<std::string, VarId> &binarizationVariables,
          std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &binarizationVariablesEquations,
          std::map<VarId, LossyFiniteAutomaton> &variablesToConstants) const {
    LossyNonCommutativePolynomial binarizedPoly = null();

      // delegate to the monomials
      for(auto &monomial: monomials_) {
          if(monomial.first.get_degree() > 2) { // only binarize productions that have at least 3 variables in them
              binarizedPoly += monomial.first.binarize(binarizationVariables, binarizationVariablesEquations, variablesToConstants);
          } else {
              InsertMonomial(binarizedPoly.monomials_, monomial.first, monomial.second);
          }
      }

      return binarizedPoly;
  }



};




#endif /* LOSSY_NON_COMMUTATIVE_POLYNOMIAL_H_ */
