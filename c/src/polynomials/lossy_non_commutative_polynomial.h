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
#include <map>
#include <string>
#include <queue>
#include <climits>


#include "non_commutative_monomial.h"
#include "non_commutative_polynomial.h"

#include "../semirings/semiring.h"
#include "../semirings/lossy-finite-automaton.h"

#include "../datastructs/var.h"

#include "../utils/string_util.h"

class LossyFiniteAutomaton;
class NonCommutativePolynomial<LossyFiniteAutomaton>;


template <>
class NonCommutativePolynomialBase<LossyFiniteAutomaton> : Semiring<NonCommutativePolynomial<LossyFiniteAutomaton>,
                                                                    Commutativity::NonCommutative,
                                                                    LossyFiniteAutomaton::GetIdempotence()> {

protected:

  std::map<NonCommutativeMonomial<LossyFiniteAutomaton>,std::uint_fast16_t> monomials_;

  static void InsertMonomial(
      std::map<NonCommutativeMonomial<LossyFiniteAutomaton>, std::uint_fast16_t> &monomials,
      NonCommutativeMonomial<LossyFiniteAutomaton> monomial,
      std::uint_fast16_t coeff)
  {
    auto iter = monomials.find(monomial);
    if(iter == monomials.end()) {
      monomials.insert({monomial, 1});
    } else {
      auto tmp = *iter;
      monomials.erase(iter);
      monomials.insert({tmp.first, tmp.second + coeff});
    }
  }

  static void InsertMonomial(
      std::map<NonCommutativeMonomial<LossyFiniteAutomaton>, std::uint_fast16_t> &monomials,
      std::vector<std::pair<elemType,int>> &idx,
      std::vector<VarId> &variables,
      std::vector<SR> &srs
  )
  {
    auto tmp_monomial = NonCommutativeMonomial<LossyFiniteAutomaton>(idx, variables, srs);
    InsertMonomial(monomials, tmp_monomial, 1);
  };


public:

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
   void findConstantsInNonterminalMonomials(std::set<LossyFiniteAutomaton> &constants, bool checkLinearTerms) const {

     // delegate to the monomials
     for(auto &monomial: monomials_) {
       monomial.first.findConstantsInNonterminalMonomials(constants, checkLinearTerms);
     }
   }


   LossyFiniteAutomaton sumOfConstantMonomials() const {
     LossyFiniteAutomaton sum = LossyFiniteAutomaton::null();

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
   NonCommutativePolynomial<LossyFiniteAutomaton> replaceConstants(std::map<LossyFiniteAutomaton, VarId> &constantsToVariables,
       bool replaceInLinearMonomials) const {
     NonCommutativePolynomial<LossyFiniteAutomaton> temp = null();

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
   NonCommutativePolynomial<LossyFiniteAutomaton> removeEpsilonMonomials() const {
     NonCommutativePolynomial<LossyFiniteAutomaton> temp = null();

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
   static void removeEpsilonProductions(std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &equations) {
     for(auto equation : equations) {
       equation.second = equation.second.removeEpsilonMonomials();
     }
   }



   /*
    * Cleans the polynomial system, i.e. removes variables that are unproductive or unreachable
    * from the set of nonterminals given in the worklist; the worklist needs to contain the
    * variables we want to start derivations from, i.e. the set of initial nonterminals.
    */
   static std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> cleanSystem
   (const std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &equations,
       std::queue<VarId> &worklist) {

     std::queue<VarId> temp; // to store the variables that still need checking
     std::map<VarId, int> encounteredVariables; // to remember whether a variable was encountered and check for
     // for being productive yet;
     // 0 means "never encountered",
     // 1 means "encountered but not checked",
     // 2 means "checked at least once"
     std::map<VarId, bool> productiveVariables; // to map variables to "is this variable known to be productive?"
     std::map<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>> productions; // to map variables to their productions

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
     std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> cleanEquations;
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
   static std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> eliminateTerminalsInNonterminalProductions
   (const std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &equations,
       std::map<VarId, LossyFiniteAutomaton> &variablesToConstants, bool eliminateInLinearTerms) {

     std::set<LossyFiniteAutomaton> constants; // will hold all terminals appearing in nonterminal productions

     // find all constants in nonterminal productions in the system
     for(auto &equation: equations) {
       equation.second.findConstantsInNonterminalMonomials(constants, eliminateInLinearTerms);
     }

     // return value
     std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> variablefiedEquations;

     // introduce a new variable for each constant found, add the respective equation to the system
     VarId var;
     std::map<LossyFiniteAutomaton, VarId> constantsToVariables; // will hold a mapping from terminals to variables that produce them
     for(auto &constant: constants) {
       var = Var::GetVarId();
       variablefiedEquations.push_back(std::make_pair(var, NonCommutativePolynomial<LossyFiniteAutomaton>(constant)));
       constantsToVariables.insert(std::make_pair(constant, var));
       variablesToConstants.insert(std::make_pair(var,constant));
     }

     // replace the constants in each nonterminal production with the new constant variables
     NonCommutativePolynomial<LossyFiniteAutomaton> allVariablesPoly;
     for(auto &equation: equations) {
       allVariablesPoly = equation.second.replaceConstants(constantsToVariables, eliminateInLinearTerms);
       variablefiedEquations.push_back(std::make_pair(equation.first, allVariablesPoly));
     }

     return variablefiedEquations;
   }

   /*
    * Binarizes nonterminal productions the way it is done when constructing a CNF of a grammar.
    */
   static std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> binarizeProductions
   (const std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &equations,
       std::map<VarId, LossyFiniteAutomaton> &variablesToConstants) {

     // stores mappings between suffixes of monomials and the variables that are introduced
     // during binarization to produce those suffixes
     std::map<std::string, VarId> binarizationVariables;

     // stores the productions (i.e. equations) that are introduced while producing the
     // Chomsky Normal Form
     std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> binarizationVariablesEquations;

     // to store the result
     std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> binarizedEquations;

     // determine the necessary new variables to give all non-terminal productions that need binarizing
     // (i.e. monomials of degree > 2) the form "X = YZ"
     NonCommutativePolynomial<LossyFiniteAutomaton> binarizedPoly;
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
       std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>>> &components,
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
       std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>>> &components,
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
    * See intersectionWithCFG(..) for details.
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
   NonCommutativePolynomial<LossyFiniteAutomaton> intersectionPolynomial(std::vector<unsigned long> &states,
       std::map<unsigned long, std::map<unsigned char, std::forward_list<unsigned long>>> &transitionTable,
       unsigned long &startState,
       unsigned long &targetState,
       std::map<unsigned long, unsigned long> &statesToIndices,
       std::map<VarId, unsigned long> &oldVariablesToIndices,
       std::vector<std::vector<std::vector<VarId>>> &newVariables) const {

     NonCommutativePolynomial<LossyFiniteAutomaton> result = NonCommutativePolynomial<LossyFiniteAutomaton>::null();
     bool epsilonAdded = false;

     // delegate to the monomials
     for(auto &monomial: monomials_) {

       // if the nonterminal can produce epsilon, then the new grammar can produce epsilon only without
       // the FA changing state (since the FA does not have epsilon transitions)
       if(monomial.first.isEpsilonMonomial()) {
         if(startState == targetState && !epsilonAdded) {
           result += NonCommutativePolynomial<LossyFiniteAutomaton>::one();
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
   static LossyFiniteAutomaton shortestWord(std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &productions, VarId &startSymbol) {

     if(productions.size() == 0) {
       return LossyFiniteAutomaton::null();
     }

     // check that the grammar contains the start symbol
     bool startSymbolIsProductive = false;
     for(auto &equation: productions) {
       if(equation.first == startSymbol) {
         startSymbolIsProductive = true;
       }
     }

     if(!startSymbolIsProductive) {
       return LossyFiniteAutomaton::null();
     }

     std::map<VarId, unsigned long> lengthOfShortestWords;
     std::map<VarId, NonCommutativeMonomial<LossyFiniteAutomaton>> productionsForShortestWords;

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
   static std::set<char> getDerivableFirstLettes(const std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &productions, const VarId &startSymbol) {
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
   NonCommutativePolynomial<LossyFiniteAutomaton> removeUnproductiveMonomials(const std::map<VarId, bool> &productiveVariables) const {
     NonCommutativePolynomial<LossyFiniteAutomaton> cleanPoly = null();

     // add all monomials to the clean polynomial that only contain productive variables
     for(auto &monomial: monomials_) {
       if(monomial.first.isProductive(productiveVariables)) {
         InsertMonomial(cleanPoly.monomials_, monomial.first, monomial.second);
       }
     }

     return cleanPoly;
   }


   /*
    * Calculates the intersection of the CFG "equations" with start symbol "oldS" and the
    * regular language represented by a finite automaton. The result is a CFG given by
    * its set of productions (return value of this method) and its start symbol "S".
    *
    * The algorithm is that of Bar-Hillel et al., "On formal properties of simple phrase structure
    * grammars.", 1961; a neat description of it can be found in Nederhof & Satta, "Probabilistic Parsing", 2008
    * which is openly accessible on Nederhof's website as of the writing of this code (spring 2014).
    *
    * libfa doesn't model epsilon transitions (see add_epsilon_trans(state *from, state *to) in fa.c of the
    * source code of libfa for details), instead replacing them by transitions on nonempty symbols and
    * adjusting the set of final states; thus the simplifying assumption of "Probabilistc Parsing" that
    * there are no epsilon transitions in the FA holds, and we only have to deal with the fact that the FA
    * may have multiple final states.
    *
    * WARNING: do not allow productions of the form "X -> empty set" in your oldGrammar, or this may break.
    */

   static std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> intersectionWithFA
       (const FiniteAutomaton &fa, VarId &newS, const VarId &oldS,
               const std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &oldGrammar) {

       // do it with a DFA, something went wrong with NFAs
       FiniteAutomaton minDfa = fa.minimize();

       // for the new grammar
       std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> resultGrammar;

       // we don't want to start out with useless stuff because the intersection grammar
       // will blow up anyway
       std::queue<VarId> worklist;
       worklist.push(oldS);
       auto workGrammar = cleanSystem(oldGrammar, worklist);

       // change the grammar to one where monomials have degree at most 2 and those monomials with degree 2
       // have the form XY
       std::map<VarId, LossyFiniteAutomaton> variablesToConstants;
       workGrammar = eliminateTerminalsInNonterminalProductions
               (workGrammar, variablesToConstants, false);
       workGrammar = binarizeProductions(workGrammar, variablesToConstants);

       // if the grammar doesn't produce anything, there is no need to do anything else; return an empty grammar
       // the same applies if the automaton represents the empty language
       if(workGrammar.size() == 0 || empty()) {
           return resultGrammar;
       }

       /*
        * assuming anyone at all will have to change something here in the future, I expect they
        * won't want to go through the hassle of sifting through the undocumented implementation
        * of libfa, so I'm gonna go ahead and build the transition table of the FA and hope it won't
        * completely destroy our memory limitations; all of this will be represented using number
        * types only, so there's no funny business or trickery to watch out for once the table is built
        */

       // to map the hash of a state "q" to a map that maps each transition symbol "c" occurring at "q"
       // to the set of hashes of target states "t", i.e. t is element of Delta(q,c);
       std::map<unsigned long, std::map<unsigned char, std::forward_list<unsigned long>>> transitionTable;

       // will hold the hashes of all states in the FA; we don't want to deal with with the internals of libfa
       // so all this stuff is represented by simple numbers - hooray for numbers
       std::vector<unsigned long> states;

       // hashes of all final states of the FA
       std::set<unsigned long> finalStates;

       // initial state of the FA
       unsigned long initialState = 0;

       minDfa.extractTransitionTable(transitionTable, states, finalStates, initialState);

       // we will need a three-dimensional lookup table for the new variables so we don't mess up
       // the assignment of polynomials to nonterminals while constructing the grammar
       std::vector<std::vector<std::vector<VarId>>> newVariables;

       // the lookup table will represent states_FA x states_FA x nonterminals_grammar;
       // initialize its dimensions accordingly
       auto numberOfStates = states.size();
       newVariables.resize(numberOfStates);
       for(unsigned long i = 0 ; i < numberOfStates ; ++i) {
           newVariables[i].resize(numberOfStates);

           for(unsigned long j = 0; j < numberOfStates; j++) {
               newVariables[i][j].resize(workGrammar.size());
           }
       }

       // generate the variables of the new grammar, name them <q,X,r> where q and r are states of the DFA and X is
       // a nonterminal of the grammar
       for(unsigned long i = 0; i < numberOfStates; i++) {
           for(unsigned long j = 0; j < numberOfStates; j++) {
               for(unsigned long k = 0; k < workGrammar.size(); k++) {
                   std::stringstream ss;
                   ss << "<" << i << "," << Var::GetVar(workGrammar[k].first).string() << "," << j << ">";
                   newVariables[i][j][k] = Var::GetVarId(ss.str());
               }
           }
       }

       // build a map from variables to indices, where the index of a variable is the line in oldGrammar
       // that contains the productions of the variable
       std::map<VarId, unsigned long> oldVariablesToIndices;
       for(unsigned long i = 0; i < workGrammar.size(); i++) {
           oldVariablesToIndices.insert(std::make_pair(workGrammar[i].first, i));
       }

       // map states to indices; this is just to avoid having to look through the states vector in linear time
       // to find out at which index it lies
       std::map<unsigned long, unsigned long> statesToIndices;
       for(unsigned long i = 0; i < states.size(); i++) {
           statesToIndices.insert(std::make_pair(states[i], i));
       }

       /*
        * now we generate the productions of the new grammar
        */

       // the indices k, i, j are intended this way; we first choose some nonterminal of the grammar
       // and then generate all productions derived from that nonterminal; this mirrors the way
       // in which the algorithm is described in the paper referred to above
       NonCommutativePolynomial<LossyFiniteAutomaton> poly;
       for(unsigned long i = 0; i < numberOfStates; i++) {
           for(unsigned long j = 0; j < numberOfStates; j++) {
               for(unsigned long k = 0; k < workGrammar.size(); k++) {
                   poly = workGrammar[k].second.intersectionPolynomial
                           (states, transitionTable, states[i], states[j], statesToIndices, oldVariablesToIndices, newVariables);
                   resultGrammar.push_back(std::make_pair(newVariables[i][j][k], poly));
               }
           }
       }

       // finally, use a new variable as initial variable and generate the productions; they are of the form
       // newS -> <q_0, oldS, q_f> where q_0 is the initial state of the FA and q_f is some final state of the FA
       NonCommutativePolynomial<LossyFiniteAutomaton> startPolynomial = NonCommutativePolynomial<LossyFiniteAutomaton>::null();

       for(auto target: finalStates) {
           startPolynomial += newVariables[statesToIndices[initialState]][statesToIndices[target]][oldVariablesToIndices[oldS]];
       }

       newS = Var::GetVarId();
       resultGrammar.push_back(std::make_pair(newS, startPolynomial));

       auto startProduction = resultGrammar[resultGrammar.size() - 1];
       resultGrammar[resultGrammar.size() - 1] = resultGrammar[0];
       resultGrammar[0] = startProduction;
       return resultGrammar;
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
   NonCommutativePolynomial<LossyFiniteAutomaton> binarize(std::map<std::string, VarId> &binarizationVariables,
       std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &binarizationVariablesEquations,
       std::map<VarId, LossyFiniteAutomaton> &variablesToConstants) const {
     NonCommutativePolynomial<LossyFiniteAutomaton> binarizedPoly = null();

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
