/*
 * lossy_non_commutative_monomial.h
 *
 *  Created on: 22.04.2015
 *      Author: schlund
 */

#ifndef LOSSY_NON_COMMUTATIVE_MONOMIAL_H_
#define LOSSY_NON_COMMUTATIVE_MONOMIAL_H_

#include <map>
#include <set>
#include <climits>
#include <forward_list>


#include "../datastructs/equations.h"
#include "non_commutative_polynomial.h"
#include "../semirings/lossy-finite-automaton.h"


template <>
class NonCommutativeMonomial<LossyFiniteAutomaton> : public NonCommutativeMonomialBase<LossyFiniteAutomaton>{


public:
  using NonCommutativeMonomialBase<LossyFiniteAutomaton>::NonCommutativeMonomialBase;

  NonCommutativeMonomial() = default;


  // constructor for implicit conversions
  NonCommutativeMonomial(const NonCommutativeMonomialBase<LossyFiniteAutomaton>& p) {
    this->idx_ = p.idx_;
    this->srs_ = p.srs_;
    this->variables_ = p.variables_;
  }

  /*
   * Extracts the terminal letters used in this monomial. Currently only recognizes alphanumeric characters.
   */
  std::set<unsigned char> get_terminals() const {
    std::set<unsigned char> terminals;

    for(auto &semiring_element: srs_) {
      if(!(semiring_element == LossyFiniteAutomaton::null()) && !(semiring_element == LossyFiniteAutomaton::one())) {
        std::string regex = semiring_element.string();

        for(int i = 0; i < regex.size(); i++) {
          if(isalnum(regex[i])) {
            terminals.insert(regex[i]);
          }
        }
      }
    }

    return terminals;
  }

  /*
   * Finds all the constants appearing in this monomial if it has degree at least 1 if checkLinearMonomials is true,
   * otherwise if it has degree at least 2.
   */
  void findConstantsInNonterminalMonomials(std::set<LossyFiniteAutomaton> &constants, bool checkLinearMonomials) const {


    // we only replace the constants in monomials of length >= 2 with new variables; we also ignore linear
    // monomials depending on checkLinearMonomials
    if(get_degree() == 0 || (!checkLinearMonomials && get_degree() == 1)
        || (get_degree() == 1 && getLeadingSR() == LossyFiniteAutomaton::one() && getTrailingSR() == LossyFiniteAutomaton::one())) {
      return;
    }

    for(int i = idx_.size() - 1; i >= 0; i--) {
      //we only want to extract the constants from the monomial; we also don't care about epsilon
      if((idx_.at(i).first == SemiringType) && !(srs_.at(idx_.at(i).second) == LossyFiniteAutomaton::one())) {

        //std::set.insert(..) does its own checking whether the element already is in the set
        constants.insert(srs_.at(idx_.at(i).second));
      }
    }

    return;
  }

  /*
   * Replaces all constants in the monomial with their respective VarId mapping.
   *
   * WARNING: This function will fail if there is an SR element in this monomial that does
   * not have a mapped VarId. Will not change the monomial consisting exclusively of epsilon.
   */
  NonCommutativePolynomialBase<LossyFiniteAutomaton> replaceConstants(std::map<LossyFiniteAutomaton, VarId> &constantsToVariables) const {
    NonCommutativePolynomialBase<LossyFiniteAutomaton> temp;

    bool firstFactorIsSR = (idx_.at(0).first == SemiringType);
    bool firstFactorIsEpsilon = false;

    if(firstFactorIsSR) {
      firstFactorIsEpsilon = (srs_.at(idx_.at(0).second) == LossyFiniteAutomaton::one());
    }

    if((idx_.size() == 1) && firstFactorIsSR && firstFactorIsEpsilon) {
      return NonCommutativePolynomialBase<LossyFiniteAutomaton>::one();
    }

    int startIndex = 1;
    if(firstFactorIsSR) {

      // we don't want unnecessary epsilon-factors
      if(firstFactorIsEpsilon) {
        temp = NonCommutativePolynomialBase<LossyFiniteAutomaton>(variables_.at(idx_.at(1).second));
        startIndex = 2;
      } else {
        temp = NonCommutativePolynomialBase<LossyFiniteAutomaton>(constantsToVariables[srs_.at(idx_.at(0).second)]);
      }
    } else {
      temp = NonCommutativePolynomialBase<LossyFiniteAutomaton>(variables_.at(idx_.at(0).second));
    }

    for(int i = startIndex; i < idx_.size(); i++) {

      // multiply as long as there was no variable encountered
      if(idx_.at(i).first == SemiringType) {

        // ignore all epsilon factors
        if(!(srs_.at(idx_.at(i).second) == LossyFiniteAutomaton::one())) {
          temp *= NonCommutativePolynomialBase<LossyFiniteAutomaton>(constantsToVariables[srs_.at(idx_.at(i).second)]);
        }
      } else {
        temp *= variables_.at(idx_.at(i).second);
      }
    }

    return temp;
  }


  /*
   * This produces the polynomial (that derives from this monomial) that is generated while
   * building the intersection grammar of a CFG and an FA.
   *
   * See NonCommutativePolynomialBase<LossyFiniteAutomaton>::intersectionPolynomial(..) for further documentation.
   *
   */
  NonCommutativePolynomialBase<LossyFiniteAutomaton> intersectionPolynomial(std::vector<unsigned long> &states,
      std::map<unsigned long, std::map<unsigned char, std::forward_list<unsigned long>>> &transitionTable,
      unsigned long startState,
      unsigned long targetState,
      std::map<unsigned long, unsigned long> &statesToIndices,
      std::map<VarId, unsigned long> &oldVariablesToIndices,
      std::vector<std::vector<std::vector<VarId>>> &newVariables) const {

    return (NonCommutativePolynomialBase<LossyFiniteAutomaton>::one() *
        generateIntersectionMonomials(states, transitionTable, startState,
            targetState, statesToIndices, oldVariablesToIndices, newVariables, 0));
  }

  /*
   * Binarizes this monomial, as one would when calculating a CNF of a grammar. Only works on monomials
   * that have no terminals (i.e. semiring elements) as factors; for all other monomials (i.e. those with
   * terminals as factors), this will return NonCommutativePolynomialBase<LossyFiniteAutomaton>::null().
   *
   * If the monomial consists of exactly one or exactly two variables, the resulting polynomial will not differ from it.
   *
   * The variables and productions that are introduced during that process will be stored in binarizationVariables
   * and binarizationVariablesEquations, respectively.
   */
  NonCommutativePolynomialBase<LossyFiniteAutomaton> binarize(std::map<std::string, VarId> &binarizationVariables,
      NCEquationsBase<LossyFiniteAutomaton> &binarizationVariablesEquations) const {

    // make sure there are no terminal factors in the monomial
    if(idx_.size() != variables_.size()) {
      return NonCommutativePolynomialBase<LossyFiniteAutomaton>::null();
    }

    NonCommutativePolynomialBase<LossyFiniteAutomaton> temp = NonCommutativePolynomialBase<LossyFiniteAutomaton>::one();
    NonCommutativePolynomialBase<LossyFiniteAutomaton> suffix;

    // we have at least two variables; we need to shorten that to exactly two by
    // introducing new variables where needed
    int variablesProcessed = 0;
    int index = 0;
    VarId var;

    // loop over all suffixes to see which need new variables introduced for them
    while(variablesProcessed < get_degree() - 1) {

      // if this is the first iteration, we start with the leftmost variable
      if(variablesProcessed == 0) {
        index = variables_.size() - 1;

        temp = NonCommutativePolynomialBase<LossyFiniteAutomaton>(variables_.at(index));
        suffix = NonCommutativePolynomialBase<LossyFiniteAutomaton>(variables_.at(index));
      } else { /* if this is not the first iteration, that means we already remembered a variable
       * to produce the current suffix of the monomial; multiply the next variable from left,
       * check if some variable already maps to the current suffix, if none does introduce
       * a new variable that produces the product of the now two variables we have, and
       * make this new variable the starting point for the next iteration; also remember
       * the suffix so far so we can check whether maybe we already have a variable that
       * produces that suffix in the next iteration
       */
        temp = NonCommutativePolynomialBase<LossyFiniteAutomaton>(variables_.at(index)) * temp; // multiplication from left is relevant -> noncommutative!
        suffix = NonCommutativePolynomialBase<LossyFiniteAutomaton>(variables_.at(index)) * suffix;

        // check if we already have a variable for the current suffix; if we do, proceed from there;
        // else, introduce a new one
        if(binarizationVariables.find(suffix.string()) != binarizationVariables.end()) {
          temp = NonCommutativePolynomialBase<LossyFiniteAutomaton>(binarizationVariables[suffix.string()]);
        } else {
          var = Var::GetVarId();
          binarizationVariablesEquations.push_back(std::make_pair(var, temp));
          binarizationVariables.insert(std::make_pair(suffix.string(), var));
          temp = NonCommutativePolynomialBase<LossyFiniteAutomaton>(var);
        }
      }

      index--; // we process the monomial from right to left
      variablesProcessed++;
    }


    return NonCommutativePolynomialBase<LossyFiniteAutomaton>(variables_.at(0)) * temp;
  }

  /*
   * Check if this monomial is a linear one that leads to a lower component.
   *
   * Used in the algorithm by Courcelle, see NonCommuatativePolynomial<LossyFiniteAutomaton>::downwardClosureCourcelle.
   */
  void findLowerLinearTerms(std::set<int> &lowerLinearTerms, std::map<VarId, int> &varToComponent, int component) const {
      if((get_degree() != 1) || (varToComponent[variables_[0]] == component)) {
          return;
      }
      lowerLinearTerms.insert(varToComponent[variables_[0]]);
  }


  /*
   * Used when finding the shortest word derivable from a CFG. This method tells us whether
   * the algorithm already knows a terminal string that derives from this monomial.
   */
  bool findLengthOfDerivableStrings(std::map<VarId, unsigned long> &lengthOfShortestTerminal,
      std::map<VarId, NonCommutativeMonomial<LossyFiniteAutomaton>> &productionsForShortestWords,
      VarId &lhsOfProduction) const {

    // if this is a terminal, remember its length if we don't already know it; otherwise, see if we already
    // know the length of some derivable string and update it if necessary
    if(variables_.size() == 0){

      // if this is the shortest terminal the lefthand side of the production can produce, update
      if(getLeadingSR() != LossyFiniteAutomaton::null() && lengthOfShortestTerminal[lhsOfProduction] > getLeadingSR().string().size()) {
        lengthOfShortestTerminal[lhsOfProduction] = getLeadingSR().string().size();
        productionsForShortestWords[lhsOfProduction] = *this;
        return true;
      }

      return false;
    } else {

      unsigned long shortestMonomialLength = 0;

      // check for all variables whether we know how to derive a terminal
      for(VarId var: variables_) {
        if(lengthOfShortestTerminal.at(var) == ULONG_MAX) {
          shortestMonomialLength = ULONG_MAX;
          break;
        }

        shortestMonomialLength += lengthOfShortestTerminal.at(var);
      }

      // if we found a shorter word than previously known for the nonterminal this production
      // belongs to, update
      if(shortestMonomialLength < ULONG_MAX) {

        for(auto terminal: srs_) {
          shortestMonomialLength += terminal.string().size();
        }

        if(shortestMonomialLength < lengthOfShortestTerminal.at(lhsOfProduction)) {
          // update the map with the lengths of shortest terminals
          productionsForShortestWords[lhsOfProduction] = *this;
          lengthOfShortestTerminal[lhsOfProduction] = shortestMonomialLength;
          return true;
        }
      }

      return false;
    }
  }

  /*
   * returns true if some information changed
   */
  bool getDerivableFirstLettes(std::map<VarId, std::set<char>> &varToDerivableFirstLetters, std::set<VarId> &varsWithEpsilon,
      const VarId &lhsOfProduction) const {

    // do not touch null-monomials
    for(int i = 0; i < srs_.size(); i++) {
      if(srs_[i] == LossyFiniteAutomaton::null()) {
        return false;
      }
    }

    bool foundNew = false;

    // first, check whether this monomial can derive epsilon
    if(varsWithEpsilon.find(lhsOfProduction) ==  varsWithEpsilon.end()) {
      bool canProduceEpsilon = true;

      for(int i = 0; i < srs_.size(); i++) {
        if(srs_[i] != LossyFiniteAutomaton::one()) {
          canProduceEpsilon = false;
          break;
        }
      }

      if(canProduceEpsilon) {
        for(int i = 0; i < variables_.size(); i++) {
          if(varsWithEpsilon.find(variables_[i]) == varsWithEpsilon.end()) {
            canProduceEpsilon = false;
            break;
          }
        }
      }

      if(canProduceEpsilon) {
        varsWithEpsilon.insert(lhsOfProduction);
        foundNew = true;
      }
    }

    // now figure out which prefixes we can derive from this monomial given the currently available information
    auto leading_sr = getLeadingSR();
    if(get_degree() == 0) {
      if(leading_sr != LossyFiniteAutomaton::null() && leading_sr != LossyFiniteAutomaton::one()) {
        foundNew |= varToDerivableFirstLetters[lhsOfProduction].insert(leading_sr.string().at(0)).second;
      }
    } else {

      for(int i = 0; i < idx_.size(); i++) {
        if(idx_[i].first == SemiringType && srs_[idx_[i].second] != LossyFiniteAutomaton::one()) {
          foundNew |= varToDerivableFirstLetters[lhsOfProduction].insert(srs_[idx_[i].second].string().at(0)).second;
          break;
        }

        if(idx_[i].first == Variable) {
          for(char letter: varToDerivableFirstLetters[variables_[idx_[i].second]]) {
            foundNew |= varToDerivableFirstLetters[lhsOfProduction].insert(letter).second;
          }

          if(varsWithEpsilon.find(variables_[idx_[i].second]) == varsWithEpsilon.end()) {
            break;
          }
        }
      }
    }

    return foundNew;
  }

  int getLength() const {
    int length = 0;

    for(VarId var: variables_) {
      length++;
    }

    for(LossyFiniteAutomaton sr: srs_) {
      if(sr != LossyFiniteAutomaton::null() && sr != LossyFiniteAutomaton::one()) {
        length += sr.string().length();
      }
    }

    return length;
  }

  /*
   * Derives the shortest derivable terminal for this monomial, assuming the info in productionsForShortestWords
   * is correct.
   */
  LossyFiniteAutomaton shortestDerivableTerminal(std::map<VarId, NonCommutativeMonomial<LossyFiniteAutomaton>> &productionsForShortestWords) const {
    LossyFiniteAutomaton terminal = LossyFiniteAutomaton::one();

    if(get_degree() == 0) {
      terminal = getLeadingSR();
    } else {
      for(auto &indexpair: idx_) {
        if(indexpair.first == SemiringType) {
          terminal *= srs_.at(indexpair.second);
        } else {
          assert(productionsForShortestWords.find(variables_.at(indexpair.second)) != productionsForShortestWords.end());
          terminal *= productionsForShortestWords[variables_.at(indexpair.second)]
                                                  .shortestDerivableTerminal(productionsForShortestWords);
        }
      }
    }

    return terminal;
  }

  /*
   * Checks whether this monomial is productive, depending on the set of variables
   * that are already known to be productive.
   */
  bool isProductive(const std::map<VarId, bool> &productiveVariables) const {

    // if this monomial has no variables, then it is productive since
    // it represents an element of the semiring
    if(get_degree() == 0) {
      return true;
    }

    // check if all variables of this monomial are productive
    for(auto &variable: variables_) {

      // if there is any variable that is not known to be productive in
      // the current monomial, then it isn't productive
      if(!productiveVariables.at(variable)) {
        return false;
      }
    }

    // if all variables are productive, then so is the monomial
    return true;
  }

  /*
   * See NonCommutativePolynomialBase<LossyFiniteAutomaton>::componentIsSquarable(..) for commentary.
   */
  bool componentIsSquarable(std::map<VarId, int> &varToComponent, int component) const {
    int count = 0;

    // if the monomial contains any two variables from the same component,
    // then we can duplicate all variables in the component
    for(auto var: variables_) {
      if(varToComponent[var] == component) {
        count++;

        if(count >= 2) {
          return true;
        }
      }
    }

    return false;
  }

  /*
   * Used in the algorithm by Courcelle, see LossyFiniteAutomaton::downwardClosureCourcelle.
   *
   * Assumes that all productions have been binarized.
   */
  void mapQuadraticLHStoRHS(std::map<int, std::map<int, std::set<int>>> &quadraticLHStoRHS,
      std::map<VarId, int> &varToComponent, int component) const {
    if(get_degree() != 2) {
      return;
    }

    int lhsComponent = varToComponent[variables_[0]];
    int rhsComponent = varToComponent[variables_[1]];

    if((lhsComponent != component) && (rhsComponent != component)) {
      quadraticLHStoRHS[component][lhsComponent].insert(rhsComponent);
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
      std::vector<NCEquationsBase<LossyFiniteAutomaton>> &components,
      std::map<VarId, int> &varToComponent,
      int component) const {

    if(get_degree() != 2) {
      return;
    }

    int lhsComponent = varToComponent[variables_[0]];
    int rhsComponent = varToComponent[variables_[1]];

    if(lhsComponent == component) {
      rhsLowerComponentVariables[component].insert(rhsComponent);
    } else if(rhsComponent == component) {
      lhsLowerComponentVariables[component].insert(lhsComponent);
    }
  }

  /*
   * Used in the algorithm by Courcelle, see LossyFiniteAutomaton::downwardClosureCourcelle.
   *
   * Assumes that all productions have been binarized.
   *
   * WARNING: will break if used with anything but NonCommutativePolynomialBase<LossyFiniteAutomaton>
   */
  void calculateSameComponentLetters(
      std::map<int, std::set<unsigned char>> &lhsSameComponentLetters,
      std::map<int, std::set<unsigned char>> &rhsSameComponentLetters,
      std::vector<NCEquationsBase<LossyFiniteAutomaton>> &components,
      std::map<VarId, int> &varToComponent,
      int component) const {

    if(get_degree() != 1) {
      return;
    }

    int varComponent = varToComponent[variables_[0]];

    if(varComponent == component) {
      std::string lhsString = getLeadingSR().string();
      std::string rhsString = getTrailingSR().string();

      for(int i = 0; i < lhsString.size(); i++) {
        if (lhsString[i] == '[') {
          unsigned char letter, startOfRange, endOfRange;
          i++;
          while(lhsString[i] != ']') {
            if(lhsString[i+1] == '-') {
              startOfRange = lhsString[i];
              endOfRange = lhsString[i+2];
              i += 3;

              for(letter = startOfRange; letter <= endOfRange; letter++) {
                lhsSameComponentLetters[component].insert(letter);
              }
            } else {
              lhsSameComponentLetters[component].insert(lhsString[i]);
              i++;
            }
          }
        } else {
          lhsSameComponentLetters[component].insert(lhsString[i]);
        }
      }

      for(int i = 0; i < rhsString.size(); i++) {
        if(isalnum(rhsString[i])) {
          rhsSameComponentLetters[component].insert(rhsString[i]);
        } else if (rhsString[i] == '[') {
          unsigned char letter, startOfRange, endOfRange;
          i++;
          while(rhsString[i] != ']') {
            if(rhsString[i+1] == '-') {
              startOfRange = rhsString[i];
              endOfRange = rhsString[i+2];
              i += 3;

              for(letter = startOfRange; letter <= endOfRange; letter++) {
                rhsSameComponentLetters[component].insert(letter);
              }
            } else {
              rhsSameComponentLetters[component].insert(rhsString[i]);
              i++;
            }
          }
        }
      }
    }
  }


private:

  /*
   * Used to recursively generate all monomials for an intersection grammar between a CFG
   * and a FiniteAutomaton. See intersectionPolynomial(..) for details.
   *
   * currentState determines the state in which the FA "currently" would be, i.e. if
   * we are currently generating all nonterminals <s_i, A, s_{i+1}> of the new grammar, then
   * currentState would refer us to s_i.
   *
   * See e.g. Nederhof & Satta, "Probabilistic Parsing" (section 3), 2008 for details on the generated nonterminals.
   *
   */
  NonCommutativePolynomialBase<LossyFiniteAutomaton> generateIntersectionMonomials(std::vector<unsigned long> &states,
      std::map<unsigned long, std::map<unsigned char, std::forward_list<unsigned long>>> &transitionTable,
      unsigned long currentState,
      unsigned long targetState,
      std::map<unsigned long, unsigned long> &statesToIndices,
      std::map<VarId, unsigned long> &oldVariablesToIndices,
      std::vector<std::vector<std::vector<VarId>>> &newVariables,
      unsigned long monomialFactorIndex) const {

    NonCommutativePolynomialBase<LossyFiniteAutomaton> result = NonCommutativePolynomialBase<LossyFiniteAutomaton>::null();

    // if the current position in the monomial is a semiring element, we see where the transitions using the element
    // starting at currentStartIndex can take us and recurse from there
    // otherwise, we replace the variable with all replacement variables and recursive from there
    if(idx_[monomialFactorIndex].first == SemiringType) { // we have a semiring element

      // check if we hit the last factor of the monomial; we stop the recursion here
      if(monomialFactorIndex == idx_.size() - 1) {

        // if the last element in the monomial is epsilon, then the only way the generated monomial is valid
        // is if the target state is equal to the one we are currently at
        if(srs_[idx_[monomialFactorIndex].second] == LossyFiniteAutomaton::one()) {
          if(currentState == targetState) {
            result = NonCommutativePolynomialBase<LossyFiniteAutomaton>::one();
          }
        } else { // otherwise, see if we can reach a target state

          // get a string representation of the semiring element, see where that string takes us in the FA
          // -----------------------------------------------------
          // -----------------------------------------------------
          // -----------------------------------------------------
          // note that due to this bit here, we can only work with grammars where we have elements of SIGMA*
          // as semiring factors; if we want to use regular expressions over SIGMA here, the following
          // needs to be adjusted
          // -----------------------------------------------------
          // -----------------------------------------------------
          // -----------------------------------------------------
          std::string transitionsToDo = srs_[idx_[monomialFactorIndex].second].string();

          std::set<unsigned long> reachable, temp;
          reachable.insert(currentState);

          for(int i = 0; i < transitionsToDo.size(); i++) {

            for(auto &intermediateState: reachable) {
              for(unsigned long &nextIntermediateState: transitionTable[intermediateState][transitionsToDo[i]]) {
                temp.insert(nextIntermediateState);
              }
            }

            reachable.swap(temp);
            temp.clear();
          }

          // if we can reach the intended state with the given transitions, we allow the semiring element
          if(reachable.count(targetState) != 0) {
            result = NonCommutativePolynomialBase<LossyFiniteAutomaton>(srs_[idx_[monomialFactorIndex].second]);
          }
        }
      } else { // if this is not the last factor of the monomial, continue the recursion

        // if the semiring element is 1, then the FA cannot make a move since it doesn't have epsilon transitions,
        // so we advance in the monomial but leave the currentState untouched
        if(srs_[idx_[monomialFactorIndex].second] == LossyFiniteAutomaton::one()) {
          result = generateIntersectionMonomials(states, transitionTable, currentState, targetState,
              statesToIndices, oldVariablesToIndices, newVariables, (monomialFactorIndex + 1));
        } else { // if the semiring element is not the 1 element, get a string representation of it
          // and see where that string takes us in the FAstateTransitions

          std::string transitionsToDo = srs_[idx_[monomialFactorIndex].second].string();
          std::set<unsigned long> reachable, temp;
          reachable.insert(currentState);

          for(int i = 0; i < transitionsToDo.size(); i++) {

            for(auto &intermediateState: reachable) {
              for(unsigned long &nextIntermediateState: transitionTable[intermediateState][transitionsToDo[i]]) {
                temp.insert(nextIntermediateState);
              }
            }

            reachable.swap(temp);
            temp.clear();
          }

          // we proceed with generating the monomial starting at the next factor and all reachable states
          NonCommutativePolynomialBase<LossyFiniteAutomaton> suffixPolynomial = NonCommutativePolynomialBase<LossyFiniteAutomaton>::null();
          for(auto &nextState: reachable) {
            suffixPolynomial += generateIntersectionMonomials(states, transitionTable, nextState, targetState,
                statesToIndices, oldVariablesToIndices, newVariables, (monomialFactorIndex + 1));
          }

          result = NonCommutativePolynomialBase<LossyFiniteAutomaton>(srs_[idx_[monomialFactorIndex].second]) * suffixPolynomial;
        }
      }
    } else { // we have a variable at the current index

      // check if we hit the last factor of the monomial; we stop the recursion here
      if(monomialFactorIndex == idx_.size() - 1) {
        result = NonCommutativePolynomialBase<LossyFiniteAutomaton>(newVariables[statesToIndices[currentState]][statesToIndices[targetState]]
                                                                                           [oldVariablesToIndices.at(variables_[idx_[monomialFactorIndex].second])]);
      } else { // if this is not the last factor of the monomial, continue the recursion
        for(auto &nextState: states) {

          // we want to replace some nonterminal X by <currentState, X, nextState> for all nextState;
          // after that we need to append all possible monomials generated from the suffix of this monomial
          // that starts one symbol after X and uses nextState
          result += NonCommutativePolynomialBase<LossyFiniteAutomaton>(newVariables[statesToIndices[currentState]][statesToIndices[nextState]]
                                                                                              [oldVariablesToIndices.at(variables_[idx_[monomialFactorIndex].second])]) *
                                                                                                  generateIntersectionMonomials(states, transitionTable, nextState, targetState,
                                                                                                      statesToIndices, oldVariablesToIndices, newVariables, (monomialFactorIndex + 1));
        }
      }
    }

    return result;
  }


};




#endif /* LOSSY_NON_COMMUTATIVE_MONOMIAL_H_ */
