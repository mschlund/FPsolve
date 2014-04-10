#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <forward_list>
extern "C" {
    #include <fa.h>
}

#include "regular_language.h"
#include "../polynomials/non_commutative_polynomial.h"
#include "../datastructs/var.h"

class FiniteAutomaton : public RegularLanguage<FiniteAutomaton> {
public:

    /*
     * FA that accepts the empty language, i.e. no strings, not even the empty string.
     */
    FiniteAutomaton() {
        automaton = fa_make_basic(fa_basic::FA_EMPTY);
    }

    /*
     * Builds a finite automaton from a POSIX regular expression.
     * Use | for choice, () for grouping, * for Kleene star; don't use spaces,
     * also concatenation does not require a symbol.
     */
    FiniteAutomaton(std::string regularExpression) {
        fa_compile(regularExpression.c_str(), regularExpression.size(), &automaton);
        assert(automaton!=NULL);
    }

    /*
     * Returns an automaton that accepts exactly the empty word.
     */
    static FiniteAutomaton epsilon() {
        return FiniteAutomaton(fa_make_basic(fa_basic::FA_EPSILON));
    }

    /*
     * True iff the language of this FA is empty.
     */
    bool empty() const {
        return fa_equals(automaton, fa_make_basic(fa_basic::FA_EMPTY)) == 1;
    }

    /*
     * Minimizes the automaton of this FA. The returned automaton is the one on which the method was invoked;
     * there is no copy involved. This is due to restrictions of the library we are using.
     */
    FiniteAutomaton minimize() const {
        fa_minimize(automaton);
        return *this;
    }

    /*
     * Complements the language of this automaton.
     */
    FiniteAutomaton complement() const {
        return FiniteAutomaton(fa_complement(automaton));
    }

    FiniteAutomaton minus(FiniteAutomaton other) const {
        return FiniteAutomaton(fa_minus(automaton, other.automaton));
    }

    /*
     * Automaton for the intersection between the languages of two automata.
     */
    FiniteAutomaton intersectionWith(FiniteAutomaton other) const {
            return FiniteAutomaton(fa_intersect(automaton, other.automaton));
    }

    /*
     * Automaton for the union of the languages of two automata.
     */
    FiniteAutomaton unionWith(FiniteAutomaton other) const {
            return FiniteAutomaton(fa_union(automaton, other.automaton));
    }

    /*
     * Automaton for the concatenation of the languages of two automata.
     *
     * The concatenation order is this.append(other).
     */
    FiniteAutomaton concatenate(FiniteAutomaton other) const {
        return FiniteAutomaton(fa_concat(automaton, other.automaton));
    }

    /*
     * Automaton for the Kleene iteration of this automaton.
     */
    FiniteAutomaton kleeneStar() const {
        return FiniteAutomaton(fa_iter(automaton, 0, -1));
    }

    /*
     * True iff L(this) and L(other) are disjoint, i.e. there is no word that is in both languages.
     */
    bool disjoint(FiniteAutomaton other) const {
        return intersectionWith(other).empty();
    }

    /**
     * True iff L(this) is a superset of L(sub).
     */
    bool contains(FiniteAutomaton sub) const {
        return (fa_contains(sub.automaton, automaton) == 1);
    }

    /*
     * True iff L(this) == L(other).
     */
    bool equals(FiniteAutomaton other) const {
        return (fa_equals(automaton, other.automaton) == 1);
    }

    /*
     * Clones this automaton, i.e. creates another automaton with the same set of states, transitions and final states
     * as this one.
     */
    FiniteAutomaton clone() const {

        /*
         * Currently this workaround is required since we don't yet have access to the function fa_clone(..)
         * that's available in libfa's sourcecode.
         */
        return concatenate(epsilon());
    }

    /*
     * A string representation of the language of this automaton. Since we lack a symbol for it,
     * the empty language will be represented by the string "__EMPTYLANGUAGE__".
     */
    std::string string() const {

        // first handle empty and epsilon separately because
        // the library doesn't handle those the way we would want it to
        if(empty()) {
            return std::string("__EMPTYLANGUAGE__");
        }

        if(equals(epsilon())) {
            return std::string("()");
        }

        char *regularExpression;
        size_t size;
        fa_as_regexp(automaton, &regularExpression, &size);
        return std::string(regularExpression);
    }

    int size() const {
        int states = 0;

        for (auto s = automaton->initial; s != NULL; s = s->next){
            states++;
        }

        return states;
    }

    /*
     * Calculates the intersection of the CFG "equations" with start symbol "oldS" and the
     * regular language represented by this finite automaton. The result is a CFG given by
     * its set of productions (return value of this method) and its start symbol "S".
     *
     * The algorithm is that of Bar-Hillel et al., "On formal properties of simple phrase structure
     * grammars.", 1961; a neat description of it can be found in Nederhof & Satta, "Probabilistic Parsing", 2008
     * which is openly accessible on Nederhof's website as of the writing of this code (spring 2014).
     *
     * libfa doesn't model epsilon transitions (see add_epsilon_trans(struct fa, struct fa) in fa.c of the
     * source code of libfa for details), instead replacing them by transitions on nonempty symbols and
     * adjusting the set of final states; thus the simplifying assumption of "Probabilistc Parsing" that
     * there are no epsilon transitions in the FA holds, and we only have to deal with the fact that the FA
     * may have multiple final states.
     *
     * WARNING: do not allow productions of the form "X -> empty set" in your oldGrammar, or this may break.
     */
    template<typename SR>
    std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> intersectionWithCFG
        (VarId &newS, const VarId &oldS,
                std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &oldGrammar) {

        // for the new grammar
        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> resultGrammar;

        // we don't want to start out with useless stuff because the intersection grammar
        // will blow up anyway
        std::queue<VarId> worklist;
        worklist.push(oldS);
        oldGrammar = NonCommutativePolynomial<SR>::cleanSystem(oldGrammar, worklist);

        // if the grammar doesn't produce anything, there is no need to do anything else;
        // return an empty grammar
        // the same applies if the automaton represents the empty language
        if(oldGrammar.size() == 0 || empty()) {
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

//        std::cout << "lets see if we can get that FA info" << std::endl;
        extractTransitionTable(transitionTable, states, finalStates, initialState);

//        std::cout << "number of final states: " << finalStates.size() << std::endl;
//        std::cout << "number of states: " << states.size() << std::endl;
//        for(auto &transitions: transitionTable){
//            std::cout << "transitions for state:\t" << transitions.first << std::endl;
//
//            for(auto &trans: transitions.second) {
//                for(auto &target: trans.second) {
//                    std::cout << "delta(" << transitions.first << "," << trans.first << ") = " << target << std::endl;
//                }
//
//            }
//        }

//        if(transitionTable.size() == 0) {
//            std::cout << "something is fishy; transition table is empty" << std::endl;
//        }
//
//        std::cout << "extraction of FA info does not cause segfault" << std::endl;

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
                newVariables[i][j].resize(oldGrammar.size());
            }
        }

        // generate the variables of the new grammar
        for(unsigned long i = 0; i < numberOfStates; i++) {
            for(unsigned long j = 0; j < numberOfStates; j++) {
                for(unsigned long k = 0; k < oldGrammar.size(); k++) {
                    std::stringstream ss;
                    ss << "<" << i << "," << Var::GetVar(oldGrammar[k].first).string() << "," << j << ">";
                    newVariables[i][j][k] = Var::GetVarId(ss.str());
//                    newVariables[i][j][k] = Var::GetVarId();
//                    std::cout << "new var introduced: " << Var::GetVar(newVariables[i][j][k]).string() << std::endl;
                }
            }
        }

//        std::cout << "neither does generating new variables" << std::endl;

        // build a map from variables to indices, where the index of a variable is the line in oldGrammar
        // that contains the productions of the variable
        std::map<VarId, unsigned long> oldVariablesToIndices;
        for(unsigned long i = 0; i < oldGrammar.size(); i++) {
            oldVariablesToIndices.insert(std::make_pair(oldGrammar[i].first, i));
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


//        std::cout << "the problem must lie in the actual intersection algorithm" << std::endl;

        // the indices k, i, j are intended this way; we first choose some nonterminal of the grammar
        // and then generate all productions derived from that nonterminal; this mirrors the way
        // in which the algorithm is described in the paper referred to above
        NonCommutativePolynomial<SR> poly;
        for(unsigned long k = 0; k < oldGrammar.size(); k++) {
            for(unsigned long i = 0; i < numberOfStates; i++) {
                for(unsigned long j = 0; j < numberOfStates; j++) {
//                    std::cout << "k, i, j: " << k << ", " << i << ", " << j << std::endl;
                    poly = oldGrammar[k].second.intersectionPolynomial
                            (states, transitionTable, states[i], states[j], statesToIndices, oldVariablesToIndices, newVariables);
                    resultGrammar.push_back(std::make_pair(newVariables[i][j][k], poly));
                }
            }
        }

//        for(auto state: states) {
//            std::cout << "state: " << state << std::endl;
//            std::cout << "statesToIndices[state]: " << statesToIndices[state] << std::endl;
//        }

        // finally, use a new variable as initial variable and generate the productions; they are of the form
        // newS -> <q_0, oldS, q_f> where q_0 is the initial state of the FA and q_f is some final state of the FA
        NonCommutativePolynomial<SR> startPolynomial = NonCommutativePolynomial<SR>::null();
//        std::cout << "start polynomial initially: " << startPolynomial.string() << std::endl;
//        std::cout << "initial state: " << initialState << std::endl;
        for(auto target: finalStates) {
//            std::cout << "current final state: " << target << std::endl;
//            std::cout << "index mapping: " << statesToIndices[target] << std::endl;
            startPolynomial += newVariables[statesToIndices[initialState]][statesToIndices[target]][oldVariablesToIndices[oldS]];
//            std::cout << "start polynomial: " << startPolynomial.string() << std::endl;
        }
        newS = Var::GetVarId(std::string("S"));
        resultGrammar.push_back(std::make_pair(newS, startPolynomial));

        auto startProduction = resultGrammar[resultGrammar.size() - 1];
        resultGrammar[resultGrammar.size() - 1] = resultGrammar[0];
        resultGrammar[0] = startProduction;
        return resultGrammar;
    }

private:
    FiniteAutomaton(struct fa *FA) {
        automaton = FA;
    }

    /*
     * Builds the transition table of this FA. The table will map the hash of a state to a map that maps transition
     * symbols to hashes of target states. Since we are having general FAs, the transition relation need not be a
     * function, and so for each pair (state, symbol), we may have more than one target state; for that reason we're
     * using a list of target states here.
     */
    void extractTransitionTable(std::map<unsigned long, std::map<unsigned char,
                std::forward_list<unsigned long>>> &transitionTable,
            std::vector<unsigned long> &states,
            std::set<unsigned long> &finalStates,
            unsigned long &initialState) {

        initialState = automaton->initial->hash;
//        std::cout << "the initial state is:\t" << initialState << std::endl;

        // iterate through all states, build the transition table state for state
        struct trans *trans;
        for (auto s = automaton->initial; s != NULL; s = s->next){

            // remember each state by its hash, also count them
            states.push_back(s->hash);
//            std::cout << "some state:\t" << s->hash << std::endl;

            // remember the state as being final if it is
            if(s->accept) {
                finalStates.insert(s->hash);
//                std::cout << "found a final state:\t" << s->hash << std::endl;
            }

            // we build a transition table for each state, i.e. what would be a line in the
            // transition table of the whole FA
            std::map<unsigned char, std::forward_list<unsigned long>> stateTable;

            // trans is an array of transitions
            trans = s->trans;

            // iterate over all transitions
            for(int i = 0; i < s->tused; i++) {

                // iterate over all symbols of this transition
                for(unsigned char c = (trans+i)->min; c <= (trans+i)->max; c++) {

                    // the map[] operator will create a new empty list of target states if none yet exists
                    stateTable[c].push_front((trans+i)->to->hash);
//                    std::cout << "a transition:\tdelta(" << s->hash << "," << c << ") = " << (trans+i)->to->hash << std::endl;
                }
            }

            transitionTable.insert(std::make_pair(s->hash, stateTable));
//            std::cout << "extracting for this state is done" << std::endl;
        }

//        std::cout << "done extracting" << std::endl;
    }

    struct fa* automaton;
};
