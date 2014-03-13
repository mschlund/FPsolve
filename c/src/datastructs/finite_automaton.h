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
        automaton = fa_make_basic(0);
        assert(automaton!=NULL);
    }

    /*
     * Builds a finite automaton from a POSIX regular expression.
     * Use | for choice, () for grouping, * for Kleene star; don't use spaces,
     * also concatenation does not require a symbol.
     */
    FiniteAutomaton(std::string regularExpression) {
        fa_compile(regularExpression.c_str(), regularExpression.size(), &automaton);
        std::cout << "automaton language:\t" + string() << std::endl;
        assert(automaton!=NULL);
    }

    // every other way I tried to build an automaton for epsilon messed something up,
    // so we're doing it the ugly way
    static FiniteAutomaton epsilon() {
        FiniteAutomaton returnMe = FiniteAutomaton();
        returnMe.automaton = fa_make_basic(1);
        return returnMe;
    }

    bool empty() {
        return fa_equals(automaton, fa_make_basic(0)) == 1;
    }

    FiniteAutomaton minimize() {
        struct fa* automatonClone = automaton;
        fa_minimize(automatonClone);

        return FiniteAutomaton(automatonClone);
    }

    FiniteAutomaton complement() {
        return FiniteAutomaton(fa_complement(automaton));
    }

    FiniteAutomaton intersectionWith(FiniteAutomaton other) {
            return FiniteAutomaton(fa_intersect(automaton, other.automaton));
    }

    FiniteAutomaton unionWith(FiniteAutomaton other) {
            return FiniteAutomaton(fa_union(automaton, other.automaton));
    }

    FiniteAutomaton concatenate(FiniteAutomaton other) {
        return FiniteAutomaton(fa_concat(automaton, other.automaton));
    }

    FiniteAutomaton kleeneStar() const {
        return FiniteAutomaton(fa_iter(automaton, 0, -1));
    }

    bool disjoint(FiniteAutomaton other) {
        return intersectionWith(other).empty();
    }

    bool containedIn(FiniteAutomaton super) {
        return (fa_contains(automaton, super.automaton) == 1);
    }

    bool contains(FiniteAutomaton sub) {
        return (fa_contains(sub.automaton, automaton) == 1);
    }

    bool equals(FiniteAutomaton other) {
        return contains(other) && containedIn(other);
    }

    std::string string() const {
        char *regularExpression;
        size_t size;
        fa_as_regexp(automaton, &regularExpression, &size);
        std::string str(regularExpression);
        return str;
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
        (VarId &S, const VarId &oldS,
                std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &oldGrammar) {

        // we don't want to start out with useless stuff because the intersection grammar
        // will blow up anyway
        std::queue<VarId> worklist;
        worklist.push(oldS);
        oldGrammar = NonCommutativePolynomial<SR>::cleanSystem(oldGrammar, worklist);

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

        // the number of states in the FA since forward_list doesn't have that information available
        long numberOfStates = 0;

        // hashes of all final states of the FA
        std::set<unsigned long> finalStates;

        // initial state of the FA
        unsigned long initialState = 0;
        extractTransitionTable(transitionTable, states, numberOfStates, finalStates, initialState);

        // we will need a three-dimensional lookup table for the new variables so we don't mess up
        // the assignment of polynomials to nonterminals while constructing the grammar
        std::vector<std::vector<std::vector<VarId>>> newVariables;

        // the lookup table will represent states_FA x states_FA x nonterminals_grammar;
        // initialize its dimensions accordingly
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
                    newVariables[i][j][k] = Var::GetVarId();
                }
            }
        }

        // build a map from variables to indices, where the index of a variable is the line in oldGrammar
        // that contains the productions of the variable
        std::map<VarId, unsigned long> oldVariablesToIndices;
        unsigned long index = 0;
        for(auto equation: oldGrammar) {
            oldVariablesToIndices.insert(std::make_pair(equation.first, index));
            index++;
        }

        /*
         * now we the productions of the new grammar
         */

        // for the new grammar
        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> resultGrammar;

        // the indices are intended this way; we first choose some nonterminal of the grammar
        // and then generate all productions derived from that nonterminal; this mirrors the way
        // in which the algorithm is described in the paper referred to above
        NonCommutativePolynomial<SR> poly;
        for(unsigned long k = 0; k < oldGrammar.size(); k++) {
            for(unsigned long i = 0; i < numberOfStates; i++) {
                for(unsigned long j = 0; j < numberOfStates; j++) {
                    poly = oldGrammar[k].second.intersectionPolynomial
                            (states, transitionTable, i, j, /*oldGrammar,*/ oldVariablesToIndices, newVariables);
                    resultGrammar.push_back(std::make_pair(newVariables[i][j][k], poly));
                }
            }
        }
    }

    //static const FiniteAutomaton UNIVERSE_AUTOMATON;
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
            long &numberOfStates,
            std::set<unsigned long> &finalStates,
            unsigned long &initialState) {

        initialState = automaton->initial->hash;

        // iterate through all states, build the transition table state for state
        struct trans *trans;
        for (auto s = automaton->initial->next; s != automaton->initial; s = s->next){

            // remember each state by its hash, also count them
            states.push_back(s->hash);
            numberOfStates++;

            // remember the state as being final if it is
            if(s->accept) {
                finalStates.insert(s->hash);
            }

            // we build a transition table for each state, i.e. what would be a line in the
            // transition table of the whole FA
            std::map<unsigned char, std::forward_list<unsigned long>> stateTable;
            transitionTable.insert(std::pair<unsigned long,
                    std::map<unsigned char, std::forward_list<unsigned long>>>(s->hash, stateTable));

            // trans is an array of transitions
            trans = s->trans;

            // iterate over all transitions
            for(int i = 0; i < s->tused; i++) {

                // iterate over all symbols of this transition
                for(unsigned char c = (trans+i)->min; c <= (trans+i)->max; c++) {

                    // the map[] operator will create a new set of target states if none yet exists
                    stateTable[c].push_front((trans+i)->to->hash);
                }
            }
        }
    }

    struct fa* automaton;
};
