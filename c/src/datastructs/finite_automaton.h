#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <forward_list>
#include <queue>

extern "C" {
    #include <fa.h>
}
#include <cassert>

#include "var.h"
#include "hash.h"

#include "regular_language.h"
#include "../polynomials/non_commutative_polynomial.h"


class FiniteAutomaton : public RegularLanguage<FiniteAutomaton> {
public:

    /*
     * FA that accepts the empty language, i.e. no strings, not even the empty string.
     */
    FiniteAutomaton() {
        automaton = fa_make_basic(fa_basic::FA_EMPTY);
        epsilon_closed = 0;
    }

    /*
     * Builds a finite automaton from a POSIX regular expression.
     * Use | for choice, () for grouping, * for Kleene star; don't use spaces,
     * also concatenation does not require a symbol.
     */
    FiniteAutomaton(std::string regularExpression) {
        fa_compile(regularExpression.c_str(), regularExpression.size(), &automaton);
        epsilon_closed = 0;
        assert(automaton!=NULL);
    }

    /*
     * Returns an automaton that accepts exactly the empty word.
     */
    static FiniteAutomaton epsilon() {
        return FiniteAutomaton(fa_make_basic(fa_basic::FA_EPSILON));
    }


    /*
     * Used in the calculation of the downward closure according to Courcelle,
     * see lossy-semiring.h#downwardClosureCourcelle.
     *
     * If you change the order in which the automata are wired up, the algorithm WILL break!
     * The problem is that libfa doesn't model epsilon-transitions explicitly
     */
    static FiniteAutomaton courcelle_construct(FiniteAutomaton lefthandAlphabetStar,
            FiniteAutomaton righthandAlphabetStar,
            std::map<int, std::set<int>> &LHStoRHS,
            FiniteAutomaton closureOfConstantMonomials,
            std::set<int> &lowerLinearTerms,
            std::map<int, FiniteAutomaton> &componentToClosureAutomaton,
            int comp) {

        /*
         * build two sets that hold the components required for the bipartite construction
         */
        std::set<int> lhComponents;
        std::set<int> rhComponents;

        for(auto it = LHStoRHS.begin(); it != LHStoRHS.end(); it++) {
            lhComponents.insert(it->first);

            for(auto rhc: it->second) {
                rhComponents.insert(rhc);
            }
        }

        /*
         * these maps will store the information regarding which states need to be linked with
         * epsilon transitions before we calculate the epsilon closure
         */
        std::set<state*> lhComponentsInitialState;
        std::map<int, state*> rhComponentInitialState;
        std::map<int, std::set<state*>> lhComponentFinalStates;
        std::set<state*> rhComponentsFinalStates;

        /*
         * once we're done, this will contain the automaton for the bipartite construction
         */
        struct fa* finalAutomaton = lefthandAlphabetStar.automaton;
        struct state* source = finalAutomaton->initial;
        assert (source->next == NULL); // otherwise, the construction doesn't work as intended

        // first, clone all the automata we need for the left hand side; remember the
        // initial and final states of the clones
        for(auto comp: lhComponents) {

            // clone
            struct fa* clone = fa_clone(componentToClosureAutomaton[comp].automaton);

            // remember initial state of clone
            lhComponentsInitialState.insert(clone->initial);

            // remember final states of clone
            for (auto s = clone->initial; s != NULL; s = s->next){
                if(s->accept == 1) {
                    lhComponentFinalStates[comp].insert(s);
                }
            }

            fa_merge(finalAutomaton, &clone);

            // erase the component from the set of linear terms that need taking care of; if a component
            // appears in the bipartite construction, the epsilon closure will make it possible to
            // also produce all terminals derivable from each linear component appearing on the left hand side
            // or right hand side of the construction
            lowerLinearTerms.erase(comp);
        }

        // now, clone all the automata we need for the right hand side; remember the
        // initial and final states of the clones
        for(auto comp: rhComponents) {

            // clone
            struct fa* clone = fa_clone(componentToClosureAutomaton[comp].automaton);

            // remember initial state of clone
            rhComponentInitialState[comp] = clone->initial;

            // remember final states of clone
            for (auto s = clone->initial; s != NULL; s = s->next){
                if(s->accept == 1) {
                    rhComponentsFinalStates.insert(s);
                }
            }

            fa_merge(finalAutomaton, &clone);

            // as above
            lowerLinearTerms.erase(comp);
        }

        // merge the right hand alphabet automaton into the construction
        struct state* target = righthandAlphabetStar.automaton->initial;
        assert (target->next == NULL); // otherwise, the construction doesn't work as intended

        fa_merge(finalAutomaton, &(righthandAlphabetStar.automaton));

        // add the epsilon transitions from the final states of the right hand automata to the
        // state of the right hand alphabet automaton
        for(auto rhFinal: rhComponentsFinalStates) {
            add_epsilon_trans(rhFinal, target);
        }

        /*
         *  add the meat of the construction: the actual bipartite construction
         */

        // iterate over all left hand automata
        for(auto lhIt = lhComponentFinalStates.begin(); lhIt != lhComponentFinalStates.end(); lhIt++) {

            // for every left hand automaton, find all right hand automata it is connected with
            for(auto rhc: LHStoRHS[lhIt->first]){

                // for every final state of the left hand automaton, add and epsilon transition to the initial state
                // of the right hand automataon
                for(auto lhFinalState: lhComponentFinalStates[lhIt->first]) {
                    add_epsilon_trans(lhFinalState, rhComponentInitialState[rhc]);
                }
            }
        }

        // add the transitions from left hand alphabet automaton to the initial states of the left
        // hand automata of the bipartite construction
        for(auto lhInitial: lhComponentsInitialState) {
            add_epsilon_trans(source, lhInitial);
        }

        // add the missing linear terms
        for(int linearComp: lowerLinearTerms) {
            struct fa* clone = fa_clone(componentToClosureAutomaton[linearComp].automaton);
            struct state* linearInitial = clone->initial;

            std::set<state*> finalStates;

            // remember final states of clone
            for (auto s = clone->initial; s != NULL; s = s->next){
                if(s->accept == 1) {
                    finalStates.insert(s);
                }
            }

            // add the linear term automaton to our construction
            fa_merge(finalAutomaton, &clone);


            // wire it up with the required transitions
            for(auto linearFinal: finalStates) {
                add_epsilon_trans(linearFinal, target);
            }
            add_epsilon_trans(source, linearInitial);
        }

        // add the constant terms
        if(!closureOfConstantMonomials.empty()) {
            struct state* constantInitial = closureOfConstantMonomials.automaton->initial;
            std::set<state*> finalsOfConstant;
            for(auto s = closureOfConstantMonomials.automaton->initial; s != NULL; s = s->next) {
                if(s->accept == 1) {
                    finalsOfConstant.insert(s);
                }
            }

            fa_merge(finalAutomaton, &(closureOfConstantMonomials.automaton));

            for(auto s: finalsOfConstant) {
                add_epsilon_trans(s, target);
            }
            add_epsilon_trans(source, constantInitial);
        }

        finalAutomaton->minimal = 0;
        finalAutomaton->deterministic = 0;

        return FiniteAutomaton(finalAutomaton);
    }

    /*
     * True iff the language of this FA is empty.
     */
    bool empty() const {
        return fa_equals(automaton, fa_make_basic(fa_basic::FA_EMPTY)) == 1;
    }

    /*
     * Minimizes the automaton of this FA; involves determinization. The returned FiniteAutomaton
     * is the one on which the method was invoked; it's just the internal representation of the
     * automaton that is being modified. This is due to restrictions of the library we are using.
     */
    FiniteAutomaton minimize() const {
        fa_minimize(automaton);
        return *this;
    }

    /*
     * Complements the language of this automaton; involves determinization.
     */
    FiniteAutomaton complement() const {
        return FiniteAutomaton(fa_complement(automaton));
    }

    /*
     * Constructs an automaton for L(this)\L(other); involves determinization.
     */
    FiniteAutomaton minus(FiniteAutomaton other) const {
        return FiniteAutomaton(fa_minus(automaton, other.automaton));
    }

    /*
     * Automaton for the intersection between the languages of two automata.
     */
    FiniteAutomaton intersectionWith(FiniteAutomaton other) const {
        FiniteAutomaton intersection = FiniteAutomaton(fa_intersect(automaton, other.automaton));

        if((other.epsilon_closed == 1) && (epsilon_closed ==1)) {
            intersection.epsilon_closed = 1;
        }

        return intersection;
    }

    /*
     * Automaton for the union of the languages of two automata.
     */
    FiniteAutomaton unionWith(FiniteAutomaton other) const {
        FiniteAutomaton unionAutomaton = FiniteAutomaton(fa_union(automaton, other.automaton));

        if((other.epsilon_closed == 1) && (epsilon_closed ==1)) {
            unionAutomaton.epsilon_closed = 1;
        }

        return unionAutomaton;
    }

    /*
     * Automaton for the concatenation of the languages of two automata.
     *
     * The concatenation order is this.append(other).
     */
    FiniteAutomaton concatenate(FiniteAutomaton other) const {
        FiniteAutomaton concatenation = FiniteAutomaton(fa_concat(automaton, other.automaton));

        if((other.epsilon_closed == 1) && (epsilon_closed ==1)) {
            concatenation.epsilon_closed = 1;
        }

        return concatenation;
    }

    /*
     * Automaton for the Kleene iteration of this automaton.
     */
    FiniteAutomaton kleeneStar() const {
        FiniteAutomaton star = FiniteAutomaton(fa_iter(automaton, 0, -1));

        if(epsilon_closed ==1) {
            star.epsilon_closed = 1;
        }

        return star;
    }

    /*
     * Calculates the epsilon closure of this automaton, i.e. for every state s
     * and every state t that is reachable from s, adds the transition d(s,epsilon) = t.
     */
    FiniteAutomaton epsilonClosure() const {

        if(empty() || epsilon_closed == 1) {
            return *this;
        }

        std::set<unsigned long> foundStates;
        std::vector<state*> bfsQueue;

        struct fa* automaton = fa_clone(this->automaton); // we dont want to alter the original automaton,
                                                    // we want a new one that is lossy
        struct trans *trans;
        struct state *next;
        for (auto s = automaton->initial; s != NULL; s = s->next){
            foundStates.clear();
            bfsQueue.clear();

            foundStates.insert(s->hash);
            bfsQueue.push_back(s);
            std::set<state*> targets;

            while(!bfsQueue.empty()) {
                next = bfsQueue.back();
                bfsQueue.pop_back();
                trans = next->trans;

                // remove this line to earn yourself hours of debugging:
                int iterations = next->tused;

                // iterate over all transitions
                for(int i = 0; i < iterations; i++) {

                    if(foundStates.count((trans+i)->to->hash) == 0) {
                        foundStates.insert((trans+i)->to->hash);
                        bfsQueue.push_back((trans+i)->to);

                        //remember the new target state
                        targets.insert((trans+i)->to);
                    }
                }
            }

            // epsilon-close
            for(auto target: targets) {
                add_epsilon_trans(s, target);
            }
        }

        FiniteAutomaton closure = FiniteAutomaton(automaton);
        closure.epsilon_closed = 1;

        // libfa doesn't set the following flags since it doesn't expect the user to
        // use any of the low level FA fiddling, so we need to do this here
        closure.automaton->deterministic = 0;
        closure.automaton->minimal = 0;

        return closure;
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

    void write_dot_file(FILE * file) const {
        fa_dot(file, automaton);
    }

    int size() const {
        int states = 0;

        for (auto s = automaton->initial; s != NULL; s = s->next){
            states++;
        }

        return states;
    }

    std::set<unsigned char> alphabet() const{
        std::set<unsigned char> alphabet;

        struct trans *trans;

        // iterate through all states
        for (auto s = automaton->initial; s != NULL; s = s->next){

            // trans is an array of transitions
            trans = s->trans;

            // iterate over all transitions
            for(int i = 0; i < s->tused; i++) {

                // iterate over all symbols of this transition
                for(unsigned char c = (trans+i)->min; c <= (trans+i)->max; c++) {
                    alphabet.insert(c);
                }
            }
        }

        return alphabet;
    }

    std::map<int, std::set<std::string>> prefixesToMaxLength(int maxLength, std::set<char> &derivableFirstLetters) const {
        std::map<int, std::map<unsigned long, std::set<std::string>>> prefixesViaStates;
        std::map<unsigned long, struct state*> hashToState;
        std::map<int, std::set<std::string>> prefixes;

        if(maxLength <= 0) {
            return prefixes;
        }

        /* find the prefixes of length 1; afterwards, to find prefixes of
         * length i > 1, use the prefixes of length i-1 > 0
         */

        struct trans *trans;

        auto s = automaton->initial;

        // trans is an array of transitions
        trans = s->trans;

        // iterate over all transitions
        for(int i = 0; i < s->tused; i++) {

            // iterate over all symbols of this transition
            for(unsigned char c = (trans+i)->min; c <= (trans+i)->max; c++) {
                std::string prefix (1, c);
                hashToState[(trans+i)->to->hash] = (trans+i)->to;
                if(derivableFirstLetters.find(prefix.at(0)) != derivableFirstLetters.end()) {
                    prefixesViaStates[1][(trans+i)->to->hash].insert(prefix);
                    prefixes[1].insert(prefix);
                }
            }
        }

        // now use prefixes of length i-1 > 0 to build those of length i > 1
        for(int length = 1; length < maxLength; length++) {

            // mapping: state to prefixes that end there
            for(auto &entry: prefixesViaStates[length]) {

                auto sourceState = hashToState[entry.first];
                trans = sourceState->trans;

                // iterate over prefixes of length i-1
                for(auto &prePrefix: entry.second) {

                    // iterate over transitions
                    for(int j = 0; j < sourceState->tused; j++) {


                        // iterate over chars encoded in this transition
                        for(unsigned char c = (trans+j)->min; c <= (trans+j)->max; c++) {

                            // append the character to the prefix of length i-1 to make one of length i
                            std::string prefix = std::string(prePrefix).append(1,c);

                            // store the results
                            hashToState[(trans+j)->to->hash] = (trans+j)->to;
                            prefixesViaStates[length+1][(trans+j)->to->hash].insert(prefix);
                            prefixes[length+1].insert(prefix);
                        }
                    }
                }
            }
        }

        return prefixes;
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
            unsigned long &initialState) const {

        initialState = automaton->initial->hash;

        // iterate through all states, build the transition table state for state
        struct trans *trans;
        for (auto s = automaton->initial; s != NULL; s = s->next) {

            // remember each state by its hash, also count them
            states.push_back(s->hash);

            // remember the state as being final if it is
            if(s->accept) {
                finalStates.insert(s->hash);
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
                }
            }

            transitionTable.insert(std::make_pair(s->hash, stateTable));
        }
    }




private:
    FiniteAutomaton(struct fa *FA) {
        automaton = FA;
        epsilon_closed = 0;
    }

    struct fa* automaton; // the actual language

    // this flag won't just save you running time, it will also save you hours of debugging time
    int epsilon_closed;

};
