#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <forward_list>
extern "C" {
    #include <fa.h>
}
#include <assert.h>

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

//        std::cout << "in FiniteAutomaton::courcelle_construct..." << std::endl;
        /*
         * build two sets that hold the components required for the bipartite construction
         */
        std::set<int> lhComponents;
        std::set<int> rhComponents;

        for(auto it = LHStoRHS.begin(); it != LHStoRHS.end(); it++) {
            lhComponents.insert(it->first);

            for(auto rhc: it->second) {
                rhComponents.insert(rhc);
//                std::cout << "lhs\t" << it->first << "\tmapped to rhs\t" << rhc << std::endl;
            }
        }

//        std::cout << "finding lhComponents and rhComponents done..." << std::endl;
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

        if(comp == 3) {
            FILE * initialFile;
            initialFile = fopen ("courcelle_3_initial.dot","w");
            fa_dot(initialFile, finalAutomaton);
            fclose (initialFile);
        }

//        std::cout << "remembered finalAutomaton and its initial state..." << std::endl;

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

        if(comp == 3) {
            FILE * mergedLHS;
            mergedLHS = fopen ("courcelle_3_merged_lhs.dot","w");
            fa_dot(mergedLHS, finalAutomaton);
            fclose (mergedLHS);
        }

//        std::cout << "merged left hand clones into construction..." << std::endl;

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


        if(comp == 3) {
            FILE * mergedRHS;
            mergedRHS = fopen ("courcelle_3_merged_rhs.dot","w");
            fa_dot(mergedRHS, finalAutomaton);
            fclose (mergedRHS);
        }

//        std::cout << "merged right hand clones into construction..." << std::endl;

        // merge the right hand alphabet automaton into the construction
        struct state* target = righthandAlphabetStar.automaton->initial;
        assert (target->next == NULL); // otherwise, the construction doesn't work as intended

        fa_merge(finalAutomaton, &(righthandAlphabetStar.automaton));
//        std::cout << "merged right hand alphabet clone into construction..." << std::endl;


        if(comp == 3) {
            FILE * mergedRAlph;
            mergedRAlph = fopen ("courcelle_3_merged_ralph.dot","w");
            fa_dot(mergedRAlph, finalAutomaton);
            fclose (mergedRAlph);
        }

        // add the epsilon transitions from the final states of the right hand automata to the
        // state of the right hand alphabet automaton
        for(auto rhFinal: rhComponentsFinalStates) {
            add_epsilon_trans(rhFinal, target);
        }

        if(comp == 3) {
            FILE * transrhs;
            transrhs = fopen ("courcelle_3_trans_rhs.dot","w");
            fa_dot(transrhs, finalAutomaton);
            fclose (transrhs);
        }

//        std::cout << "added epsilon trans rhFinal->target..." << std::endl;

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

        if(comp == 3) {
            FILE * transbip;
            transbip = fopen ("courcelle_3_trans_bip.dot","w");
            fa_dot(transbip, finalAutomaton);
            fclose (transbip);
        }

        // add the transitions from left hand alphabet automaton to the initial states of the left
        // hand automata of the bipartite construction
        for(auto lhInitial: lhComponentsInitialState) {
            add_epsilon_trans(source, lhInitial);
        }
//        std::cout << "added epsilon trans source->lhInitial..." << std::endl;

        if(comp == 3) {
            FILE * translhs;
            translhs = fopen ("courcelle_3_trans_lhs.dot","w");
            fa_dot(translhs, finalAutomaton);
            fclose (translhs);
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

        if(comp == 3) {
            FILE * linTerms;
            linTerms = fopen ("courcelle_3_lin_terms.dot","w");
            fa_dot(linTerms, finalAutomaton);
            fclose (linTerms);
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

            if(comp == 3) {
                FILE * constants;
                constants = fopen ("courcelle_3_constants.dot","w");
                fa_dot(constants, finalAutomaton);
                fclose (constants);
            }
        }

        finalAutomaton->minimal = 0;
        finalAutomaton->deterministic = 0;
//        std::cout << "leaving FiniteAutomaton::courcelle_construct..." << std::endl;

        return FiniteAutomaton(finalAutomaton);
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
//        std::cout << "doing epsilon closure of language" << std::endl;
//        std::cout << "language:\t" << string() << std::endl;
//        std::cout << "size of automaton:\t" << size() << std::endl;
        if(empty() || epsilon_closed == 1) {
            return *this;
        }
//        std::cout << "failed here if no more text\n\n\n" << std::endl;

        std::set<unsigned long> foundStates;
        std::vector<state*> bfsQueue;

        struct trans *trans;
        struct state *next;
        for (auto s = automaton->initial; s != NULL; s = s->next){
            foundStates.clear();
            bfsQueue.clear();

            foundStates.insert(s->hash);
            bfsQueue.push_back(s);
            std::set<state*> targets;

            while(!bfsQueue.empty()) {
//                std::cout << "getting last element of queue" << std::endl;
                next = bfsQueue.back();
//                std::cout << "\"next\" state hash:"<< std::endl;
//                std::cout <<  next->hash << std::endl;
//                std::cout << "popping queue" << std::endl;
                bfsQueue.pop_back();
//                std::cout << "getting transition array of \"next\" state" << std::endl;
                trans = next->trans;
//                std::cout << "number of transitions of \"next\":\t" << next->tused << std::endl;

                // remove this line to earn yourself hours of debugging:
                int iterations = next->tused;
                std::cout << "next->tused = iterations = " << iterations << std::endl;


                // iterate over all transitions
                for(int i = 0; i < iterations; i++) {
//                    std::cout << "iteration:\t" << i << "\thash of \"next\":\t" << next->hash << std::endl;

                    if(foundStates.count((trans+i)->to->hash) == 0) {
                        foundStates.insert((trans+i)->to->hash);
//                        std::cout << "inserted into found states:" << std::endl;
                        bfsQueue.push_back((trans+i)->to);
//                        std::cout << "pushed into queue" << std::endl;

//                        std::cout << "hash target state:" << std::endl;
//                        std::cout << (trans+i)->to->hash << std::endl;


                        // debugging for the thousandth time
                        FILE * file;
                        file = fopen ("analyze.dot","w");
                        fa_dot(file, automaton);
                        fclose (file);
                        printf("%p to %p\n", next, (trans+i)->to);

                        //remember the new target state
                        targets.insert((trans+i)->to);
                    }
//                    std::cout << " " << std::endl;
                }

//                std::cout << "\n\n\n\n" << std::endl;
            }

            // epsilon-close
            for(auto target: targets) {
                add_epsilon_trans(s, target);
            }
        }
//        std::cout << "nope, wasnt that, new language:\t" << string() << "\n\n\n" << std::endl;

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
//                    std::cout << "letter in alphabet:\t" << c << std::endl;
                }
            }
        }

        return alphabet;
    }

    std::map<int, std::set<std::string>> prefixesToMaxLength(int maxLength) const {
        std::map<int, std::map<unsigned long, std::set<std::string>>> prefixesViaStates;
        std::map<unsigned long, struct state*> hashToState;
        std::map<int, std::set<std::string>> prefixes;

//        std::cout << "automaton:\t" << string() << std::endl;
//
//        std::cout << "length:\t" << maxLength << std::endl;

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
                prefixesViaStates[1][(trans+i)->to->hash].insert(prefix);
                prefixes[1].insert(prefix);
//                std::cout << "new prefix length 1:\t" << prefix << std::endl;
            }
        }
//
//        std::cout << "prefixes of length 1:" << std::endl;
//        for(auto &prefix: prefixes[1]) {
//            std::cout << prefix << std::endl;
//        }
//        std::cout << "end prefixes of length 1" << std::endl;
//        std::cout << "prefixesViaStates:" << std::endl;
//
//        for(auto &outerEntry: prefixesViaStates) {
//            for(auto &innerEntry: outerEntry.second) {
//                for(auto &prefix: innerEntry.second) {
//                    std::cout << "prefixesViaStates[" << outerEntry.first << "][" << innerEntry.first << "]:\t" << prefix << std::endl;
//                }
//            }
//        }
//
//       std::cout << "prefixesViaStates end." << std::endl;

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
     * Calculates the intersection of the CFG "equations" with start symbol "oldS" and the
     * regular language represented by this finite automaton. The result is a CFG given by
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
    template<typename SR>
    std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> intersectionWithCFG
        (VarId &newS, const VarId &oldS,
                const std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> &oldGrammar) const {

        // for the new grammar
        std::vector<std::pair<VarId, NonCommutativePolynomial<SR>>> resultGrammar;

        // we don't want to start out with useless stuff because the intersection grammar
        // will blow up anyway
        std::queue<VarId> worklist;
        worklist.push(oldS);
        auto workGrammar = NonCommutativePolynomial<SR>::cleanSystem(oldGrammar, worklist);


        // change the grammar to one where monomials have degree at most 2 and those monomials with degree 2
        // have the form XY
        std::map<VarId, SR> variablesToConstants;
        workGrammar = NonCommutativePolynomial<SR>::eliminateTerminalsInNonterminalProductions
                (workGrammar, variablesToConstants, false);
        workGrammar = NonCommutativePolynomial<SR>::binarizeProductions(workGrammar, variablesToConstants);

        // if the grammar doesn't produce anything, there is no need to do anything else;
        // return an empty grammar
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
                newVariables[i][j].resize(workGrammar.size());
            }
        }

        // generate the variables of the new grammar
        for(unsigned long i = 0; i < numberOfStates; i++) {
            for(unsigned long j = 0; j < numberOfStates; j++) {
                for(unsigned long k = 0; k < workGrammar.size(); k++) {
                    std::stringstream ss;
                    ss << "<" << i << "," << Var::GetVar(workGrammar[k].first).string() << "," << j << ">";
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


//        std::cout << "the problem must lie in the actual intersection algorithm" << std::endl;

        // the indices k, i, j are intended this way; we first choose some nonterminal of the grammar
        // and then generate all productions derived from that nonterminal; this mirrors the way
        // in which the algorithm is described in the paper referred to above
        NonCommutativePolynomial<SR> poly;
        for(unsigned long k = 0; k < workGrammar.size(); k++) {
            for(unsigned long i = 0; i < numberOfStates; i++) {
                for(unsigned long j = 0; j < numberOfStates; j++) {
//                    std::cout << "k, i, j: " << k << ", " << i << ", " << j << std::endl;
                    poly = workGrammar[k].second.intersectionPolynomial
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
        epsilon_closed = 0;
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

    // this flag won't just save you running time, it will also save you hours of debugging time
    int epsilon_closed;
};
