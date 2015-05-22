#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <string>
#include <ctype.h>
#include <algorithm>


#include "../datastructs/var.h"
#include "../datastructs/finite_automaton.h"
#include "../semirings/semiring.h"


/*
 * Models regular languages with the option to make them lossy, i.e. give each symbol of the alphabet the property
 * a = a | epsilon. The languages that result from calling a constructor of LossyFiniteAutomaton do NOT have that
 * property originally; if you want to add the property, call lossify() - this will not change the automaton
 * on which you call the method, but access another automaton that is the "lossified" version of the original one.
 *
 * Note: equality checks on these languages are with respect to the lossified versions only.
 */
class LossyFiniteAutomaton: public StarableSemiring<LossyFiniteAutomaton, Commutativity::NonCommutative, Idempotence::Idempotent> {
public:


    /* Default constructor creates one (multiplicative neutral) element. */
    LossyFiniteAutomaton() {
        language = EPSILON.language;
    }

    /*
     * Put in a POSIX regular expression and this will give you a shiny new LossyFiniteAutomaton representing
     * the lossyfied language of that expression. Use | for choice, () for grouping, * for Kleene star;
     * don't use spaces, also concatenation does not require a symbol.
     */
    LossyFiniteAutomaton(const std::string regularExpression) {
        language = FiniteAutomaton(regularExpression).minimize();
    }

    /*
     * Same as LossyFiniteAutomaton(string), really.
     */
    LossyFiniteAutomaton(const VarId var) {
        language = FiniteAutomaton(Var::GetVar(var).string()).minimize();
    }

    LossyFiniteAutomaton lossify() const {
        if (language.empty()) {
            return *this;
        }

        return LossyFiniteAutomaton(language.epsilonClosure()).minimize();
    }

    static LossyFiniteAutomaton null() {
        return LossyFiniteAutomaton(FiniteAutomaton());
    }

    static LossyFiniteAutomaton one() {
        return LossyFiniteAutomaton(FiniteAutomaton::epsilon());
    }

    /*
     * Lossifies a valid regular expression, i.e. every alphabet symbol a will be replaced by (a|()).
     */
    static std::string lossifiedRegex(std::string regex) {

        std::stringstream ss;

        for(int i = 0; i < regex.size(); i++) {

            // we want to be able to lossify strings that are given by libfa; those are extended
            // POSIX regexes, which may contain [] as character set and {m,n} as iteration, and those
            // need to be treated slightly differently
            if(regex[i] == '{') { // iteration groups will not require lossification themselves since
                                  // the expression they iterate will be lossified already
                while(regex[i] != '}') {
                    ss << regex[i];
                    i++;
                }

                ss << regex[i]; // don't forget the '}'
            } else if (regex[i] == '[') {
                ss << "("; // we can't simply lossify every element of the character set since that
                           // would just put the elements '(', ')' and '|' into the set without
                           // giving the possibility to produce epsilon; so instead, we lossify
                           // the whole set

                while(regex[i] != ']') { // leave the character set itself unchanged
                    ss << regex[i];
                    i++;
                }

                ss << regex[i] << "|())"; // don't forget the ']'
            } else if(isalnum(regex[i])) {
                ss << "(" << regex[i] << "|())"; // lossify alphabet characters
            } else {
                ss << regex[i]; // leave "syntax characters" unchanged
            }
        }

        auto lossyRegex = ss.str();
        return lossyRegex;
    }

    static LossyFiniteAutomaton courcelle_construct(LossyFiniteAutomaton lefthandAlphabetStar,
            LossyFiniteAutomaton righthandAlphabetStar,
            std::unordered_map<int, std::set<int>> &LHStoRHS,
            LossyFiniteAutomaton closureOfConstantMonomials,
            std::set<int> &lowerLinearTerms,
            std::unordered_map<int, LossyFiniteAutomaton> &componentToClosure,
            int comp) {

        // convert the map so FiniteAutomaton doesn't need to know about LossyFiniteAutomaton
        std::unordered_map<int, FiniteAutomaton> componentToClosureAutomaton;
        for(auto it = componentToClosure.begin(); it != componentToClosure.end(); it++) {
            componentToClosureAutomaton[it->first] = (it->second).language;
        }

        return LossyFiniteAutomaton(FiniteAutomaton::courcelle_construct(lefthandAlphabetStar.language,
                righthandAlphabetStar.language, LHStoRHS, closureOfConstantMonomials.language,
                lowerLinearTerms, componentToClosureAutomaton, comp));
    }

    LossyFiniteAutomaton minimize() {
        return LossyFiniteAutomaton(language.minimize());
    }

    LossyFiniteAutomaton star() const {
        return LossyFiniteAutomaton(language.kleeneStar().minimize());
    }

    LossyFiniteAutomaton complement() const {
        return LossyFiniteAutomaton(language.complement().minimize());
    }

    LossyFiniteAutomaton minus(const LossyFiniteAutomaton &x) const {
        return LossyFiniteAutomaton(language.minus(x.language).minimize());
    }

    bool contains(const LossyFiniteAutomaton &other) const {
        return language.contains(other.language);
    }


    LossyFiniteAutomaton operator+(const LossyFiniteAutomaton &x) {
        return LossyFiniteAutomaton (language.unionWith(x.language).minimize());
    }

    LossyFiniteAutomaton& operator+=(const LossyFiniteAutomaton &x) {
        language = language.unionWith(x.language).minimize();
        return *this;
    }

    LossyFiniteAutomaton operator*(const LossyFiniteAutomaton &x) {
        return LossyFiniteAutomaton (language.concatenate(x.language).minimize());
    }

    LossyFiniteAutomaton& operator*=(const LossyFiniteAutomaton &x) {
        language = language.concatenate(x.language).minimize();
        return *this;
    }

    /*
     * Checks if the LOSSY languages represented by the two automata are equal.
     */
    bool operator==(const LossyFiniteAutomaton &x) const {

        // we only do this because there is no special symbol for the empty language, and we don't want to
        // lossify the string representation of the empty language: "__EMPTYLANGUAGE__"
        if(x.language.empty() && language.empty()) {
            return true;
        } else if(x.language.empty() || language.empty()) {
            return false;
        }

        // this is to avoid unnecessary lossification of the languages
        if(x.language.equals(FiniteAutomaton::epsilon()) && language.equals(FiniteAutomaton::epsilon())) {
            return true;
        } else if(x.language.equals(FiniteAutomaton::epsilon()) || language.equals(FiniteAutomaton::epsilon())) {
            return false;
        }

        return lossify().language.equals(x.lossify().language);
    }

    std::string string() const {
        return language.string();
    }

    void write_dot_file(FILE * file) const {
        language.write_dot_file(file);
    }

    std::string lossyString() const {
        return lossifiedRegex(string());
    }

    int size() const {
        return language.size();
    }

    FiniteAutomaton getFA() {
      return language;
    }

    /*
     * Gives you the alphabet used in this automaton.
     */
    std::set<unsigned char> alphabet() const {
        return language.alphabet();
    }

    std::unordered_map<int, std::set<std::string>> prefixesToMaxLength(int maxLength, std::set<char> &derivableFirstLetters) const {
        return language.prefixesToMaxLength(maxLength, derivableFirstLetters);
    }



private:
    FiniteAutomaton language;

    LossyFiniteAutomaton(FiniteAutomaton fa) : language(fa){}

    static LossyFiniteAutomaton EMPTY;
    static LossyFiniteAutomaton EPSILON;

    friend struct std::hash<LossyFiniteAutomaton>;

};
