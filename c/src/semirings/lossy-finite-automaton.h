#pragma once

#include <string>
#include <memory>
#include <unordered_map>
#include <queue>
#include <string>
#include <ctype.h>
#include <algorithm>
extern "C" {
    #include <fa.h>
}

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/free-structure.h"
#include "../datastructs/finite_automaton.h"
#include "../polynomials/non_commutative_polynomial.h"

#include "lossy-semiring.h"

template<typename SR>
class Evaluator;

/*
 * Models regular languages with the option to make them lossy, i.e. give each symbol of the alphabet the property
 * a = a | epsilon. The languages that result from calling a constructor of LossyFiniteAutomaton do NOT have that
 * property originally; if you want to add the property, call lossify() - this will not change the automaton
 * on which you call the method, but construct another automaton that is the "lossified" version of the original one.
 */
class LossyFiniteAutomaton: public LossySemiring<LossyFiniteAutomaton>  {
public:


    static int maxStates;

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


//        FiniteAutomaton fa, temp;
//        fa = FiniteAutomaton::epsilon();
//
//        // iterate over the string; for each character, construct an automaton that accepts exactly that character
//        // or epsilon; concatenate all automata thus constructed
//        for(int i = 0; i < regularExpression.size(); i++) {
//
//            // the "union with {epsilon}" is crucial here since it models lossyness
//            temp = FiniteAutomaton(std::string(1, regularExpression[i])).unionWith(FiniteAutomaton::epsilon());
////            std::cout << "temp:\t" + temp.string() << std::endl;
//            fa = fa.concatenate(temp).minimize();
////            std::cout << "fa:\t" + fa.string() << std::endl;
//        }

//        std::cout << "LFA(str) begin" << std::endl;
        language = FiniteAutomaton(regularExpression).minimize();

//        std::cout << "LFA(str):\t" + language.string() << std::endl;

//        if(!language.empty()) {
//            std::cout << "lossifiedRegex:\t" + lossifiedRegex(regularExpression) << std::endl;
//        }
    }

    /*
     * Same as LossyFiniteAutomaton(string), really.
     */
    LossyFiniteAutomaton(const VarId var) {
//        std::string regularExpression = Var::GetVar(var).string();
//        FiniteAutomaton fa, temp;

        // iterate over the string; for each character, construct an automaton that accepts exactly that character
        // or epsilon; concatenate all automata thus constructed
//        for(int i = 0; i < regularExpression.size(); i++) {
//            // the "union with {epsilon}" is crucial here since it models lossyness
//            temp = FiniteAutomaton(std::string(1, regularExpression[i])).unionWith(FiniteAutomaton::epsilon());
////            std::cout << "temp:\t" + temp.string() << std::endl;
//            fa = fa.concatenate(temp);
////            std::cout << "fa:\t" + fa.string() << std::endl;
//        }
//        std::cout << "LFA(var) begin" << std::endl;
        FiniteAutomaton fa = FiniteAutomaton(Var::GetVar(var).string());
        maxStates = std::max(maxStates,fa.size());
        language = fa.minimize();
//        std::cout << "LFA(var):\t" + language.string() << std::endl;
    }

    LossyFiniteAutomaton lossify() const {
        if ((*this) == null()) {
            return *this;
        }

        return LossyFiniteAutomaton(lossifiedRegex(string())).minimize();
    }

    static LossyFiniteAutomaton null() {
//        std::cout << "empty:\t" + EMPTY.string() << std::endl;
        return EMPTY;
    }

    static LossyFiniteAutomaton one() {
//        std::cout << "epsilon:\t" + EPSILON.string() << std::endl;
        return EPSILON;
    }

    /*
     * Lossifies a valid regular expression, i.e. every alphabet symbol a will be replaced by (a|()).
     */
    static std::string lossifiedRegex(std::string regex) {

//        std::cout << "regex before lossification:\t\"" + regex +"\""<< std::endl;
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
//        std::cout << "regex after lossification:\t" + lossyRegex << std::endl;
        return lossyRegex;
    }

    LossyFiniteAutomaton minimize() {
        auto newLanguage = language.minimize();
        maxStates = std::max(maxStates,newLanguage.size());
        return LossyFiniteAutomaton{newLanguage};
    }

    LossyFiniteAutomaton star() const {
        auto newLanguage = language.kleeneStar();
        maxStates = std::max(maxStates,newLanguage.size());
        return LossyFiniteAutomaton(newLanguage.minimize());
    }

    LossyFiniteAutomaton complement() const {
        return LossyFiniteAutomaton(language.complement().minimize());
    }

    LossyFiniteAutomaton minus(const LossyFiniteAutomaton &x) const {
        auto newLanguage = language.minus(x.language);
        maxStates = std::max(maxStates,newLanguage.size());
        return LossyFiniteAutomaton(newLanguage).minimize();
    }

    bool contains(const LossyFiniteAutomaton &other) const {
        return language.contains(other.language);
    }


    LossyFiniteAutomaton operator+(const LossyFiniteAutomaton &x) {
        auto newLanguage = language.unionWith(x.language);
        maxStates = std::max(maxStates,newLanguage.size());
        return LossyFiniteAutomaton {newLanguage.minimize()};
    }

    LossyFiniteAutomaton& operator+=(const LossyFiniteAutomaton &x) {
        language = language.unionWith(x.language);
        maxStates = std::max(maxStates,language.size());
        language = language.minimize();
        return *this;
    }

    LossyFiniteAutomaton operator*(const LossyFiniteAutomaton &x) {
        auto newLanguage = language.concatenate(x.language);
        maxStates = std::max(maxStates,newLanguage.size());
        return LossyFiniteAutomaton {newLanguage.minimize()};
    }

    LossyFiniteAutomaton& operator*=(const LossyFiniteAutomaton &x) {
        language = language.concatenate(x.language);
        maxStates = std::max(maxStates,language.size());
        language = language.minimize();
        return *this;
    }

    /*
     * Checks if the LOSSY languages represented by the two automata are equal.
     */
    bool operator==(const LossyFiniteAutomaton &x) const {
//        std::cout << "LFA==" << std::endl;

        // we only do this because there is no special symbol for the empty language, and we don't want to
        // lossify the string representation of the empty language: "__EMPTYLANGUAGE__"
        if(x.language.empty() && language.empty()) {
//            std::cout << "both empty" << std::endl;
//            std::cout << "LFA== both empty" << std::endl;
            return true;
        } else if(x.language.empty() || language.empty()) {
//            std::cout << "exactly one empty" << std::endl;
//            std::cout << "LFA== exactly one empty" << std::endl;
            return false;
        }

        // this is to avoid unnecessary lossification of the languages
        if(x.language.equals(FiniteAutomaton::epsilon()) && language.equals(FiniteAutomaton::epsilon())) {
//            std::cout << "LFA== both epsilon" << std::endl;
            return true;
        } else if(x.language.equals(FiniteAutomaton::epsilon()) || language.equals(FiniteAutomaton::epsilon())) {
//            std::cout << "LFA== exactly one epsilon" << std::endl;
            return false;
        }

//        std::cout << "we actually need to lossify, languages:" << std::endl;
//        std::cout << "this:\t" << string() << std::endl;
//        std::cout << "other:\t" << x.string() << std::endl;
//        std::cout << "and here goes lossification..." << std::endl;
//        std::cout << "LFA==: lossifying first" << std::endl;
//        auto fa = lossify();
//        std::cout << "LFA==: lossifying second" << std::endl;
//        auto other = x.lossify();
//        std::cout << "LFA== done lossifying" << std::endl;
//        std::cout << "language of [(a|())(b|())]:\t";
//        LossyFiniteAutomaton lossyAorB = LossyFiniteAutomaton("[(a|())(b|())]");
//
//        std::cout << lossyAorB.string() << std::endl;
//        std::cout << "language of (a|()){1,15}:\t";
//        LossyFiniteAutomaton lossyAlimitedIteration = LossyFiniteAutomaton("(a|()){1,15}");
//        std::cout << lossyAlimitedIteration.string() << std::endl;
//        std::cout << "lossify a{1,15}:\t";
//        std::cout << lossifiedRegex("a{1,15}") << std::endl;
//        std::cout << "language of [axb]{1,5} automaton:" << std::endl;
//        LossyFiniteAutomaton axbchoice = LossyFiniteAutomaton("[axb]{1,5}");
//        std::cout << axbchoice.string() << std::endl;
//        std::cout << "lossify [axb]{1,15}:\t";
//        std::cout << lossifiedRegex("[axb]{1,5}") << std::endl;
//        std::cout << "language of lossified [axb]{1,5} automaton:" << std::endl;
//        LossyFiniteAutomaton lossyaxbchoice = LossyFiniteAutomaton(lossifiedRegex("[axb]{1,5}"));
//        std::cout << lossyaxbchoice.string() << std::endl;
//        std::string hi("hi");
//        std::cout << "hi[1] == 'i':\t" << (hi[1] == 'i') << std::endl;

//        std::cout << "end of LFA==" << std::endl;
        return lossify().language.equals(x.lossify().language);
    }

    std::string string() const {
        return language.string();
    }

    std::string lossyString() const {
        return lossifiedRegex(string());
    }

    int size() const {
        return language.size();
    }

    /*
     * Calculates the intersection of the CFG "equations" with start symbol "oldS" and the
     * regular language represented by this finite automaton. The result is returned as a CFG given by
     * its set of productions while its start symbol is stored in "newS".
     */
    std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> intersectionWithCFG
        (VarId &newS, const VarId &oldS, const std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &equations) const {
        return language.intersectionWithCFG(newS, oldS, equations);
    }

    /*
     * Gives you the alphabet used in this automaton.
     */
    std::set<unsigned char> alphabet() const {
        return language.alphabet();
    }

    std::map<int, std::set<std::string>> prefixesToMaxLength(int maxLength) const {
        return language.prefixesToMaxLength(maxLength);
    }

    /*
     * Quick test stuff. On input of grammar <S> ::= "a"<S> | "b"; this should give you "shortest word in intersection: aaab".
     */
    static void intersectionTest(std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> &equations) {
        VarId oldS;

        for(auto &equation: equations) {
            if(Var::GetVar(equation.first).string() == "S") {
                oldS = equation.first;
//                std::cout << "found S:\t" + Var::GetVar(oldS).string() << std::endl;
                break;
            }
        }

        std::cout << "Grammar:" << std::endl;
        for(auto &equation: equations) {
            std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
        }

        LossyFiniteAutomaton dummy = LossyFiniteAutomaton("a([efg]*|bo)(l|k)c*|h*");
        std::cout << "regular language:\t" << dummy.string() << std::endl;
        VarId newS;
        std::vector<std::pair<VarId, NonCommutativePolynomial<LossyFiniteAutomaton>>> intersection =
                dummy.intersectionWithCFG(newS, oldS, equations);

//        std::cout << "raw intersection system, start symbol: " << Var::GetVar(newS).string() << std::endl;
//        for(auto &equation: intersection) {
//            std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//        }

        std::queue<VarId> worklist;
        worklist.push(newS);
        auto cleanIntersection = NonCommutativePolynomial<LossyFiniteAutomaton>::cleanSystem(intersection, worklist);

        std::cout << "clean intersection system, start symbol: " << Var::GetVar(newS).string() << std::endl;
        for(auto &equation: cleanIntersection) {
            std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
        }

        LossyFiniteAutomaton shortestWord = NonCommutativePolynomial<LossyFiniteAutomaton>::shortestWord(cleanIntersection, newS);
        std::cout << "shortest word in intersection:\t" + shortestWord.string() << std::endl;

        std::cout << "Testing regex []*:" << std::endl;
        LossyFiniteAutomaton testLFA = LossyFiniteAutomaton("[]*");
        std::cout << testLFA.string() << std::endl;
        std::cout << "Testing regex []:" << std::endl;
        LossyFiniteAutomaton testLFA2 = LossyFiniteAutomaton("[]*");
        std::cout << testLFA2.string() << std::endl;
    }

private:
    FiniteAutomaton language;

    LossyFiniteAutomaton(FiniteAutomaton fa) : language(fa){}

    static LossyFiniteAutomaton EMPTY;
    static LossyFiniteAutomaton EPSILON;
    //static LossyFiniteAutomaton UNIVERSE;

    friend struct std::hash<LossyFiniteAutomaton>;
};
