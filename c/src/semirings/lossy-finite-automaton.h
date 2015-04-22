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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graphviz.hpp>

template<typename LFA>
class Evaluator;

template <typename LFA>
struct Vertex {
    VarId var; // var and rex combines the equations in the vertex
    NonCommutativePolynomial<LFA> rex;
};

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
            std::map<int, std::set<int>> &LHStoRHS,
            LossyFiniteAutomaton closureOfConstantMonomials,
            std::set<int> &lowerLinearTerms,
            std::map<int, LossyFiniteAutomaton> &componentToClosure,
            int comp) {

        // convert the map so FiniteAutomaton doesn't need to know about LossyFiniteAutomaton
        std::map<int, FiniteAutomaton> componentToClosureAutomaton;
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

    /*
     * Calculates the intersection of the CFG "equations" with start symbol "oldS" and the
     * regular language represented by this finite automaton. The result is returned as a CFG given by
     * its set of productions while its start symbol is stored in "newS".
     */
    std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> intersectionWithCFG
        (VarId &newS, const VarId &oldS, const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations) const {
        return language.intersectionWithCFG(newS, oldS, equations);
    }

    /*
     * Gives you the alphabet used in this automaton.
     */
    std::set<unsigned char> alphabet() const {
        return language.alphabet();
    }

    std::map<int, std::set<std::string>> prefixesToMaxLength(int maxLength, std::set<char> &derivableFirstLetters) const {
        return language.prefixesToMaxLength(maxLength, derivableFirstLetters);
    }

    /*
     * Solves a polynomial system over the lossy semiring. For derivation of
     * the algorithm and proof of correctness, see "Esparza, Kiefer, Luttenberger:
     * Derivation Tree Analysis for Accelerated Fixed-Point Computation".
     *
     * If the system is not clean, i.e. there are polynomials in it with a fixpoint of 0, then you must set
     * systemNeedsCleaning to true in order for the algorithm to work as intended.
     *
     * If the system is not in lossy QNF, then you need to set systemNeedsLossyQNF to true.
     *
     * NOTE: not currently wired up anywhere.
     */
    static  ValuationMap<LossyFiniteAutomaton> solvePolynomialSystem(
            std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations,
            bool systemNeedsCleaning,
            bool systemNeedsLossyQNF) {

        ValuationMap<LossyFiniteAutomaton> solution;


        // start out with a clean system
        if(systemNeedsCleaning) {
            std::queue<VarId> worklist;
            equations = cleanSystem(equations, worklist);
        }

        // bring the clean system into qnf, if necessary
        std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> qnf;

        if(systemNeedsLossyQNF) {
            qnf = quadraticNormalForm(equations, false);
            lossifyQNF(qnf);
        } else {
            qnf = equations;
        }

        // build a zero vector for the evaluation of f
        ValuationMap<LossyFiniteAutomaton> zeroSystem;
        for(auto &equation: qnf) {
            zeroSystem.insert({equation.first, LossyFiniteAutomaton::null()});
        }

        // build the vectors f(0) and f^n(0) where n is the number of equations in the system
        int times = 1;
        ValuationMap<LossyFiniteAutomaton> f_0 = evaluateSystem(qnf, times, zeroSystem);
        times = qnf.size() - 1;
        ValuationMap<LossyFiniteAutomaton> f_n_0 = evaluateSystem(qnf, times, f_0);

        // find the LossySemiring element in the "middle" of the expression
        LossyFiniteAutomaton middle = LossyFiniteAutomaton::null();
        for(auto &elem_mapping: f_0) {
            middle += elem_mapping.second;
        }

        // build the differential of the system
        ValuationMap<LossyNonCommutativePolynomial> differential;
        for(auto &equation: qnf) {
                differential.insert(std::make_pair(equation.first, equation.second.differential_at(f_n_0)));
        }

        // sum all polynomials of the differential of the system; we don't need attribution
        // of each polynomial to the respective variable since we only care about the
        // leading and trailing coefficients of each monomial
        LossyNonCommutativePolynomial differential_sum = LossyNonCommutativePolynomial::null();
        for(auto &equation: differential) {
            differential_sum += equation.second;
        }

        // get the lefthand and righthand semiring element of the fixpoint
        LossyFiniteAutomaton lefthandSum = differential_sum.getSumOfLeadingFactors();

        LossyFiniteAutomaton righthandSum = differential_sum.getSumOfTrailingFactors();
        LossyFiniteAutomaton fixpoint = lefthandSum.star() * middle * righthandSum.star();

        for(auto &variable_mapping: qnf) {
            solution.insert(std::make_pair(variable_mapping.first, fixpoint));
        }

        return solution;
    }

    /*
     * Takes two grammars and refines them by intersecting them with the prefix languages of their shared downward closure.
     */
    static LossyFiniteAutomaton refineCourcelle(const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations_1,
            const VarId &S_1, const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations_2,
            const VarId &S_2, int maxLengthOfPrefixes) {

        LossyFiniteAutomaton approx_1 = downwardClosureCourcelle(equations_1, S_1);
        LossyFiniteAutomaton approx_2 = downwardClosureCourcelle(equations_2, S_2);

        bool A1_subset_A2 = approx_2.contains(approx_1);
        bool A2_subset_A1 = approx_1.contains(approx_2);

        if(A1_subset_A2 && A2_subset_A1) {
            LossyFiniteAutomaton difference = LossyFiniteAutomaton::null();

            if(!(approx_1 == LossyFiniteAutomaton::null()) && !(approx_1 == LossyFiniteAutomaton::one())) {

                // build the automaton Sigma*
                std::string regexAlphabetStar = "[";
                for (auto c: approx_1.alphabet()) {
                    regexAlphabetStar.append(1, c);
                }
                regexAlphabetStar += "]*";

                // build the map of (length, prefixes of that length)
                std::set<char> derivableFirstLetters = LossyNonCommutativePolynomial::getDerivableFirstLettes(equations_1, S_1);
                std::set<char> derivableFirstLetters2 = LossyNonCommutativePolynomial::getDerivableFirstLettes(equations_2, S_2);

                for(char letter: derivableFirstLetters2) {
                    derivableFirstLetters.insert(letter);
                }

                std::unordered_map<int, std::set<std::string>> prefixesPerLength = approx_1.prefixesToMaxLength(maxLengthOfPrefixes, derivableFirstLetters);

                // iterate over the prefixes
                bool noDifference = true;
                std::queue<VarId> worklist;
                for(int i = 1; (i <= maxLengthOfPrefixes) && noDifference; i++) {
                    std::cout << "refining " << i << "..." << std::endl;

                    for(auto &prefix: prefixesPerLength[i]) {
                        std::string prefixAlphaStar = std::string(prefix) + regexAlphabetStar;
                        LossyFiniteAutomaton prefixAutomaton(prefixAlphaStar);

                        VarId S_1_partition, S_2_partition;

                        // generate the subset of language 1
                        auto equations_1_Partition = prefixAutomaton.intersectionWithCFG(S_1_partition, S_1, equations_1);

                        while(!worklist.empty()) {
                            worklist.pop();
                        }

                        worklist.push(S_1_partition);
                        equations_1_Partition = cleanSystem(equations_1_Partition, worklist);

                        // generate the subset of language 2
                        auto equations_2_Partition = prefixAutomaton.intersectionWithCFG(S_2_partition, S_2, equations_2);

                        while(!worklist.empty()) {
                            worklist.pop();
                        }

                        worklist.push(S_2_partition);
                        equations_2_Partition = cleanSystem(equations_2_Partition, worklist);

                        // approximate the grammars
                        LossyFiniteAutomaton approx_1_part = downwardClosureCourcelle(equations_1_Partition, S_1_partition);
                        LossyFiniteAutomaton approx_2_part = downwardClosureCourcelle(equations_2_Partition, S_2_partition);

                        LossyFiniteAutomaton tempDiff = compareClosures(equations_1_Partition, S_1_partition, equations_2_Partition, S_2_partition, approx_1_part, approx_2_part);

                        if(!(tempDiff == LossyFiniteAutomaton::null())) {
                            difference = tempDiff;
                            noDifference = false;
                            std::cout << "different " << i << std::endl;
                            break;
                        }
                    }
                }
            } else { // both languages empty
                std::cout << "equal 0" << std::endl;
                return difference;
            }

            if(difference == LossyFiniteAutomaton::null()) {
                std::cout << "maybe_equal " << maxLengthOfPrefixes << std::endl;
            }

            return difference;
        } else {
            LossyFiniteAutomaton difference = compareClosures(equations_1, S_1, equations_2, S_2, approx_1, approx_2);
            std::cout << "different 0" << std::endl;

            return difference;
        }
    }

    /*
     * Approximates the language defined by the grammar starting at S given in the system of equations by calculating
     * its downward closure.
     *
     * NOTE: not currently wired up anywhere.
     */
    static LossyFiniteAutomaton downwardClosureDerivationTrees(const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations,
            const VarId &S) {

        std::queue<VarId> worklist;
        worklist.push(S);

        auto cleanEquations = cleanSystem(equations, worklist);

        // if the start symbol is unproductive, then it doesn't generate anything;
        // the downward closure of the empty set is the empty set
        if(!variableHasProduction(cleanEquations, S)) {
            return LossyFiniteAutomaton::null();
        }

        auto cleanQNF = quadraticNormalForm(cleanEquations, false);
        lossifyQNF(cleanQNF);
        auto components = group_by_scc(cleanQNF, false);

        ValuationMap<LossyFiniteAutomaton> knownValuations;
        for(int i = 0; i < components.size(); i++) {

            if(i != 0) {
                for(auto &equation: components[i]) {
                    equation.second = equation.second.partial_eval(knownValuations);
                }
            }

            auto newValuations = solvePolynomialSystem(components[i], false, false);

            for(auto &value: newValuations) {
                knownValuations.insert(std::make_pair(value.first, value.second));
            }
        }

        return knownValuations[S].lossify();
    }

    /*
     * Calculates the downward closure of the grammar defined by "equations" starting at "S"; the algorithm
     * derives from the "On Constructing Obstruction Sets of Words", B. Courcelle, Bulletin of EATCS 1991.
     *
     * WARNING: This function will break if you use a LossyFiniteAutomaton that does not have a constructor that takes POSIX
     * regex string or if there are any terminals in the grammar that contain non-alphanumeric letters,
     * see non_commutative_monomial.get_terminals().
     */
    static LossyFiniteAutomaton downwardClosureCourcelle(const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations,
            const VarId &S) {

        std::queue<VarId> worklist;
        worklist.push(S);

        auto cleanEquations = cleanSystem(equations, worklist);

        // if the start symbol is unproductive, then it doesn't generate anything;
        // the downward closure of the empty set is the empty set
        if(!variableHasProduction(cleanEquations, S)) {
            return LossyFiniteAutomaton::null();
        }

        auto cleanQNF = quadraticNormalForm(cleanEquations, true);
        auto components = group_by_scc(cleanQNF, false);

        // will hold the downward closure of all components
        std::map<int, LossyFiniteAutomaton> componentToClosure;

        // map variables to their components
        std::map<VarId, int> varToComponent = mapVariablesToComponents(components);

        // these components will have the star of the letters reachable from the variables
        // of the respective component as their closure
        std::set<int> squarableComponents = findSquarableComponents(components, varToComponent);

        // we need this map for the case where we can duplicate a variable
        std::map<int, std::set<unsigned char>> componentToReachableLetters =
                findReachableLetters(components, varToComponent);

        // we can already calculate the downward closures of squarable components
        calculateClosuresOfSquarableComponents(componentToReachableLetters, squarableComponents, componentToClosure);

        // for each component, this map holds all pairs AB where A and B are in a lower component;
        // we use the downward closure for those pairs and concatenate accordingly to construct the
        // "middle" part of the downward closure the way Courcelle calculates it
        std::map<int, std::map<int, std::set<int>>> quadraticLHStoRHS = mapQuadraticLHStoRHS(components, varToComponent, squarableComponents);

        // get the components of variables Y_l and Y_r such that there is a monomial in the component of any nonsquarable B
        // that has the form Y_l*B or B*Y_r;
        // their closures will be part of the lefthand/righthand terms of the calculation by Courcelle
        std::map<int, std::set<int>> lhsLowerComponentVariables;
        std::map<int, std::set<int>> rhsLowerComponentVariables;
        calculateLowerComponentVariables(lhsLowerComponentVariables, rhsLowerComponentVariables,
                components, varToComponent, squarableComponents);

        // for the middle part of each component, we also need the sum of the closures of all constant monomials
        // appearing in that component
        std::map<int, LossyFiniteAutomaton> closuresOfConstantMonomials;
        calculateClosuresOfConstantMonomials(closuresOfConstantMonomials, components, squarableComponents);

        // we have everything we can prepare beforehand; everything else will need to be dealt with while
        // putting the closures together
        calculateClosuresOfNonsquarableComponents(componentToClosure, components,
                squarableComponents, quadraticLHStoRHS,
                lhsLowerComponentVariables, rhsLowerComponentVariables,
                closuresOfConstantMonomials, varToComponent, componentToReachableLetters);

        return componentToClosure[varToComponent[S]];
    }

private:
    FiniteAutomaton language;

    LossyFiniteAutomaton(FiniteAutomaton fa) : language(fa){}

    static LossyFiniteAutomaton EMPTY;
    static LossyFiniteAutomaton EPSILON;

    friend struct std::hash<LossyFiniteAutomaton>;

    /*
     * Cleans the polynomial system, i.e. removes variables that are unproductive for all the grammars we get
     * if exactly the symbols in "startSymbols" are considered as initial symbols of a grammar.
     *
     * In the context of fixpoints, if you only want to eliminate symbols that have a 0 fixpoint, hand this function an
     * empty queue "startSymbols".
     */
    static std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> cleanSystem
        (const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations, std::queue<VarId> &startSymbols) {

        if(startSymbols.empty()) {
            std::cout << "start symbols is empty" << std::endl;
            for(auto &equation: equations) {
                startSymbols.push(equation.first);
            }
        } else {
            std::queue<VarId> temp;

            while(!startSymbols.empty()) {
                auto var = startSymbols.front();
                startSymbols.pop();
                temp.push(var);
            }

            swap(temp,startSymbols);
        }

        // just hand the system to the cleaning procedure for grammars and use the set of all
        // variables for which there exist productions as the set of initial nonterminals
        return LossyNonCommutativePolynomial::cleanSystem(equations, startSymbols);
    }

    /*
     * Brings a system into quadratic normal form.
     */
    static std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> quadraticNormalForm
            (const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations,
                    bool eliminateTerminalsInLinearProductions) {
            std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> systemBeingProcessed;

            // for all nonterminal productions, introduce new variables for the terminal factors that
            // appear in them as one would when calculating the CNF of a CFG; do not do this for linear terms since
            // QNF allows for linear terms
            std::map<VarId, LossyFiniteAutomaton> variablesToConstants;
            systemBeingProcessed = LossyNonCommutativePolynomial::eliminateTerminalsInNonterminalProductions
                    (equations, variablesToConstants, eliminateTerminalsInLinearProductions);

            // binarize all productions, i.e. nonterminal productions of length at least 3 will be split up until
            // there are only productions of length <= 2 left
            systemBeingProcessed = LossyNonCommutativePolynomial::binarizeProductions(systemBeingProcessed, variablesToConstants);

            return systemBeingProcessed;
    }

    /*
     * Lossifies a given QNF; this means that every variable occurring in a polynomial is added to it as
     */
    static void lossifyQNF(std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations) {

        // add the monomials of degree 1; since we are in an idempotent semiring, we have a+a = a and so
        // we don't need to worry about the number of occurrences of each variable
        std::set<VarId> vars;
        LossyNonCommutativePolynomial monomialsOfDegreeOne;
        for(auto &equation: equations) {
            vars = equation.second.get_variables_quadratic_monomials();
            monomialsOfDegreeOne = LossyNonCommutativePolynomial::null();

            for(VarId var: vars) {
                monomialsOfDegreeOne += LossyNonCommutativePolynomial(var);
            }

            equation.second = equation.second + monomialsOfDegreeOne;
        }

        // add 1 to each equation
        for(auto &equation: equations) {
            equation.second = equation.second + LossyNonCommutativePolynomial::one();
        }
    }

    /*
     * Checks if a given variable is productive in a given clean system.
     */
    static bool variableHasProduction(const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &cleanEquations,
            const VarId &S) {
        bool ShasProduction = false;

        for(auto &equation: cleanEquations) {
            if(equation.first == S) {
                ShasProduction = true;
            }
        }

        return ShasProduction;
    }

    /*
     * Evaluates a polynomial system at a given vector.
     *
     * Not passing valuation as a reference is purposeful, since we need to copy it anyway as we don't want to change it.
     */
    static std::map<VarId, LossyFiniteAutomaton> evaluateSystem (std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations,
        int &times, ValuationMap<LossyFiniteAutomaton> valuation) {

        ValuationMap<LossyFiniteAutomaton> tempValuation;
        VarId variable;

        // iterate the desired number of times
        for(int i = 0; i < times; i++) {
            time_t rawtime;
            time (&rawtime);

            // evaluate each polynomial and map the appropriate variable to the result
            for(auto &equation: equations) {
                tempValuation.insert(std::make_pair(equation.first, equation.second.eval(valuation)));
            }

            // prepare next iteration
            valuation.swap(tempValuation);

            if(i + 1 != times) {
                tempValuation.clear();
            }
        }

        return valuation;
    }

    /*
     * Does what it says.
     */
    static std::map<VarId, int> mapVariablesToComponents(
            std::vector<std::vector<std::pair<VarId, LossyNonCommutativePolynomial>>> &components) {

        std::map<VarId, int> varToComponent;

        for(int i = 0; i < components.size(); i++) {
            for(auto &equation: components[i]) {
                varToComponent.insert(std::make_pair(equation.first, i));
            }
        }

        return varToComponent;
    }


    /*
     * Finds components where we can duplicate literals, i.e. where we have derivations
     * of the form A ->* _A_A_; in the paper by Courcelle, those are the components where A <_2 A.
     */
    static std::set<int> findSquarableComponents(
            std::vector<std::vector<std::pair<VarId, LossyNonCommutativePolynomial>>> &components,
            std::map<VarId, int> &varToComponent) {

        std::set<int> squarableComponents;

        for(int i = 0; i < components.size(); i++) {
            bool squarable = false;

            for(auto &equation: components[i]) {
                squarable |= equation.second.componentIsSquarable(varToComponent, i);

                if(squarable) {
                    squarableComponents.insert(i);
                    break;
                }
            }
        }

        return squarableComponents;
    }

    /*
     * Find for each component the components which are reachable from it.
     */
    static std::map<int, std::set<int>> findReachableComponents(
            std::vector<std::vector<std::pair<VarId, LossyNonCommutativePolynomial>>> &components,
            std::map<VarId, int> &varToComponent) {
        std::map<int, std::set<int>> reachabilityMap;

        // iterate over all components
        for(int i = 0; i < components.size(); i++) {
            std::set<int> reachableComponents;
            reachableComponents.insert(i);

            // build the set of reachable components: for each variable in the component,
            // check which other variables it can produce and add their components to the set of
            // reachable components
            for(auto &equation: components[i]) {
                auto tmp = equation.second.get_variables();

                for(auto var: tmp) {

                    // don't re-add the stuff from within the current component
                    if(varToComponent[var] < i) {
                        auto tmpReachable = reachabilityMap[varToComponent[var]];
                        reachableComponents.insert(tmpReachable.begin(), tmpReachable.end());
                    }
                }
            }

            reachabilityMap.insert(std::make_pair(i, reachableComponents));
        }

        return reachabilityMap;
    }

    /*
     * Find the letters reachable from each component
     *
     * The function assumes that "components" is sorted in reverse topological order, i.e. the component
     * that depends on no other component comes first.
     */
    static std::map<int, std::set<unsigned char>> findReachableLetters(
            std::vector<std::vector<std::pair<VarId, LossyNonCommutativePolynomial>>> &components,
            std::map<VarId, int> &varToComponent) {

        std::map<int, std::set<unsigned char>> componentToReachableLetters;
        std::map<int, std::set<int>> reachableComponents = findReachableComponents(components, varToComponent);

        for (int i = 0; i < components.size(); i++) {
            std::set<unsigned char> reachableLetters;

            // iterate over all components reachable from component i
            for(auto reachableComp: reachableComponents[i]) {

                // for every lower component, we have already calculated the set of reachable letters;
                // just add it to the letters reachable from i
                if(reachableComp != i) {
                    auto letters = componentToReachableLetters[reachableComp];
                    reachableLetters.insert(letters.begin(), letters.end());
                } else { // add the letters appearing in some production that stays within i itself
                    for(auto &equation: components[i]) {
                        std::set<unsigned char> letters = equation.second.get_terminals();

                        reachableLetters.insert(letters.begin(), letters.end());
                    }
                }
            }

            // once we're done, remember the set of reachable letters
            componentToReachableLetters.insert(std::make_pair(i, reachableLetters));
        }

        return componentToReachableLetters;
    }

    static void calculateClosuresOfSquarableComponents(
            std::map<int, std::set<unsigned char>> &componentToReachableLetters,
            std::set<int> &squarableComponents,
            std::map<int, LossyFiniteAutomaton> &componentToClosure) {

        // for squarable components, we only need the set of letters that appear in the strings derivable from that component
        // and star that set; the set is always nonempty since otherwise, the variable would be unproductive and would have
        // been eliminated while cleaning the system
        for(auto i: squarableComponents) {
            std::set<unsigned char> reachableLetters = componentToReachableLetters[i];

            if(reachableLetters.empty()) {
                componentToClosure.insert(std::make_pair(i, LossyFiniteAutomaton::one()));
            } else {
                std::stringstream ss;
                ss << "[";

                for(unsigned char letter: reachableLetters) {
                    ss << letter;
                }

                ss << "]*";
                componentToClosure.insert(std::make_pair(i, LossyFiniteAutomaton(ss.str())));
            }
        }
    }

    /*
     * Finds monomials of degree two where both variables are in a different and therefore lower
     * scc than the axiom of the production.
     *
     * Assumes that "components" is sorted in reverse topological order.
     */
    static std::map<int, std::map<int, std::set<int>>>  mapQuadraticLHStoRHS(
            std::vector<std::vector<std::pair<VarId, LossyNonCommutativePolynomial>>> &components,
            std::map<VarId, int> &varToComponent,
            std::set<int> &squarableComponents) {

        std::map<int, std::map<int, std::set<int>>> quadraticLHStoRHS;

        for(int i = 0; i < components.size(); i++) {
            if(squarableComponents.count(i) == 0) {
                for(auto &equation: components[i]) {
                    equation.second.mapQuadraticLHStoRHS(quadraticLHStoRHS, varToComponent, i);
                }
            }
        }

        return quadraticLHStoRHS;
    }

    static void calculateLowerComponentVariables(
            std::map<int, std::set<int>> &lhsLowerComponentVariables,
            std::map<int, std::set<int>> &rhsLowerComponentVariables,
            std::vector<std::vector<std::pair<VarId, LossyNonCommutativePolynomial>>> &components,
            std::map<VarId, int> &varToComponent,
            std::set<int> &squarableComponents) {

        for(int i = 0; i < components.size(); i++) {
            if(squarableComponents.count(i) == 0) {
                for(auto &equation: components[i]) {
                    equation.second.calculateLowerComponentVariables(lhsLowerComponentVariables, rhsLowerComponentVariables,
                            components, varToComponent, i);
                }
            }
        }
    }

    static void calculateSameComponentLetters(
            std::map<int, std::set<unsigned char>> &lhsSameComponentLetters,
            std::map<int, std::set<unsigned char>> &rhsSameComponentLetters,
            std::vector<std::vector<std::pair<VarId, LossyNonCommutativePolynomial>>> &components,
            std::map<VarId, int> &varToComponent,
            std::set<int> &squarableComponents) {

        for(int i = 0; i < components.size(); i++) {
            if(squarableComponents.count(i) == 0) {
                for(auto &equation: components[i]) {
                    equation.second.calculateSameComponentLetters(lhsSameComponentLetters, rhsSameComponentLetters,
                            components, varToComponent, i);
                }
            }
        }
    }

    static LossyFiniteAutomaton compareClosures(const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations_1,
            const VarId &S_1, const std::vector<std::pair<VarId, LossyNonCommutativePolynomial>> &equations_2,
            const VarId &S_2, const LossyFiniteAutomaton approx_1, const LossyFiniteAutomaton approx_2) {

        bool A1_subset_A2 = approx_2.contains(approx_1);
        bool A2_subset_A1 = approx_1.contains(approx_2);

        if(A1_subset_A2 && A2_subset_A1) {
            return LossyFiniteAutomaton::null();
        } else {
            LossyFiniteAutomaton L1_intersect_A2c = LossyFiniteAutomaton::null();
            bool L1_intersect_A2c_changed = false;
            LossyFiniteAutomaton L2_intersect_A1c = LossyFiniteAutomaton::null();
            bool L2_intersect_A1c_changed = false;

            if(!A1_subset_A2) {
                VarId startSymbol_1_2;
                auto A2c = approx_1.minus(approx_2);
                auto intersectionGrammar = A2c.intersectionWithCFG(startSymbol_1_2, S_1, equations_1);

                std::queue<VarId> worklist;
                worklist.push(startSymbol_1_2);
                intersectionGrammar = LossyNonCommutativePolynomial::cleanSystem(intersectionGrammar, worklist);

                L1_intersect_A2c = LossyNonCommutativePolynomial::shortestWord
                        (intersectionGrammar, startSymbol_1_2);
                L1_intersect_A2c_changed = true;
            }

            if(!A2_subset_A1) {
                VarId startSymbol_2_1;
                auto A1c = approx_2.minus(approx_1);
                auto intersectionGrammar_2 = A1c.intersectionWithCFG(startSymbol_2_1, S_2, equations_2);

                std::queue<VarId> worklist;
                worklist.push(startSymbol_2_1);
                intersectionGrammar_2 = LossyNonCommutativePolynomial::cleanSystem(intersectionGrammar_2, worklist);

                L2_intersect_A1c = LossyNonCommutativePolynomial::shortestWord
                        (intersectionGrammar_2, startSymbol_2_1);
                L2_intersect_A1c_changed = true;
            }

            if(L1_intersect_A2c_changed && L2_intersect_A1c_changed) {
                if(L1_intersect_A2c.size() <= L2_intersect_A1c.size() && L1_intersect_A2c != LossyFiniteAutomaton::null()) {
                    return L1_intersect_A2c;
                } else {
                    return L2_intersect_A1c;
                }
            } else {
                if(L1_intersect_A2c_changed) {
                    return L1_intersect_A2c;
                } else {
                    return L2_intersect_A1c;
                }
            }
        }
    }
};
