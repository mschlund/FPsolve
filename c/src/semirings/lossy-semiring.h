#ifndef LOSSY_SEMIRING_H_
#define LOSSY_SEMIRING_H_

#include <string>
#include <memory>
#include <unordered_map>
#include <queue>
#include <time.h>
#include <stdint.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <algorithm>



#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graphviz.hpp>

#include "../datastructs/hash.h"
#include "../datastructs/matrix.h"
#include "../datastructs/var.h"
#include "../datastructs/free-structure.h"
#include "../polynomials/non_commutative_polynomial.h"
#include "../polynomials/non_commutative_monomial.h"
#include "../utils/timer.h"

#include "semiring.h"

template<typename LSR>
class Evaluator;

template <typename LSR>
struct Vertex {
    VarId var;         // var and rex combines the equations in the vertex
    NonCommutativePolynomial<LSR> rex;
};

template <typename LSR>
class LossySemiring: public StarableSemiring<LSR,
        Commutativity::NonCommutative, Idempotence::Idempotent> {
public:
    virtual ~LossySemiring(){};
    virtual LSR lossify() const = 0;

    /*
     * Solves a polynomial system over the lossy semiring. For derivation of
     * the algorithm and proof of correctness, see "Esparza, Kiefer, Luttenberger:
     * Derivation Tree Analysis for Accelerated Fixed-Point Computation".
     *
     * If the system is not clean, i.e. there are polynomials in it with a fixpoint of 0, then you must set
     * systemNeedsCleaning to true in order for the algorithm to work as intended.
     *
     * If the system is not in lossy QNF, then you need to set systemNeedsLossyQNF to true.
     */
    static std::map<VarId, LSR> solvePolynomialSystem(
            std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations,
            bool systemNeedsCleaning,
            bool systemNeedsLossyQNF) {

        std::map<VarId, LSR> solution;

//        auto components = group_by_scc(equations);
//        for(int i = 0; i < components.size(); i++) {
//            std::cout << "component " << i << ":\t" << std::endl;
//
//            for(auto &equation: components[i]) {
//                std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//            }
//        }

        // start out with a clean system
        if(systemNeedsCleaning) {
            std::queue<VarId> worklist;
            equations = cleanSystem(equations, worklist);
        }

        // bring the clean system into qnf, if necessary
        std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> qnf;

        if(systemNeedsLossyQNF) {
            qnf = quadraticNormalForm(equations);
            lossifyQNF(qnf);
        } else {
            qnf = equations;
        }

        // build a zero vector for the evaluation of f
        std::map<VarId, LSR> zeroSystem;
        for(auto &equation: qnf) {
            zeroSystem.insert(std::make_pair(equation.first, LSR::null()));
        }

        // build the vectors f(0) and f^n(0) where n is the number of equations in the system
        int times = 1;
//        std::cout << "f(0), 1 iteration:" << std::endl;
        std::map<VarId, LSR> f_0 = evaluateSystem(qnf, times, zeroSystem);
        times = qnf.size() - 1;
//        std::cout << "f^n(0), "<< times << " iterations total:" << std::endl;
        std::map<VarId, LSR> f_n_0 = evaluateSystem(qnf, times, f_0);

        // find the LossySemiring element in the "middle" of the expression
        LSR middle = LSR::null();
        for(auto &elem_mapping: f_0) {
            middle += elem_mapping.second;
        }

        LSR::maxStates = std::max(LSR::maxStates,middle.size());

        // build the differential of the system
        std::map<VarId, NonCommutativePolynomial<LSR>> differential;
        for(auto &equation: qnf) {
                differential.insert(std::make_pair(equation.first, equation.second.differential_at(f_n_0)));
        }

//                std::cout << "differential:" << std::endl;
//                for(auto &mapping: differential) {
//                    std::cout << Var::GetVar(mapping.first).string() + " -> " + mapping.second.string() << std::endl;
//                }
//
//                std::cout << "f(0):" << std::endl;
//                for(auto &mapping: f_0) {
//                    std::cout << Var::GetVar(mapping.first).string() + ":\t\t" + mapping.second.string() << std::endl;
//                }
//                std::cout << "f^n(0):" << std::endl;
//                for(auto &mapping: f_n_0) {
//                    std::cout << Var::GetVar(mapping.first).string() + ":\t\t" + mapping.second.string() << std::endl;
//                }

        // sum all polynomials of the differential of the system; we don't need attribution
        // of each polynomial to the respective variable since we only care about the
        // leading and trailing coefficients of each monomial
        NonCommutativePolynomial<LSR> differential_sum = NonCommutativePolynomial<LSR>::null();
        for(auto &equation: differential) {
            differential_sum += equation.second;
        }

//        // add the constant summands appearing in any component of the differential
//        middle += differential_sum.eval(zeroSystem);
//                            std::cout << "middle: " << middle.string() << std::endl;

        // get the lefthand and righthand semiring element of the fixpoint
        LSR lefthandSum = differential_sum.getSumOfLeadingFactors();

        LSR::maxStates = std::max(LSR::maxStates,lefthandSum.size());
//                            std::cout << "LHS: " << lefthandSum.string() << std::endl;
        LSR righthandSum = differential_sum.getSumOfTrailingFactors();
        LSR::maxStates = std::max(LSR::maxStates,righthandSum.size());
//                            std::cout << "RHS: " << righthandSum.string() << std::endl;
        // fixpoint element
        LSR fixpoint = lefthandSum.star() * middle * righthandSum.star();
        LSR::maxStates = std::max(LSR::maxStates,fixpoint.size());

        for(auto &variable_mapping: qnf) {
            solution.insert(std::make_pair(variable_mapping.first, fixpoint));
        }

        return solution;
    }

    /*
     * Takes two grammars and refines them by intersecting them with the prefix languages of their shared downward closure.
     */
    static LSR refineCourcelle(const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations_1,
            const VarId &S_1, const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations_2,
            const VarId &S_2, int maxLengthOfPrefixes) {
        int who = RUSAGE_SELF;
        struct rusage usage;
        int ret;
        struct timeval startT, endT;
        long mtime, seconds, useconds;
        gettimeofday(&startT, NULL);
        Timer timer;
        timer.Start();
        int largestRawGrammar = 0;
        int largestCleanGrammar = 0;
        int iterationsBeforeSuccess = 0;
        int outputLength = 0;

        LSR approx_1 = downwardClosureCourcelle(equations_1, S_1);
        LSR approx_2 = downwardClosureCourcelle(equations_2, S_2);

//        std::cout << "approx 1: " << approx_1.string() << std::endl;
//        std::cout << "approx 2: " << approx_1.string() << std::endl;
        bool A1_subset_A2 = approx_2.contains(approx_1);
        bool A2_subset_A1 = approx_1.contains(approx_2);
        ret = getrusage(who, &usage);
        auto memoryUsage = usage.ru_maxrss;
//        std::cout << memoryUsage << ",";

        if(A1_subset_A2 && A2_subset_A1) {
            LSR difference = LSR::null();

            if(!(approx_1 == LSR::null()) && !(approx_1 == LSR::one())) {

                // build the automaton Sigma*
                std::string regexAlphabetStar = "[";
                for (auto c: approx_1.alphabet()) {
                    regexAlphabetStar.append(1, c);
                }
                regexAlphabetStar += "]*";

                // build the map of (length, prefixes of that length)
                std::map<int, std::set<std::string>> prefixesPerLength = approx_1.prefixesToMaxLength(maxLengthOfPrefixes);

                // iterate over the prefixes
                bool noDifference = true;
                std::queue<VarId> worklist;
                for(int i = 1; (i <= maxLengthOfPrefixes) && noDifference; i++) {
                    iterationsBeforeSuccess++;
                    for(auto &prefix: prefixesPerLength[i]) {
                        std::string prefixAlphaStar = std::string(prefix) + regexAlphabetStar;
//                        std::cout << "prefixalphastar\t" << prefixAlphaStar << std::endl;
                        LSR prefixAutomaton(prefixAlphaStar);
//                        std::cout << "done prefixalphastar" << std::endl;
//                        std::cout << "prefix automaton:\t" << prefixAutomaton.string() << std::endl;

                        VarId S_1_partition, S_2_partition;

                        // generate the subset of language 1
                        auto equations_1_Partition = prefixAutomaton.intersectionWithCFG(S_1_partition, S_1, equations_1);
                        largestRawGrammar = equations_1_Partition.size() > largestRawGrammar ? equations_1_Partition.size() : largestRawGrammar;
                        while(!worklist.empty()) {
                            worklist.pop();
                        }
                        worklist.push(S_1_partition);
                        equations_1_Partition = cleanSystem(equations_1_Partition, worklist);
                        largestCleanGrammar = equations_1_Partition.size() > largestCleanGrammar ? equations_1_Partition.size() : largestCleanGrammar;

                        // output grammar 1
//                        std::cout << "equ. 1 part:" << std::endl;
//                        for(auto &equation: equations_1_Partition) {
//                            std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//                        }

                        // generate the subset of language 2
                        auto equations_2_Partition = prefixAutomaton.intersectionWithCFG(S_2_partition, S_2, equations_2);
                        largestRawGrammar = equations_2_Partition.size() > largestRawGrammar ? equations_2_Partition.size() : largestRawGrammar;
                        while(!worklist.empty()) {
                            worklist.pop();
                        }
                        worklist.push(S_2_partition);
                        equations_2_Partition = cleanSystem(equations_2_Partition, worklist);
                        largestCleanGrammar = equations_2_Partition.size() > largestCleanGrammar ? equations_2_Partition.size() : largestCleanGrammar;

                        // output grammar 2
//                        std::cout << "equ. 2 part:" << std::endl;
//                        for(auto &equation: equations_2_Partition) {
//                            std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//                        }

                        // approximate the grammars
                        LSR approx_1_part = downwardClosureCourcelle(equations_1_Partition, S_1_partition);
                        LSR approx_2_part = downwardClosureCourcelle(equations_2_Partition, S_2_partition);

                        LSR tempDiff = compareClosures(equations_1_Partition, S_1_partition, equations_2_Partition, S_2_partition, approx_1_part, approx_2_part, largestRawGrammar, largestCleanGrammar);

                        if(!(tempDiff == LSR::null())) {
                            difference = tempDiff;
                            noDifference = false;
                            break;
                        }
                    }
                }
            }


            timer.Stop();
            ret = getrusage(who, &usage);
            auto memoryUsage = usage.ru_maxrss;
            if(difference != LSR::null()) {
                outputLength = difference.string().size();
            }
            std::cout << timer.GetMilliseconds().count() << "," <<  memoryUsage << "," << largestRawGrammar << "," << largestCleanGrammar << "," << iterationsBeforeSuccess << "," << outputLength << ",";
            return difference;
        } else {
            LSR difference = compareClosures(equations_1, S_1, equations_2, S_2, approx_1, approx_2, largestRawGrammar, largestCleanGrammar);
            timer.Stop();
            ret = getrusage(who, &usage);
            auto memoryUsage = usage.ru_maxrss;
            if(difference != LSR::null()) {
                outputLength = difference.string().size();
            }
            std::cout << timer.GetMilliseconds().count() << "," << memoryUsage << "," << largestRawGrammar << "," << largestCleanGrammar << "," << iterationsBeforeSuccess << "," << outputLength << ",";
            return difference;
        }
    }

    /*
     * Approximates the language defined by the grammar starting at S given in the system of equations by calculating
     * its downward closure.
     */
    static LSR downwardClosureDerivationTrees(const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations,
            const VarId &S) {

        int who = RUSAGE_SELF;
        struct rusage usage;
        int ret;
        struct timeval startT, endT;
        long mtime, seconds, useconds;
        gettimeofday(&startT, NULL);
        Timer timer;
        timer.Start();

        std::queue<VarId> worklist;
        worklist.push(S);

//        std::cout << "system:" << std::endl;
//        for(auto &equation: equations) {
//            std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
//        }

        auto cleanEquations = cleanSystem(equations, worklist);

        // if the start symbol is unproductive, then it doesn't generate anything;
        // the downward closure of the empty set is the empty set
        if(!variableHasProduction(cleanEquations, S)) {

            ret = getrusage(who, &usage);
            gettimeofday(&endT, NULL);seconds  = endT.tv_sec  - startT.tv_sec;
            useconds = endT.tv_usec - startT.tv_usec;
            mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
            auto memoryUsage = usage.ru_maxrss;
            timer.Stop();
//            std::cout << /*"\t time:\t" << */mtime << ","/*<< "\ttimer time: " << timer.GetMilliseconds().count()*/ /* << "\tmemory used: "*/ << usage.ru_maxrss <<","/*<< "KB\t" << "max number of states: " */<< LSR::maxStates<<","/* << "\tstates before lossification: " */<< "0" << std::endl;
            return LSR::null();
        }

        auto cleanQNF = quadraticNormalForm(cleanEquations);
        lossifyQNF(cleanQNF);
        auto components = group_by_scc(cleanQNF);

//        std::cout << "derivation trees components" << std::endl;
//        for(int i = 0; i < components.size(); i++) {
//            std::cout << "component: " << i << std::endl;
//            for(auto &equation: components[i]) {
//                std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//            }
//        }

        std::map<VarId, LSR> knownValuations;
        for(int i = 0; i < components.size(); i++) {
//            std::cout << "component: " << i << "\t size: " << components[i].size() << std::endl;
            // if this is not the bottom component, substitute known valuations
            if(i != 0) {
                for(auto &equation: components[i]) {
                    equation.second = equation.second.partial_eval(knownValuations);
                }
            }

            auto newValuations = solvePolynomialSystem(components[i], false, false);

            for(auto &value: newValuations) {
                knownValuations.insert(std::make_pair(value.first, value.second));
                LSR::maxStates = std::max(LSR::maxStates,knownValuations[value.first].size());
            }
        }

        ret = getrusage(who, &usage);
        gettimeofday(&endT, NULL);seconds  = endT.tv_sec  - startT.tv_sec;
        useconds = endT.tv_usec - startT.tv_usec;
        mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
        auto memoryUsage = usage.ru_maxrss;
        timer.Stop();
//        std::cout << /*"\t time:\t" << */mtime << ","/*<< "\ttimer time: " << timer.GetMilliseconds().count()*/ /* << "\tmemory used: "*/ << usage.ru_maxrss <<","/*<< "KB\t" << "max number of states: " */<< LSR::maxStates<<","/* << "\tstates before lossification: " */<< knownValuations[S].size() << std::endl;
//        std::cout << "closure before lossification: " << knownValuations[S].string() << "\t";
//        std::cout << "size of closure automaton before lossification:\t" << knownValuations[S].size() << std::endl;
        return knownValuations[S].lossify();
    }

    /*
     * Calculates the downward closure of the grammar defined by "equations" starting at "S"; the algorithm
     * derives from the "On Constructing Obstruction Sets of Words", B. Courcelle, Bulletin of EATCS 1991.
     *
     * WARNING: This function will break if you use an LSR that does not have a constructor that takes POSIX regex string.
     */
    static LSR downwardClosureCourcelle(const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations,
            const VarId &S) {

        int who = RUSAGE_SELF;
        struct rusage usage;
        int ret;
        struct timeval startT, endT;
        long mtime, seconds, useconds;
        gettimeofday(&startT, NULL);

        std::queue<VarId> worklist;
        worklist.push(S);

//        std::cout << "system:" << std::endl;
//        for(auto &equation: equations) {
//            std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
//        }

        auto cleanEquations = cleanSystem(equations, worklist);

        // if the start symbol is unproductive, then it doesn't generate anything;
        // the downward closure of the empty set is the empty set
        if(!variableHasProduction(cleanEquations, S)) {
            ret = getrusage(who, &usage);
            gettimeofday(&endT, NULL);seconds  = endT.tv_sec  - startT.tv_sec;
            useconds = endT.tv_usec - startT.tv_usec;
            mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
            auto memoryUsage = usage.ru_maxrss;
//            std::cout << /*"\t time:\t" << */mtime << ","/*<< "\ttimer time: " << timer.GetMilliseconds().count()*/ /* << "\tmemory used: "*/ << usage.ru_maxrss <<","/*<< "KB\t" << "max number of states: " */<< LSR::maxStates<<","/* << "\tstates before lossification: " */<< "0";
            return LSR::null();
        }

        auto cleanQNF = quadraticNormalForm(cleanEquations);
        auto components = group_by_scc(cleanQNF);

//        for(int i = 0; i < components.size(); i++) {
//            std::cout << "component: " << i << std::endl;
//            for(auto &equation: components[i]) {
//                std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//            }
//        }

        // will hold the downward closure of all components
        std::map<int, LSR> componentToClosure;

        // map variables to their components
        std::map<VarId, int> varToComponent = mapVariablesToComponents(components);

//        std::cout << "mapped variables to their components" << std::endl;

        // these components will have as their closure simply the star of the letters reachable
        // from the variables of the component
        std::set<int> squarableComponents = findSquarableComponents(components, varToComponent);

//        std::cout << "found squarable components" << std::endl;

        // we need this map for the case where we can duplicate a variable
        std::map<int, std::set<unsigned char>> componentToReachableLetters =
                findReachableLetters(components, varToComponent);

//        std::cout << "mapped components to reachable letters" << std::endl;

        // we can already calculate the downward closures of squarable components
        calculateClosuresOfSquarableComponents(componentToReachableLetters, squarableComponents, componentToClosure);

//        std::cout << "closures of squarable components done" << std::endl;

        // this map holds for each component all pairs AB where A and B are in a lower component;
        // we use the downward closure for those pairs and concatenate accordingly to construct the
        // "middle" part of the downward closure the way Courcelle calculates it
        std::map<int, std::map<int, std::set<int>>> quadraticLHStoRHS = mapQuadraticLHStoRHS(components, varToComponent, squarableComponents);
//        std::cout << "quadraticLHStoRHS[4][3] contains 3: " << quadraticLHStoRHS[4][3].count(3) << std::endl;
//        std::cout << "quadraticLHStoRHS[5][1] contains 4: " << quadraticLHStoRHS[5][1].count(4) << std::endl;
//        std::cout << "quadraticLHStoRHS is not null, size: " << quadraticLHStoRHS.size() << std::endl;
//        for(auto &entry: quadraticLHStoRHS) {
//            std::cout << "quadraticLHStoRHS component: " << entry.first << std::endl;

//            for(auto &mapping: entry.second) {
//                for(auto target: mapping.second) {
//                    std::cout << "mapping: " << mapping.first << " to " << target << std::endl;
//                }
//            }
//        }



//        std::cout << "mapped quadratic lhs to rhs" << std::endl;

        // get the components of variables Y_l and Y_r such that there is a monomial in the component of any nonsquarable B
        // that has the form Y_l*B or B*Y_r;
        // their closures will be part of the lefthand/righthand terms of the calculation by Courcelle
        std::map<int, std::set<int>> lhsLowerComponentVariables;
        std::map<int, std::set<int>> rhsLowerComponentVariables;
        calculateLowerComponentVariables(lhsLowerComponentVariables, rhsLowerComponentVariables,
                components, varToComponent, squarableComponents);

//        std::cout << "got lower component variables" << std::endl;

        // get the lefthand/righthand letters for each nonsquarable component, i.e. the letters occurring in terms xBy
        // where x, y are strings over the terminal alphabet and B is a variable from the same component as the current one
        std::map<int, std::set<unsigned char>> lhsSameComponentLetters;
        std::map<int, std::set<unsigned char>> rhsSameComponentLetters;
        calculateSameComponentLetters(lhsSameComponentLetters, rhsSameComponentLetters,
                components, varToComponent, squarableComponents);

//        std::cout << "got same component letters" << std::endl;

        // for the middle part of each component, we also need the sum of the closures of all constant monomials
        // appearing in that component
        std::map<int, LSR> closuresOfConstantMonomials;
        calculateClosuresOfConstantMonomials(closuresOfConstantMonomials, components, squarableComponents);

//        std::cout << "closure of constant monomials done" << std::endl;

        // we have everything we can prepare beforehand; everything else will need to be dealt with while
        // putting the closures together
        calculateClosuresOfNonsquarableComponents(componentToClosure, components,
                squarableComponents, quadraticLHStoRHS,
                lhsLowerComponentVariables, rhsLowerComponentVariables,
                lhsSameComponentLetters, rhsSameComponentLetters,
                closuresOfConstantMonomials, varToComponent, componentToReachableLetters);

//        std::cout << "closures of nonsquarable components done" << std::endl;

        ret = getrusage(who, &usage);
        gettimeofday(&endT, NULL);seconds  = endT.tv_sec  - startT.tv_sec;
        useconds = endT.tv_usec - startT.tv_usec;
        mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
        auto memoryUsage = usage.ru_maxrss;
//        std::cout  /*<< "\t time:\t" */<< mtime <<"," /*<< "\tmemory used: " */<< usage.ru_maxrss<<"," /*<< "KB\t" << "max number of states: " */<< LSR::maxStates <<","/*<< "\tstates of closure: " */<< componentToClosure[varToComponent[S]].size();

        return componentToClosure[varToComponent[S]];
    }

private:

    /*
     * Cleans the polynomial system, i.e. removes variables that are unproductive for all the grammars we get
     * if exactly the symbols in "startSymbols" are considered as initial symbols of a grammar.
     *
     * In the context of fixpoints, if you only want to eliminate symbols that have a 0 fixpoint, hand this function an
     * empty queue "startSymbols".
     */
    static std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> cleanSystem
        (const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations, std::queue<VarId> &startSymbols) {

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
        return NonCommutativePolynomial<LSR>::cleanSystem(equations, startSymbols);
    }

    /*
     * Brings a system into quadratic normal form.
     */
    static std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> quadraticNormalForm
            (const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations) {
            std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> systemBeingProcessed;

            // for all nonterminal productions, introduce new variables for the terminal factors that
            // appear in them as one would when calculating the CNF of a CFG; do not do this for linear terms since
            // QNF allows for linear terms
            std::map<VarId, LSR> variablesToConstants;
            systemBeingProcessed = NonCommutativePolynomial<LSR>::eliminateTerminalsInNonterminalProductions
                    (equations, variablesToConstants, false);

//            std::cout << "system after constant elimination:" << std::endl;
//            for(auto &equation: systemBeingProcessed) {
//                std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//            }

            // binarize all productions, i.e. nonterminal productions of length at least 3 will be split up until
            // there are only productions of length <= 2 left
            systemBeingProcessed = NonCommutativePolynomial<LSR>::binarizeProductions(systemBeingProcessed, variablesToConstants);

//            std::cout << "system after binarization:" << std::endl;
//            for(auto &equation: systemBeingProcessed) {
//                std::cout << Var::GetVar(equation.first).string() << " -> " << equation.second.string() << std::endl;
//            }

            return systemBeingProcessed;
    }

    /*
     * Lossifies a given QNF; this means that every variable occurring in a polynomial is added to it as
     */
    static void lossifyQNF(std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations) {

        // add the monomials of degree 1; since we are in an idempotent semiring, we have a+a = a and so
        // we don't need to worry about the number of occurrences of each variable
        std::set<VarId> vars;
        NonCommutativePolynomial<LSR> monomialsOfDegreeOne;
        for(auto &equation: equations) {
            vars = equation.second.get_variables_quadratic_monomials();
            monomialsOfDegreeOne = NonCommutativePolynomial<LSR>::null();

            for(VarId var: vars) {
                monomialsOfDegreeOne += NonCommutativePolynomial<LSR>(var);
            }

            equation.second = equation.second + monomialsOfDegreeOne;
        }

        // add 1 to each equation
        for(auto &equation: equations) {
            equation.second = equation.second + NonCommutativePolynomial<LSR>::one();
        }

//        std::cout << "lossy QNF:" << std::endl;
//        for(auto &equation: equations) {
//            std::cout << Var::GetVar(equation.first).string() + " -> " + equation.second.string() << std::endl;
//        }
    }

    /*
     * Used to decompose a system into its strongly connected components.
     *
     * The resulting vector will have the components bottom-up, i.e. if the graph is {A,B,C} with edges A->B, B->C, C->B,
     * then the vector will have component {B,C} at its first entry, the second entry will be {A}. Components between
     * which there is no directed path may appear in any order.
     */
    static std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> group_by_scc
        (const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations) {

        // create map of variables to [0..n]. this is used to enumerate important variables
        // in a clean way from 0 to n during graph construction
        std::unordered_map<VarId, int> var_key;

        // build the graph
        boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, Vertex<LSR>> graph(equations.size());

        for (const auto &eq : equations) {
            // if the variable is not yet in the map insert it together with the size of
            // the map. this way we get a unique identifier for the variables counting from 0 to n
            if (var_key.find(eq.first) == var_key.end()) {
                var_key.insert(std::make_pair(eq.first, var_key.size()));
            }

            int a = var_key.find(eq.first)->second; // variable key
            graph[a].var = eq.first; // store VarId to the vertex
            graph[a].rex = eq.second; // store the regular expression to the vertex

            auto variables = eq.second.get_variables(); // all variables of this rule;

            for (const auto &var : variables) {
                if (var_key.find(var) == var_key.end()) { // variable is not yet in the map
                    var_key.insert(var_key.begin(), std::pair<VarId,int>(var,var_key.size()));
                }

                int b = var_key.find(var)->second; // variable key
                boost::add_edge(a, b, graph);
            }
        }

        // calculate strongly connected components and store them in 'component'
        std::vector<int> component(boost::num_vertices(graph));
        boost::strong_components(graph,&component[0]);

        // group neccessary equations together
        int num_comp = *std::max_element(component.begin(), component.end()) + 1; // find the number of components
        std::vector<std::vector<std::pair<VarId,NonCommutativePolynomial<LSR>>>> grouped_equations(num_comp);

        // iterate over all vertices (0 to n)
        // collect the necessary variables + equations for every component
        for (std::size_t j = 0; j != component.size(); ++j) {
            grouped_equations[component[j]].push_back(std::pair<VarId,NonCommutativePolynomial<LSR>>(graph[j].var, graph[j].rex));
        }

        return grouped_equations;
    }

    /*
     * Checks if a given variable is productive in a given clean system.
     */
    static bool variableHasProduction(const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &cleanEquations,
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
    static std::map<VarId, LSR> evaluateSystem (std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations,
        int &times, std::map<VarId, LSR> valuation) {

        std::map<VarId, LSR> tempValuation;
        VarId variable;

        // iterate the desired number of times
        for(int i = 0; i < times; i++) {
            time_t rawtime;
            time (&rawtime);
//            std::cout << "system evaluation: iteration\t" << (i+1) << "\t at\t" << ctime (&rawtime);
            // evaluate each polynomial and map the appropriate variable to the result
            for(auto &equation: equations) {
                tempValuation.insert(std::make_pair(equation.first, equation.second.eval(valuation)));
                LSR::maxStates = std::max(LSR::maxStates,tempValuation[equation.first].size());
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
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> &components) {

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
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> &components,
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
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> &components,
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
     * Find the letters reachable from any variables that can be duplicated (in the paper
     * by Courcelle: for all A such that A <_2 A). We need those sets to immediately generate the
     * downward closures of those components.
     *
     * The function assumes that "components" is sorted in reverse topological order.
     */
    static std::map<int, std::set<unsigned char>> findReachableLetters(
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> &components,
            std::map<VarId, int> &varToComponent) {

        std::map<int, std::set<unsigned char>> componentToReachableLetters;
        std::map<int, std::set<int>> reachability = findReachableComponents(components, varToComponent);

        for (int i = 0; i < components.size(); i++) {
            std::set<unsigned char> reachableLetters;

            for(auto reachableComp: reachability[i]) {

                if(reachableComp != i) {
                    auto letters = componentToReachableLetters[reachableComp];
                    reachableLetters.insert(letters.begin(), letters.end());
                } else {
                    for(auto &equation: components[i]) {
                        auto letters = equation.second.get_terminals();
                        reachableLetters.insert(letters.begin(), letters.end());
                    }
                }
            }

            componentToReachableLetters.insert(std::make_pair(i, reachableLetters));
        }

        return componentToReachableLetters;
    }

    static void calculateClosuresOfSquarableComponents(
            std::map<int, std::set<unsigned char>> &componentToReachableLetters,
            std::set<int> &squarableComponents,
            std::map<int, LSR> &componentToClosure) {

        // for squarable components, we only need the set of letters that appear in the strings derivable from that component
        // and star that set; the set is always nonempty since otherwise, the variable would be unproductive and would have
        // been eliminated while cleaning the system
        for(auto i: squarableComponents) {
            std::set<unsigned char> reachableLetters = componentToReachableLetters[i];
            std::stringstream ss;
            ss << "[";

            for(unsigned char letter: reachableLetters) {
                ss << letter;
            }

            ss << "]*";
            componentToClosure.insert(std::make_pair(i, LSR(ss.str())));
            LSR::maxStates = std::max(LSR::maxStates,componentToClosure[i].size());
//            std::cout << "reachable letters for component " << i << ": " << LSR(ss.str()).string() << std::endl;
        }
    }

    /*
     * Finds monomials of degree two where both variables are in a different scc than the axiom of the production.
     */
    static std::map<int, std::map<int, std::set<int>>>  mapQuadraticLHStoRHS(
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> &components,
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
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> &components,
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
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> &components,
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

    static void calculateClosuresOfConstantMonomials(
            std::map<int, LSR> &closuresOfConstantMonomials,
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> &components,
            std::set<int> &squarableComponents) {

        for(int i = 0; i < components.size(); i++) {
            if(squarableComponents.count(i) == 0) {
                LSR closure = LSR::null();

                for(auto &equation: components[i]) {
                    closure = closure + equation.second.sumOfConstantMonomials();
                }


                LSR::maxStates = std::max(LSR::maxStates,closuresOfConstantMonomials[i].size());
                closuresOfConstantMonomials[i] = closure.lossify();
//                std::cout << "closure of constant monomials component " << i << ": " << closure.string() << std::endl;
            }
        }
    }

    static void calculateClosuresOfNonsquarableComponents(
            std::map<int, LSR> &componentToClosure,
            std::vector<std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>>> &components,
            std::set<int> &squarableComponents,
            std::map<int, std::map<int, std::set<int>>> &quadraticLHStoRHS,
            std::map<int, std::set<int>> &lhsLowerComponentVariables,
            std::map<int, std::set<int>> &rhsLowerComponentVariables,
            std::map<int, std::set<unsigned char>> &lhsSameComponentLetters,
            std::map<int, std::set<unsigned char>> &rhsSameComponentLetters,
            std::map<int, LSR> &closuresOfConstantMonomials,
            std::map<VarId, int> &varToComponent,
            std::map<int, std::set<unsigned char>> &componentToReachableLetters) {

        for(int i = 0; i < components.size(); i++) {
            if(squarableComponents.count(i) == 0){
//                std::cout << "component " << i << " of " << components.size() << std::endl;
                /*
                 * closures of elements in the middle
                 */

                //linear terms with variables in lower components
                std::set<NonCommutativeMonomial<LSR>> lowerLinearTerms;
                for(auto &equation: components[i]) {
                    equation.second.findLowerLinearTerms(lowerLinearTerms, varToComponent, i);
                }

//                std::cout << "found lower linear terms; total of " << lowerLinearTerms.size() << ". Calculating their closure" << std::endl;
                LSR closureOfLowerLinearTerms = findClosureOfLowerLinearTerms(lowerLinearTerms, varToComponent, componentToClosure);

                LSR::maxStates = std::max(LSR::maxStates,closureOfLowerLinearTerms.size());
//                std::cout << "closure of lower linear terms done" << std::endl;

                // quadratic terms with variables in lower components
                LSR closureOfLowerQuadraticTerms = LSR::null();
                for(auto &entry: quadraticLHStoRHS[i]) {
                    for(auto rhs: entry.second) {
//                        std::cout << entry.first << " * " << rhs << " in quadratic stuff of component " << i << std::endl;
                        closureOfLowerQuadraticTerms += componentToClosure[entry.first] * componentToClosure[rhs];
                        LSR::maxStates = std::max(LSR::maxStates,closureOfLowerQuadraticTerms.size());
                    }
                }



//                std::cout << "closure of lower quadratic terms done" << std::endl;

                // sum of the closures of: constant monomials in this component, linear and quadratic monomials with variables
                // exclusively in lower components
                LSR middle = closureOfLowerQuadraticTerms + closureOfLowerLinearTerms + closuresOfConstantMonomials[i];

                LSR::maxStates = std::max(LSR::maxStates,middle.size());
//                std::cout << "middle calculated" << std::endl;

                /*
                 * lefthand and righthand side reachable alphabets starred
                 */

                // find the reachable letters
                std::set<unsigned char> lefthandReachableLetters;
                lefthandReachableLetters.insert(lhsSameComponentLetters[i].begin(), lhsSameComponentLetters[i].end());

                for(auto variable: lhsLowerComponentVariables[i]) {
                    lefthandReachableLetters.insert(componentToReachableLetters[variable].begin(),
                            componentToReachableLetters[variable].end());
                }

//                std::cout << "lefthand reachable letters done" << std::endl;

                std::set<unsigned char> righthandReachableLetters;
                righthandReachableLetters.insert(rhsSameComponentLetters[i].begin(), rhsSameComponentLetters[i].end());

                for(auto variable: rhsLowerComponentVariables[i]) {
                    righthandReachableLetters.insert(componentToReachableLetters[variable].begin(),
                            componentToReachableLetters[variable].end());
                }

//                std::cout << "righthand reachable letters done" << std::endl;

                // build the automata
                LSR lefthandAlphabetStar = LSR::one();
                if(lefthandReachableLetters.size() != 0) {
                    std::stringstream ssLHS;
                    ssLHS << "[";

                    for(auto letter: lefthandReachableLetters) {
                        if(isalnum(letter)) {
                            ssLHS << letter;
                        }
                    }

                    ssLHS << "]*";
                    lefthandAlphabetStar = LSR(ssLHS.str());

                    LSR::maxStates = std::max(LSR::maxStates,lefthandAlphabetStar.size());
                }

//                std::cout << "lefthand automaton: " << lefthandAlphabetStar.string() << std::endl;

                LSR righthandAlphabetStar = LSR::one();
                if(righthandReachableLetters.size() != 0) {
                    std::stringstream ssRHS;
                    ssRHS << "[";

                    for(auto letter: righthandReachableLetters) {
                        if(isalnum(letter)) {
                            ssRHS << letter;
                        }
                    }

                    ssRHS << "]*";
                    righthandAlphabetStar = LSR(ssRHS.str());
                    LSR::maxStates = std::max(LSR::maxStates,righthandAlphabetStar.size());
                }

//                std::cout << "righthand automaton: " << righthandAlphabetStar.string() << std::endl;

                LSR closure = lefthandAlphabetStar * middle * righthandAlphabetStar;
                LSR::maxStates = std::max(LSR::maxStates,closure.size());

//                std::cout << "component closure: " << closure.string() << std::endl;
                componentToClosure[i] = closure;
            }
        }
    }

    static LSR findClosureOfLowerLinearTerms(std::set<NonCommutativeMonomial<LSR>> &lowerLinearTerms,
            std::map<VarId, int> &varToComponent,
            std::map<int, LSR> &componentToClosure) {
        LSR closure = LSR::null();

        for(auto &monomial: lowerLinearTerms) {

//            std::string leadingFactor = monomial.getLeadingSR().string();
//            std::string trailingFactor = monomial.getTrailingSR().string();
//
//            std::stringstream ssLHS;
//            ssLHS << "[";
//            for(int i = 0; i < leadingFactor.size(); i++) {
//                if(isalnum(leadingFactor[i])) {
//                    ssLHS << leadingFactor[i];
//                }
//            }
//            ssLHS << "]";
//
//            std::stringstream ssRHS;
//            ssRHS << "[";
//            for(int i = 0; i < trailingFactor.size(); i++) {
//                if(isalnum(trailingFactor[i])) {
//                    ssRHS << trailingFactor[i];
//                }
//            }
//            ssRHS << "]";

            LSR variableClosure = componentToClosure[varToComponent[*(monomial.get_variables().begin())]];

//            std::string lhsString = ssLHS.str();
//            std::string rhsString = ssRHS.str();
//
//            LSR lhsLSR;
//            if(lhsString.size() > 2) {
//                lhsLSR = LSR(lhsString);
//            } else {
//                lhsLSR = LSR::one();
//            }
//
//            LSR rhsLSR;
//            if(rhsString.size() > 2) {
//                rhsLSR = LSR(rhsString);
//            } else {
//                rhsLSR = LSR::one();
//            }
            auto lhsLSR = monomial.getLeadingSR().lossify();
            auto rhsLSR = monomial.getTrailingSR().lossify();

            closure += lhsLSR * variableClosure * rhsLSR;

            LSR::maxStates = std::max(LSR::maxStates,closure.size());
        }

        return closure;
    }

    static LSR compareClosures(const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations_1,
            const VarId &S_1, const std::vector<std::pair<VarId, NonCommutativePolynomial<LSR>>> &equations_2,
            const VarId &S_2, const LSR approx_1, const LSR approx_2, int &largestRawGrammar, int &largestCleanGrammar) {
        int who = RUSAGE_SELF;
        struct rusage usage;
        int ret;
        struct timeval startT, endT;
        long mtime, seconds, useconds;
        gettimeofday(&startT, NULL);

        bool A1_subset_A2 = approx_2.contains(approx_1);
        bool A2_subset_A1 = approx_1.contains(approx_2);

        ret = getrusage(who, &usage);
        auto memoryUsage = usage.ru_maxrss;
//        std::cout << memoryUsage << ",";

        if(A1_subset_A2 && A2_subset_A1) {
            return LSR::null();
//            std::cout << memoryUsage << "," << "0,0,0,0" << std::endl;
        } else {
            Timer timer;
            timer.Start();
            int L21raw = 0, L12raw = 0, L21 = 0, L12 = 0;
            LSR L1_intersect_A2c = LSR::null();
            bool L1_intersect_A2c_changed = false;
            LSR L2_intersect_A1c = LSR::null();
            bool L2_intersect_A1c_changed = false;

            if(!A1_subset_A2) {
                VarId startSymbol_1_2;
                auto A2c = approx_1.minus(approx_2);
                auto intersectionGrammar = A2c.intersectionWithCFG(startSymbol_1_2, S_1, equations_1);
                largestRawGrammar = intersectionGrammar.size() > largestRawGrammar ? intersectionGrammar.size() : largestRawGrammar;
                L12raw = intersectionGrammar.size();
                std::queue<VarId> worklist;
                worklist.push(startSymbol_1_2);
                intersectionGrammar = NonCommutativePolynomial<LSR>::cleanSystem(intersectionGrammar, worklist);
                largestCleanGrammar = intersectionGrammar.size() > largestCleanGrammar ? intersectionGrammar.size() : largestCleanGrammar;
                L12 = intersectionGrammar.size();

                L1_intersect_A2c = NonCommutativePolynomial<LSR>::shortestWord
                        (intersectionGrammar, startSymbol_1_2);
                L1_intersect_A2c_changed = true;
            }

            if(!A2_subset_A1) {
                VarId startSymbol_2_1;
                auto A1c = approx_2.minus(approx_1);
                auto intersectionGrammar_2 = A1c.intersectionWithCFG(startSymbol_2_1, S_2, equations_2);
                largestRawGrammar = intersectionGrammar_2.size() > largestRawGrammar ? intersectionGrammar_2.size() : largestRawGrammar;
                L21raw = intersectionGrammar_2.size();
                std::queue<VarId> worklist;
                worklist.push(startSymbol_2_1);
                intersectionGrammar_2 = NonCommutativePolynomial<LSR>::cleanSystem(intersectionGrammar_2, worklist);
                largestCleanGrammar = intersectionGrammar_2.size() > largestCleanGrammar ? intersectionGrammar_2.size() : largestCleanGrammar;
                L21 = intersectionGrammar_2.size();

                L2_intersect_A1c = NonCommutativePolynomial<LSR>::shortestWord
                        (intersectionGrammar_2, startSymbol_2_1);
                L2_intersect_A1c_changed = true;
            }
            timer.Stop();
            ret = getrusage(who, &usage);
            auto memoryUsage = usage.ru_maxrss;
//            std::cout << memoryUsage << ",";
            if(L1_intersect_A2c_changed && L2_intersect_A1c_changed) {
                auto first = L1_intersect_A2c.string();
                auto second = L2_intersect_A1c.string();

                if(first.size() <= second.size()) {
                    return L1_intersect_A2c;
//                    std::cout << timer.GetMilliseconds().count() << "," << L12raw << "," << L12 << "," << first.size() /*<< "\t" << first*/ << std::endl;
                } else {
                    return L2_intersect_A1c;
//                    std::cout << timer.GetMilliseconds().count() << "," << L21raw << "," << L21 << "," << second.size() /*<< "\t" << second */<< std::endl;
                }
            } else {
                if(L1_intersect_A2c_changed) {
                    return L1_intersect_A2c;
//                    std::cout << timer.GetMilliseconds().count() << "," << L12raw << "," << L12 << "," << L1_intersect_A2c.string().size() /*<< "\t" << L1_intersect_A2c.string() */<< std::endl;
                } else {
                    return L2_intersect_A1c;
//                    std::cout << timer.GetMilliseconds().count() << "," << L21raw << "," << L21 << "," << L2_intersect_A1c.string().size() /*<< "\t" << L2_intersect_A1c.string() */<< std::endl;
                }
            }
        }
    }
};

#endif
