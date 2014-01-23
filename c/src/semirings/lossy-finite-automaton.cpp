#include "lossy-finite-automaton.h"

LossyFiniteAutomaton LossyFiniteAutomaton::EMPTY = LossyFiniteAutomaton(FiniteAutomaton::EMPTY_AUTOMATON);
LossyFiniteAutomaton LossyFiniteAutomaton::EPSILON = LossyFiniteAutomaton(FiniteAutomaton::EPSILON_AUTOMATON);
LossyFiniteAutomaton LossyFiniteAutomaton::UNIVERSE = LossyFiniteAutomaton(FiniteAutomaton::UNIVERSE_AUTOMATON);
