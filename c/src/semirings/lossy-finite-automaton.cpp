#include "lossy-finite-automaton.h"

const LossyFiniteAutomaton LossyFiniteAutomaton::EMPTY = LossyFiniteAutomaton{FiniteAutomaton::EMPTY_AUTOMATON};
const LossyFiniteAutomaton LossyFiniteAutomaton::EPSILON = LossyFiniteAutomaton{FiniteAutomaton::EPSILON_AUTOMATON};
const LossyFiniteAutomaton LossyFiniteAutomaton::UNIVERSE = LossyFiniteAutomaton{FiniteAutomaton::UNIVERSE_AUTOMATON};
