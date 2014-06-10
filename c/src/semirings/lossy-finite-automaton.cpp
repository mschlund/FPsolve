#include "lossy-finite-automaton.h"

int LossyFiniteAutomaton::maxStates = 0;

LossyFiniteAutomaton LossyFiniteAutomaton::EMPTY = LossyFiniteAutomaton(FiniteAutomaton());
LossyFiniteAutomaton LossyFiniteAutomaton::EPSILON = LossyFiniteAutomaton(FiniteAutomaton::epsilon());
