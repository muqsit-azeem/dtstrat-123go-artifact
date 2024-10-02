//
// Created by steffi on 05.01.22.
//
#include "storm/modelchecker/explorationDT/heuristics/NonHeuristic.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType>
bool NonHeuristic<StateType>::useHeuristic() {
    return false;
}

template<typename StateType>
unsigned int NonHeuristic<StateType>::getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) {
    return -1;
}

template<typename StateType>
std::vector<uint> NonHeuristic<StateType>::getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) {
    return {};
}

template<typename StateType>
mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* NonHeuristic<StateType>::getDT(){
    return nullptr;
};

template<typename StateType>
void NonHeuristic<StateType>::printHeuristic(std::string const& filename) const{
    std::ofstream ofs(filename + ".txt");
    ofs << ".";
    ofs.close();
};

template<typename StateType>
void NonHeuristic<StateType>::setExplorationEpsilon(double const& explEps) {
    // Intentionally left empty
}

template class NonHeuristic<uint32_t>;
}
}
}