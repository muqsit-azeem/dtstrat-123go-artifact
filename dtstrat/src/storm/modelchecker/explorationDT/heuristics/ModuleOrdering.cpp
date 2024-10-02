//
// Created by steffi on 12.01.22.
//

#include <algorithm>
#include <random>

#include "storm/modelchecker/explorationDT/heuristics/ModuleOrdering.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType>
ModuleOrdering<StateType>::ModuleOrdering() {}

template<typename StateType>
ModuleOrdering<StateType>::ModuleOrdering(const uint &numModules, bool const& random) {
    for (uint i=0; i<numModules; i++) moduleOrdering.push_back(i);
    if (random) {
        auto rng = std::default_random_engine {};
        std::shuffle(std::begin(moduleOrdering), std::end(moduleOrdering), rng);
    }
}

template<typename StateType>
ModuleOrdering<StateType>::ModuleOrdering(std::vector<uint> moduleOrdering) : moduleOrdering(moduleOrdering) {
}

template<typename StateType>
void ModuleOrdering<StateType>::printHeuristic(const std::string &filename) const {
    std::ofstream ofs(filename);
    ofs << "Module Ordering Heuristic\n";
    for (auto v : moduleOrdering) ofs << v << " ";
    ofs.close();
}

template<typename StateType>
bool ModuleOrdering<StateType>::useHeuristic() {
    if (!moduleOrdering.empty()) return true;
    else return false;
}

template<typename StateType>
std::vector<uint> ModuleOrdering<StateType>::getAllDecisions(const StateValuationDT &inputValues, const StateType &currentStateId) {
    return moduleOrdering;
}

template<typename StateType>
unsigned int ModuleOrdering<StateType>::getDecision(const StateValuationDT &inputValues, const StateType &currentStateId, bool &worked) {
    if (!moduleOrdering.empty()) {
        worked=true;
        return moduleOrdering[0];
    } else return -1;
}

template<typename StateType>
uint ModuleOrdering<StateType>::getDecision() const {
    if(!moduleOrdering.empty()) return moduleOrdering[0];
    else return 0;
}

template class ModuleOrdering<uint32_t>;

}
}
}