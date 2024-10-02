//
// Created by steffi on 12.01.22.
//

#include "storm/modelchecker/explorationDT/heuristics/DecisionTreeForModules.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
DecisionTreeForModules<StateType, ValueType>::DecisionTreeForModules(exploration_detail_dt::BRTDPInformation<StateType, ValueType> & brtdpInfo, const std::string &fileName) : brtdpInfo(&brtdpInfo) {
    decisionTree = new storm::storage::DecisionTree<uint>(fileName);
}

template<typename StateType, typename ValueType>
unsigned int DecisionTreeForModules<StateType, ValueType>::getDecision(const StateValuationDT &inputValues, const StateType &currentStateId, bool &worked) {
    auto allModules = getAllDecisions(inputValues, currentStateId);
    if (allModules.empty()) {
        worked = false;
        return -1;
    } else {
        worked = true;
        return allModules[0];
    }
}

template<typename StateType, typename ValueType>
std::vector<uint> DecisionTreeForModules<StateType, ValueType>::getAllDecisions(const StateValuationDT &inputValues, const StateType &currentStateId) {
    auto allModules = decisionTree->getDecisions(inputValues);
    return allModules;
}

template<typename StateType, typename ValueType>
bool DecisionTreeForModules<StateType, ValueType>::useHeuristic() {
    if (decisionTree) return true;
    else return false;
}

template<typename StateType, typename ValueType>
void DecisionTreeForModules<StateType, ValueType>::printHeuristic(const std::string &filename) const {
    std::ofstream ofs(filename + ".dot");
    ofs << decisionTree->getDotString();
    ofs.close();
}

template<typename StateType, typename ValueType>
DecisionTreeForModules<StateType, ValueType>::~DecisionTreeForModules() {
    delete decisionTree;
}

template<typename StateType, typename ValueType>
DecisionTreeForModules<StateType, ValueType>& DecisionTreeForModules<StateType, ValueType>::operator=(const DecisionTreeForModules<StateType, ValueType> &copy) {
    if (this == &copy) return *this;
    else {
        delete decisionTree;
        decisionTree = copy.decisionTree;
        brtdpInfo = copy.brtdpInfo;
        return *this;
    }
}

template class DecisionTreeForModules<uint32_t, double>;
}
}
}