//
// Created by steffi on 12.01.22.
//

#include "storm/modelchecker/explorationDT/heuristics/MLDecisionTreeForModules.h"

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
MLDecisionTreeForModules<StateType, ValueType>::MLDecisionTreeForModules(exploration_detail_dt::BRTDPInformation<StateType, ValueType> & brtdpInfo, arma::Mat<int>& data, arma::Row<long unsigned int>& modules) : brtdpInfo(&brtdpInfo) {
    decisionTree = new mlpack::tree::DecisionTree<mlpack::tree::GiniGain>(data, modules, brtdpInfo.getNumberOfModules(), 1);
}

template<typename StateType, typename ValueType>
unsigned int MLDecisionTreeForModules<StateType, ValueType>::getDecision(const StateValuationDT &inputValues, const StateType &currentStateId, bool &worked) {
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
std::vector<uint> MLDecisionTreeForModules<StateType, ValueType>::getAllDecisions(const StateValuationDT &inputValues, const StateType &currentStateId) {
    auto matrix = inputValues.transform();
    arma::Row<long unsigned int> predictions;
    arma::mat probabilities;
    decisionTree->Classify(matrix, predictions, probabilities);
    std::vector<double> probs(probabilities.n_rows);
    for(uint i=0; i<probs.size(); i++) {
        probs[i] = probabilities(i);
    }
    std::vector<uint> result(probs.size());
    std::iota(result.begin(), result.end(),0);
    sort(result.begin(), result.end(), [&](uint i, uint j){return probs[i]>probs[j];});
    return result;
}

template<typename StateType, typename ValueType>
bool MLDecisionTreeForModules<StateType, ValueType>::useHeuristic() {
    if (decisionTree) return true;
    else return false;
}

template<typename StateType, typename ValueType>
void MLDecisionTreeForModules<StateType, ValueType>::printHeuristic(const std::string &filename) const {
    std::ofstream ofs(filename + ".dot");
    ofs << decisionTree->getDot();
    ofs.close();
}

template<typename StateType, typename ValueType>
MLDecisionTreeForModules<StateType, ValueType>::~MLDecisionTreeForModules() {
    delete decisionTree;
}

template<typename StateType, typename ValueType>
MLDecisionTreeForModules<StateType, ValueType>& MLDecisionTreeForModules<StateType, ValueType>::operator=(const MLDecisionTreeForModules<StateType, ValueType> &copy) {
    if (this == &copy) return *this;
    else {
        delete decisionTree;
        decisionTree = copy.decisionTree;
        brtdpInfo = copy.brtdpInfo;
        return *this;
    }
}

template class MLDecisionTreeForModules<uint32_t, double>;
}
}
}