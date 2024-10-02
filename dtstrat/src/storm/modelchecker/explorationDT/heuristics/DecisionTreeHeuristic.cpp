//
// Created by steffi on 05.01.22.
//

#include "storm/modelchecker/explorationDT/heuristics/DecisionTreeHeuristic.h"
#include <random>

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
DecisionTreeHeuristic<StateType, ValueType>::DecisionTreeHeuristic(
    std::string const& file, exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>* explorationInformation, double const& explorationEpsilon,
    bool const& json)
    : explorationInformation(explorationInformation), explorationEpsilon(explorationEpsilon) {
    updateDT(file, json);
}

template<typename StateType, typename ValueType>
bool DecisionTreeHeuristic<StateType, ValueType>::useHeuristic() {
    if (decisionTree) {
        std::random_device rd;   // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0.0, 1.0);
        double decisionEpsilon = dis(gen);
        if (decisionEpsilon > explorationEpsilon)
            return true;
        else
            return false;
    } else return false;
}

template<typename StateType, typename ValueType>
unsigned int DecisionTreeHeuristic<StateType, ValueType>::getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) {
    auto allActions = getAllDecisions(inputValues, currentStateId);
    if (allActions.empty()) return 0;
    else {
        worked = true;
        return allActions[0];
    }
}

template<typename StateType, typename ValueType>
std::vector<uint> DecisionTreeHeuristic<StateType, ValueType>::getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) {
    auto allActions = decisionTree->getDecisions(inputValues);
    //if (!allActions.empty()) STORM_PRINT(" Action: " << allActions[0] << std::endl);
    return explorationInformation->getStormActionFromDTStrat(allActions, currentStateId);
}

template<typename StateType, typename ValueType>
mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* DecisionTreeHeuristic<StateType, ValueType>::getDT(){
    return nullptr;
};


template<typename StateType, typename ValueType>
void DecisionTreeHeuristic<StateType, ValueType>::printHeuristic(std::string const& filename) const {
    std::ofstream ofs(filename + ".dot");
    ofs << decisionTree->getDotString();
    ofs.close();
}
template<typename StateType, typename ValueType>
DecisionTreeHeuristic<StateType, ValueType>::DecisionTreeHeuristic(
    exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>* explorationInformation, const double& explorationEpsilon) : explorationInformation(explorationInformation), explorationEpsilon(explorationEpsilon){
}

template<typename StateType, typename ValueType>
void DecisionTreeHeuristic<StateType, ValueType>::updateDT(const std::string& someString, const bool& fromJson) {
    if (fromJson) {
        nlohmann2::json jsonTree = nlohmann2::json::parse(someString);
        decisionTree = new storage::DecisionTree<uint>(jsonTree);
    } else {
        decisionTree = new storage::DecisionTree<uint>(someString);
    }
}

template<typename StateType, typename ValueType>
DecisionTreeHeuristic<StateType, ValueType>::~DecisionTreeHeuristic() {
    delete decisionTree;
}

template<typename StateType, typename ValueType>
DecisionTreeHeuristic<StateType, ValueType> DecisionTreeHeuristic<StateType, ValueType>::operator=(const DecisionTreeHeuristic<StateType, ValueType>& copy) {
    if (this == &copy) return *this;
    else {
        delete decisionTree;
        decisionTree = copy.decisionTree;
        explorationInformation = copy.explorationInformation;
        explorationEpsilon = copy.explorationEpsilon;
        return *this;
    }
}
template<typename StateType, typename ValueType>
void DecisionTreeHeuristic<StateType, ValueType>::setExplorationEpsilon(double const& explEps) {
    explorationEpsilon = explEps;
}

template class DecisionTreeHeuristic<uint32_t, double>;
}
}
}