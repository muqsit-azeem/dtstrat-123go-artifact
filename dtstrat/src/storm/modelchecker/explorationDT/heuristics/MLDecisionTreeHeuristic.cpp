//
// Created by steffi on 05.01.22.
//

#include "storm/modelchecker/explorationDT/heuristics/MLDecisionTreeHeuristic.h"
#include <random>

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
MLDecisionTreeHeuristic<StateType, ValueType>::MLDecisionTreeHeuristic(
    arma::mat& data, arma::Row<long unsigned int>& labels, uint numberOfActions,
    exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>& explorationInformation, double const& explorationEpsilon)
    : explorationInformation(&explorationInformation), explorationEpsilon(explorationEpsilon) {
    updateDT(data, labels, numberOfActions);
}

template<typename StateType, typename ValueType>
bool MLDecisionTreeHeuristic<StateType, ValueType>::useHeuristic() {
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
unsigned int MLDecisionTreeHeuristic<StateType, ValueType>::getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) {
    auto allActions = getAllDecisions(inputValues, currentStateId);
    if (allActions.empty()) return 0;
    else {
        worked = true;
        return allActions[0];
    }
}

template<typename StateType, typename ValueType>
std::vector<uint> MLDecisionTreeHeuristic<StateType, ValueType>::getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) {
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
    return explorationInformation->getStormActionFromDTStrat(result, currentStateId);
}
template<typename StateType, typename ValueType>
mlpack::tree::DecisionTree<mlpack::tree::GiniGain>* MLDecisionTreeHeuristic<StateType, ValueType>::getDT() {
    return this->decisionTree;
}

template<typename StateType, typename ValueType>
void MLDecisionTreeHeuristic<StateType, ValueType>::printHeuristic(std::string const& filename) const {
    std::ofstream ofs(filename + ".dot");
    ofs << decisionTree->getDot();
    ofs.close();
}

template<typename StateType, typename ValueType>
MLDecisionTreeHeuristic<StateType, ValueType>::MLDecisionTreeHeuristic(
    exploration_detail_dt::ExplorationInformationDT<StateType, ValueType>* explorationInformation, const double& explorationEpsilon) : explorationInformation(explorationInformation), explorationEpsilon(explorationEpsilon){
}

template<typename StateType, typename ValueType>
void MLDecisionTreeHeuristic<StateType, ValueType>::updateDT(arma::mat& data, arma::Row<long unsigned int>& labels, uint numberOfActions) {
    decisionTree =  new mlpack::tree::DecisionTree<mlpack::tree::GiniGain>(data, labels, numberOfActions, 1);
    // for debugging purposes - check if the prediction is the same as the true value
    /*for (uint i=0; i<data.n_cols; i++) {
        auto v = data.col(i);
        auto predicted = decisionTree->Classify(v);
        STORM_PRINT("true: " << labels(i) << " predicted: " << predicted << std::endl);
    }*/
}

template<typename StateType, typename ValueType>
MLDecisionTreeHeuristic<StateType, ValueType>::~MLDecisionTreeHeuristic() {
    delete decisionTree;
}

template<typename StateType, typename ValueType>
MLDecisionTreeHeuristic<StateType, ValueType> MLDecisionTreeHeuristic<StateType, ValueType>::operator=(const MLDecisionTreeHeuristic<StateType, ValueType>& copy) {
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
void MLDecisionTreeHeuristic<StateType, ValueType>::setExplorationEpsilon(double const& explEps) {
    explorationEpsilon = explEps;
}

template class MLDecisionTreeHeuristic<uint32_t, double>;
}
}
}