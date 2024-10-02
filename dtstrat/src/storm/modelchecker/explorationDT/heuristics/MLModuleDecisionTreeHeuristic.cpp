//
// Created by steffi on 08.01.22.
//

#include "storm/modelchecker/explorationDT/heuristics/MLModuleDecisionTreeHeuristic.h"
#include "storm/modelchecker/explorationDT/heuristics/MLDecisionTreeForModules.h"
#include "storm/modelchecker/explorationDT/heuristics/ModuleOrdering.h"
#include "storm/settings/modules/HeuristicLearningSettings.h"
#include "storm/settings/SettingsManager.h"
#include "storm/exceptions/NotImplementedException.h"
#include "storm/exceptions/NotSupportedException.h"
#include <random>
#include <utility>

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
MLModuleDecisionTreeHeuristic<StateType, ValueType>::MLModuleDecisionTreeHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo,
                                                                                   double const& explorationEpsilon,
                                                                                   bool const& useOnlyStateVariablesOfModuleItself,
                                                                                   arma::Row<long unsigned int>& labels,
                                                                                   arma::Mat<int>& data,
                                                                                   arma::Row<long unsigned int>& modules)
    : brtdpInfo(brtdpInfo),
      useOnlyStateVariablesOfModuleItself(useOnlyStateVariablesOfModuleItself),
      MLDecisionTreeHeuristic<StateType, ValueType>(&(brtdpInfo->getExplorationInformation()), explorationEpsilon) {
    // get a decision tree for each module
    decisionTrees = getTrees(labels, data, modules, brtdpInfo->getNumberOfModules());

    // get the ordering
    const auto& configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    auto moduleHeuristicType = configSettings.getModuleHeuristicType();
    moduleHeuristic = getModuleHeuristic(data, modules, moduleHeuristicType);

    // get module to dt mapping
    getModuleToDTMapping();
}

template<typename StateType, typename ValueType>
MLModuleDecisionTreeHeuristic<StateType, ValueType>::MLModuleDecisionTreeHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo,
                                                                               const double& explorationEpsilon,
                                                                               const bool& useOnlyStateVariablesOfModuleItself,
                                                                               moduleHeuristicTypes moduleHeuristicType, std::vector<uint> moduleOrdering,
                                                                               std::map<uint, uint> mapping,
                                                                               std::vector<const mlpack::tree::DecisionTree<mlpack::tree::GiniGain>*>& treesOfThisRound)
    : brtdpInfo(brtdpInfo),
      useOnlyStateVariablesOfModuleItself(useOnlyStateVariablesOfModuleItself),
      MLDecisionTreeHeuristic<StateType, ValueType>(&(brtdpInfo->getExplorationInformation()), explorationEpsilon),
      decisionTrees(treesOfThisRound) {
    // get a decision tree for each module
    //moduleHeuristic = getModuleHeuristic(moduleHeuristicType, moduleOrdering);
    moduleDTMapping = std::move(mapping);
}

template<typename StateType, typename ValueType>
void MLModuleDecisionTreeHeuristic<StateType, ValueType>::getModuleToDTMapping() {
    const auto& configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    auto moduleHeuristcType = configSettings.getModuleHeuristicType();
    if (moduleHeuristcType == decisiontree) {
        for (uint i = 0; i < decisionTrees.size(); i++) {
            moduleDTMapping[i] = i;
        }
    } else {
        for (uint i = 0; i < decisionTrees.size(); i++) {
            moduleDTMapping[i] = moduleHeuristic->getDecision();
        }
    }
}

template<typename StateType, typename ValueType>
ModuleHeuristic<StateType>* MLModuleDecisionTreeHeuristic<StateType, ValueType>::getModuleHeuristic(arma::Mat<int>& data,
                                                                                                    arma::Row<long unsigned int>& modules,
                                                                                                    moduleHeuristicTypes moduleHeuristicType, std::vector<uint> ordering) {
    if (moduleHeuristicType == decisiontree) {
        return new MLDecisionTreeForModules<StateType, ValueType>(*brtdpInfo, data, modules);
    } else if (!ordering.empty()) {
        return new ModuleOrdering<StateType>(ordering);
    } else if (moduleHeuristicType == orderingtree) {
        std::vector<std::pair<uint, uint>> indexAndSize;
        for(uint i=0; i<decisionTrees.size(); i++) {
            indexAndSize.push_back({i,decisionTrees[i]->size()});
        }
        std::sort(indexAndSize.begin(), indexAndSize.end(),
                  [](std::pair<uint,double> const &a,
                     std::pair<uint,double> const &b) {
                      return a.second > b.second;
                  });
        std::vector<uint> result;
        for (auto v : indexAndSize) {
            result.push_back(v.first);
        }
        return new ModuleOrdering<StateType>(result);
    } else if (moduleHeuristicType == orderingrandom) {
        return new ModuleOrdering<StateType>(brtdpInfo->getNumberOfModules(), true);
    } else return new ModuleOrdering<StateType>(brtdpInfo->getNumberOfModules(), false);
    return new ModuleOrdering<StateType>(brtdpInfo->getNumberOfModules(), false);
}

template<typename StateType, typename ValueType>
std::vector<const mlpack::tree::DecisionTree<mlpack::tree::GiniGain>*> MLModuleDecisionTreeHeuristic<StateType, ValueType>::getTrees(arma::Row<long unsigned int>& labels,
                                                                                                                                     arma::Mat<int>& data,
                                                                                                                                     arma::Row<long unsigned int>& modules,
                                                                                                                                     uint numberOfModules) {
    std::vector<const mlpack::tree::DecisionTree<mlpack::tree::GiniGain>*> result;
    auto assignment = brtdpInfo->getAssignment();
    for (uint module=0; module<numberOfModules; module++) {
        arma::Mat<int> res = data.cols(arma::find(modules==module));
        arma::Mat<int> res2 = res.rows(arma::find(assignment==module || assignment==-1));
        arma::Row<long unsigned int> labels2 = labels.cols(arma::find(modules==module));
        auto new_tree = new mlpack::tree::DecisionTree<mlpack::tree::GiniGain>(res2, labels2, brtdpInfo->getNumberOfActions(), 1);
        result.push_back(new_tree);
    }
    return result;
}

template<typename StateType, typename ValueType>
StateValuationDT MLModuleDecisionTreeHeuristic<StateType, ValueType>::reorder(StateValuationDT const& inputValues, uint const& module) {
    auto assignment = brtdpInfo->getStateGeneration().getStateVariableModuleAssociation();
    std::vector<int> newValues;
    std::vector<int> temp;
    for (uint i=0; i<inputValues.size(); i++) {
        if (assignment[i]==-1 || assignment[i]==module) {
            newValues.push_back(inputValues[i]);
        } else {
            temp.push_back(inputValues[i]);
        }
    }
    newValues.insert(newValues.end(), temp.begin(), temp.end());
    std::vector<bool> emp1;
    std::vector<storm::RationalNumber> emp2;
    StateValuationDT result = StateValuationDT(newValues);
    return result;
}

template<typename StateType, typename ValueType>
std::vector<uint> MLModuleDecisionTreeHeuristic<StateType, ValueType>::getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) {
    std::vector<uint> result;
    auto allModules = MLDecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getAllModules();
    StateValuationDT correctInput = inputValues;
    std::vector<uint> moduleOrdering = moduleHeuristic->getAllDecisions(inputValues, currentStateId);
    if (moduleOrdering.size()<MLDecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getNumberOfModules()) {
        auto temp = ModuleOrdering<StateType>(MLDecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getNumberOfModules(), true);
        auto more = temp.getAllDecisions(inputValues, currentStateId);
        moduleOrdering.insert(moduleOrdering.end(), more.begin(), more.end());
    }
    for (auto const &module: moduleOrdering) {
        if (useOnlyStateVariablesOfModuleItself) {
            // reorder input values
            correctInput = reorder(inputValues, module);
        }
        auto matrix = correctInput.transform();
        arma::Row<long unsigned int> predictions;
        arma::mat probabilities;
        decisionTrees[moduleDTMapping.at(module)]->Classify(matrix, predictions, probabilities);
        std::vector<double> probs(probabilities.n_rows);
        for(uint i=0; i<probs.size(); i++) {
            probs[i] = probabilities(i);
        }
        std::vector<uint> decisions(probs.size());
        std::iota(decisions.begin(), decisions.end(),0);
        sort(decisions.begin(), decisions.end(), [&](uint i, uint j){return probs[i]>probs[j];});
        std::string moduleName = allModules[module].getName();
        std::vector<uint> possibleActions = MLDecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getDTStratActionFromModuleAction(moduleName, decisions);
        result.insert(result.end(), possibleActions.begin(),
                                 possibleActions.end());
    }
    //if (!result.empty()) STORM_PRINT(result[0] << std::endl);
    return MLDecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getStormActionFromDTStrat(result, currentStateId);
}

template<typename StateType, typename ValueType>
unsigned int MLModuleDecisionTreeHeuristic<StateType, ValueType>::getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) {
    auto decisions = getAllDecisions(inputValues, currentStateId);
    if (decisions.empty()) return 0;
    else {
        worked=true;
        return decisions[0];
    }
}

template<typename StateType, typename ValueType>
void MLModuleDecisionTreeHeuristic<StateType, ValueType>::printHeuristic(const std::string& filename) const {
    for (uint i=0; i<decisionTrees.size(); i++) {
        std::string fileName = brtdpInfo->getModelName() + "_" + std::to_string(i) + ".dot";
        std::ofstream ofs(filename);
        ofs << decisionTrees[i]->getDot();
        ofs.close();
    }
}

template<typename StateType, typename ValueType>
bool MLModuleDecisionTreeHeuristic<StateType, ValueType>::useHeuristic() {
    if (!decisionTrees.empty()) {
        std::random_device rd;   // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0.0, 1.0);
        double decisionEpsilon = dis(gen);
        if (decisionEpsilon > MLDecisionTreeHeuristic<StateType, ValueType>::explorationEpsilon)
            return true;
        else
            return false;
    } else return false;
}

template<typename StateType, typename ValueType>
MLModuleDecisionTreeHeuristic<StateType, ValueType>::~MLModuleDecisionTreeHeuristic() {
    for (auto tree : decisionTrees) delete tree;
    delete moduleHeuristic;
}

template<typename StateType, typename ValueType>
MLModuleDecisionTreeHeuristic<StateType, ValueType>& MLModuleDecisionTreeHeuristic<StateType, ValueType>::operator=(const MLModuleDecisionTreeHeuristic<StateType, ValueType>& copy) {
    if (this == &copy) return *this;
    else {
        for (auto tree : decisionTrees) delete tree;
        delete moduleHeuristic;
        brtdpInfo = copy.brtdpInfo;
        decisionTrees = copy.decisionTrees;
        moduleDTMapping = copy.moduleDTMapping;
        moduleHeuristic = copy.moduleHeuristic;
        useOnlyStateVariablesOfModuleItself = copy.useOnlyStateVariablesOfModuleItself;
        return *this;
    }
}

template class MLModuleDecisionTreeHeuristic<uint32_t, double>;

}
}
}