//
// Created by steffi on 08.01.22.
//

#include "storm/modelchecker/explorationDT/heuristics/ModuleDecisionTreeHeuristic.h"
#include "storm/modelchecker/explorationDT/heuristics/DecisionTreeForModules.h"
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
ModuleDecisionTreeHeuristic<StateType, ValueType>::ModuleDecisionTreeHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo,
                                                                               double const& explorationEpsilon,
                                                                               bool const& useOnlyStateVariablesOfModuleItself)
    : brtdpInfo(brtdpInfo),
      useOnlyStateVariablesOfModuleItself(useOnlyStateVariablesOfModuleItself),
      DecisionTreeHeuristic<StateType, ValueType>(&(brtdpInfo->getExplorationInformation()), explorationEpsilon) {
    // get a decision tree for each module
    decisionTrees = getTrees(this->brtdpInfo->getNumberOfModules(), this->brtdpInfo->getModelName());

    // get the ordering
    const auto& configSettings = storm::settings::getModule<storm::settings::modules::HeuristicLearningSettings>();
    auto moduleHeuristicType = configSettings.getModuleHeuristicType();
    moduleHeuristic = getModuleHeuristic(moduleHeuristicType);

    // get module to dt mapping
    getModuleToDTMapping();
}

template<typename StateType, typename ValueType>
ModuleDecisionTreeHeuristic<StateType, ValueType>::ModuleDecisionTreeHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo,
                                                                               const double& explorationEpsilon,
                                                                               const bool& useOnlyStateVariablesOfModuleItself,
                                                                               moduleHeuristicTypes moduleHeuristicType, std::vector<uint> moduleOrdering,
                                                                               std::map<uint, uint> mapping,
                                                                               std::vector<const storm::storage::DecisionTree<uint>*>& treesOfThisRound)
    : brtdpInfo(brtdpInfo),
      useOnlyStateVariablesOfModuleItself(useOnlyStateVariablesOfModuleItself),
      DecisionTreeHeuristic<StateType, ValueType>(&(brtdpInfo->getExplorationInformation()), explorationEpsilon),
      decisionTrees(treesOfThisRound) {
    // get a decision tree for each module
    moduleHeuristic = getModuleHeuristic(moduleHeuristicType, moduleOrdering);
    moduleDTMapping = std::move(mapping);
}

template<typename StateType, typename ValueType>
void ModuleDecisionTreeHeuristic<StateType, ValueType>::getModuleToDTMapping() {
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
ModuleHeuristic<StateType>* ModuleDecisionTreeHeuristic<StateType, ValueType>::getModuleHeuristic(moduleHeuristicTypes moduleHeuristicType, std::vector<uint> ordering) {
    if (moduleHeuristicType == decisiontree) {
        std::string filename = brtdpInfo->getModelName() + "_moduleDT.dot";
        return new DecisionTreeForModules<StateType, ValueType>(*brtdpInfo, filename);
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
}

template<typename StateType, typename ValueType>
std::vector<const storm::storage::DecisionTree<uint>*> ModuleDecisionTreeHeuristic<StateType, ValueType>::getTrees(uint const& numberOfModules, std::string const& modelName) {
    std::vector<const storm::storage::DecisionTree<uint>*> result;
    for (uint i = 0; i < numberOfModules; i++) {
        try {
            std::string filename = modelName + "_moduleDT_" + std::to_string(i) + ".dot";
            std::ifstream f(filename.c_str());
            if (f.good()) {
                auto* first = new storm::storage::DecisionTree<uint>(filename);
                result.push_back(first);
            } else {
                auto* first = new storm::storage::DecisionTree<uint>;
                result.push_back(first);
            }
        } catch (const std::exception& e) {
            STORM_LOG_THROW(true, storm::exceptions::NotSupportedException, "Problem!!!");
        }
    }
    return result;
}

template<typename StateType, typename ValueType>
StateValuationDT ModuleDecisionTreeHeuristic<StateType, ValueType>::reorder(StateValuationDT const& inputValues, uint const& module) {
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
std::vector<uint> ModuleDecisionTreeHeuristic<StateType, ValueType>::getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) {
    std::vector<uint> result;
    auto allModules = DecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getAllModules();
    StateValuationDT correctInput = inputValues;
    std::vector<uint> moduleOrdering = moduleHeuristic->getAllDecisions(inputValues, currentStateId);
    if (moduleOrdering.size()<DecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getNumberOfModules()) {
        auto temp = ModuleOrdering<StateType>(DecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getNumberOfModules(), true);
        auto more = temp.getAllDecisions(inputValues, currentStateId);
        moduleOrdering.insert(moduleOrdering.end(), more.begin(), more.end());
    }
    for (auto const &module: moduleOrdering) {
        if (useOnlyStateVariablesOfModuleItself) {
            // reorder input values
            correctInput = reorder(inputValues, module);
        }
        std::vector<uint32_t> decisions = decisionTrees[moduleDTMapping.at(module)]->getDecisions(correctInput);
        std::string moduleName = allModules[module].getName();
        std::vector<uint> possibleActions = DecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getDTStratActionFromModuleAction(moduleName, decisions);
        result.insert(result.end(), possibleActions.begin(),
                                 possibleActions.end());
    }
    //if (!result.empty()) STORM_PRINT(result[0] << std::endl);
    return DecisionTreeHeuristic<StateType, ValueType>::explorationInformation->getStormActionFromDTStrat(result, currentStateId);
}

template<typename StateType, typename ValueType>
unsigned int ModuleDecisionTreeHeuristic<StateType, ValueType>::getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) {
    auto decisions = getAllDecisions(inputValues, currentStateId);
    if (decisions.empty()) return 0;
    else {
        worked=true;
        return decisions[0];
    }
}

template<typename StateType, typename ValueType>
void ModuleDecisionTreeHeuristic<StateType, ValueType>::printHeuristic(const std::string& filename) const {
    for (uint i=0; i<decisionTrees.size(); i++) {
        std::string fileName = brtdpInfo->getModelName() + "_" + std::to_string(i) + ".dot";
        std::ofstream ofs(filename);
        ofs << decisionTrees[i]->getDotString();
        ofs.close();
    }
}

template<typename StateType, typename ValueType>
void ModuleDecisionTreeHeuristic<StateType, ValueType>::updateDT(const std::string& someString, const bool& fromJson) {
    STORM_LOG_THROW(false, exceptions::NotImplementedException, "Not yet implemented.");
}

template<typename StateType, typename ValueType>
bool ModuleDecisionTreeHeuristic<StateType, ValueType>::useHeuristic() {
    if (!decisionTrees.empty()) {
        std::random_device rd;   // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0.0, 1.0);
        double decisionEpsilon = dis(gen);
        if (decisionEpsilon > DecisionTreeHeuristic<StateType, ValueType>::explorationEpsilon)
            return true;
        else
            return false;
    } else return false;
}

template<typename StateType, typename ValueType>
ModuleDecisionTreeHeuristic<StateType, ValueType>::~ModuleDecisionTreeHeuristic() {
    for (auto tree : decisionTrees) delete tree;
    delete moduleHeuristic;
}

template<typename StateType, typename ValueType>
ModuleDecisionTreeHeuristic<StateType, ValueType>& ModuleDecisionTreeHeuristic<StateType, ValueType>::operator=(const ModuleDecisionTreeHeuristic<StateType, ValueType>& copy) {
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

template class ModuleDecisionTreeHeuristic<uint32_t, double>;

}
}
}