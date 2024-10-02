//
// Created by steffi on 08.01.22.
//

#ifndef STORM_MLModuleDecisionTreeHeuristic_H
#define STORM_MLModuleDecisionTreeHeuristic_H

#include "storm/modelchecker/explorationDT/heuristics/MLDecisionTreeHeuristic.h"
#include "storm/modelchecker/explorationDT/BRTDPInformation.h"
#include "storm/modelchecker/explorationDT/heuristics/ModuleHeuristic.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/decision_tree/decision_tree.hpp>
#include <mlpack/methods/decision_tree/gini_gain.hpp>

namespace storm {
namespace modelchecker {
namespace dtstrat {

template<typename StateType, typename ValueType>
class MLModuleDecisionTreeHeuristic : public MLDecisionTreeHeuristic<StateType, ValueType> {
   public:
    MLModuleDecisionTreeHeuristic() {};
    MLModuleDecisionTreeHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo, double const& explorationEpsilon, bool const& useOnlyStateVariablesOfModuleItself, arma::Row<long unsigned int>& labels,
                                  arma::Mat<int>& data, arma::Row<long unsigned int>& modules);
    MLModuleDecisionTreeHeuristic(exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo, double const& explorationEpsilon,
                                bool const& useOnlyStateVariablesOfModuleItself, moduleHeuristicTypes moduleHeuristic, std::vector<uint> moduleOrdering,
                                std::map<uint, uint> mapping, std::vector<const mlpack::tree::DecisionTree<mlpack::tree::GiniGain>*>& treesOfThisRound);

    ~MLModuleDecisionTreeHeuristic();

    MLModuleDecisionTreeHeuristic& operator=(MLModuleDecisionTreeHeuristic<StateType, ValueType> const& copy);

    std::vector<const mlpack::tree::DecisionTree<mlpack::tree::GiniGain>*> getTrees(arma::Row<long unsigned int>& labels,
                                                                                           arma::Mat<int>& data,
                                                                                           arma::Row<long unsigned int>& modules,
                                                                                           uint numberOfModules);

    ModuleHeuristic<StateType>* getModuleHeuristic(arma::Mat<int>& data, arma::Row<long unsigned int>& modules, moduleHeuristicTypes moduleHeuristic, std::vector<uint> ordering={});
    void getModuleToDTMapping();

    unsigned int getDecision(StateValuationDT const& inputValues, StateType const& currentStateId, bool& worked) override;
    std::vector<uint> getAllDecisions(StateValuationDT const& inputValues, StateType const& currentStateId) override;
    void printHeuristic(const std::string &filename) const override;
    bool useHeuristic() override;

   protected:
    exploration_detail_dt::BRTDPInformation<StateType, ValueType>* brtdpInfo;
    std::vector<const mlpack::tree::DecisionTree<mlpack::tree::GiniGain>*> decisionTrees;
    std::map<uint, uint> moduleDTMapping;
    ModuleHeuristic<StateType>* moduleHeuristic;
    bool useOnlyStateVariablesOfModuleItself = true;

   private:
    StateValuationDT reorder(StateValuationDT const& inputValues, uint const& module);
};

}
}
}

#endif  // STORM_MLModuleDecisionTreeHeuristic_H
